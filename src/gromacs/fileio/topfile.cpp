/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright 2026- The GROMACS Authors
 * and the project initiators Erik Lindahl, Berk Hess and David van der Spoel.
 * Consult the AUTHORS/COPYING files and https://www.gromacs.org for details.
 *
 * GROMACS is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public License
 * as published by the Free Software Foundation; either version 2.1
 * of the License, or (at your option) any later version.
 *
 * GROMACS is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with GROMACS; if not, see
 * https://www.gnu.org/licenses, or write to the Free Software Foundation,
 * Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA.
 *
 * If you want to redistribute modifications to GROMACS, please
 * consider that scientific software is very special. Version
 * control is crucial - bugs must be traceable. We will be happy to
 * consider code for inclusion in the official distribution, but
 * derived work must not be called official GROMACS. Details are found
 * in the README & COPYING files - if they are missing, get the
 * official version at https://www.gromacs.org.
 *
 * To help us fund GROMACS development, we humbly ask that you cite
 * the research papers on the package. Check out https://www.gromacs.org.
 */
/*! \internal \file
 * \brief
 * Defines gmx::TopFile.
 *
 * \author Amr Alhossary <amr_alhossary@hotmail.com>
 * \ingroup module_fileio
 */
#include "gmxpre.h"

#include "gromacs/fileio/topfile.h"

#include <cstdio>
#include <cstring>

#include <memory>
#include <optional>

#include "gromacs/utility/cstringutil.h"
#include "gromacs/utility/exceptions.h"
#include "gromacs/utility/futil.h"
#include "gromacs/utility/stringutil.h"
#include "gromacs/utility/textreader.h"

namespace gmx
{

namespace
{

//! Returns the name of the section that \p line starts, if it is a section header.
std::optional<std::string> sectionName(const std::string& line)
{
    const std::string text = stripString(line);
    if (text.size() < 2 || text.front() != '[' || text.back() != ']')
    {
        return std::nullopt;
    }
    return stripString(text.substr(1, text.size() - 2));
}

//! Returns whether \p line holds data, rather than being empty or a comment.
bool isDataLine(const std::string& line)
{
    const std::string text = stripString(line);
    return !text.empty() && text.front() != ';';
}

//! Returns the molecule type that a [ molecules ] line names.
std::string moleculeType(const std::string& line)
{
    return splitString(line).front();
}

//! Returns the number of molecules that a [ molecules ] line gives.
int moleculeCount(const std::string& line)
{
    int count = 0;
    std::sscanf(line.c_str(), "%*s %d", &count);
    return count;
}

//! Returns the [ molecules ] line for \p count molecules of \p name.
std::string moleculeLine(const std::string& name, int count)
{
    return formatString("%-15s %5d\n", name.c_str(), count);
}

} // namespace

TopFile::TopFile(std::filesystem::path path) : path_(std::move(path)) {}

TopFile::~TopFile() = default;

void TopFile::read()
{
    const std::string text = TextReader::readFileToString(path_);
    lines_.clear();
    size_t start = 0;
    while (start < text.size())
    {
        const size_t end  = text.find('\n', start);
        const size_t next = (end == std::string::npos) ? text.size() : end + 1;
        lines_.push_back(text.substr(start, next - start));
        start = next;
    }
}

void TopFile::write()
{
    char temporaryName[STRLEN];
    std::strncpy(temporaryName, "temp.topXXXXXX", STRLEN);
    FILE* out = gmx_fopen_temporary(temporaryName);
    for (const std::string& line : lines_)
    {
        std::fputs(line.c_str(), out);
    }
    gmx_ffclose(out);
    make_backup(path_);
    gmx_file_rename(temporaryName, path_);
}

std::pair<size_t, size_t> TopFile::section(const char* name) const
{
    for (size_t header = 0; header < lines_.size(); header++)
    {
        const std::optional<std::string> found = sectionName(lines_[header]);
        if (found && gmx_strcasecmp(found->c_str(), name) == 0)
        {
            size_t end = header + 1;
            while (end < lines_.size() && !sectionName(lines_[end]))
            {
                end++;
            }
            return { header + 1, end };
        }
    }
    return { lines_.size(), lines_.size() };
}

std::vector<std::string> TopFile::getMoleculeTypes() const
{
    std::vector<std::string> types;
    const auto [begin, end] = section("molecules");
    for (size_t i = begin; i < end; i++)
    {
        if (!isDataLine(lines_[i]))
        {
            continue;
        }
        const std::string type  = moleculeType(lines_[i]);
        bool              known = false;
        for (const std::string& knownType : types)
        {
            known = known || gmx_strcasecmp(knownType.c_str(), type.c_str()) == 0;
        }
        if (!known)
        {
            types.push_back(type);
        }
    }
    return types;
}

int TopFile::getNumMolecules(const std::string& moleculeTypeName) const
{
    int count               = 0;
    const auto [begin, end] = section("molecules");
    for (size_t i = begin; i < end; i++)
    {
        if (isDataLine(lines_[i])
            && gmx_strcasecmp(moleculeType(lines_[i]).c_str(), moleculeTypeName.c_str()) == 0)
        {
            count += moleculeCount(lines_[i]);
        }
    }
    return count;
}

size_t TopFile::systemNameLine() const
{
    const auto [begin, end] = section("system");
    for (size_t i = begin; i < end; i++)
    {
        if (isDataLine(lines_[i]))
        {
            return i;
        }
    }
    return lines_.size();
}

std::string TopFile::getSystemName() const
{
    const size_t line = systemNameLine();
    return line < lines_.size() ? stripString(lines_[line]) : std::string();
}

void TopFile::setSystemName(const std::string& name)
{
    const size_t line = systemNameLine();
    if (line == lines_.size())
    {
        GMX_THROW(InvalidInputError(formatString("File '%s' has no name in a [ system ] section",
                                                 path_.string().c_str())));
    }
    lines_[line] = name + "\n";
}

void TopFile::addMolecules(const std::string& moleculeTypeName, int count)
{
    const size_t end = section("molecules").second;
    if (end > 0 && !lines_[end - 1].empty() && lines_[end - 1].back() != '\n')
    {
        lines_[end - 1] += '\n';
    }
    lines_.insert(lines_.begin() + end, moleculeLine(moleculeTypeName, count));
}

void TopFile::replaceMolecules(const std::string&                          moleculeTypeName,
                               int                                         count,
                               ArrayRef<const std::pair<std::string, int>> replacements)
{
    const auto [begin, end] = section("molecules");
    size_t last             = end;
    for (size_t i = begin; i < end; i++)
    {
        if (isDataLine(lines_[i])
            && gmx_strcasecmp(moleculeType(lines_[i]).c_str(), moleculeTypeName.c_str()) == 0)
        {
            last = i;
        }
    }
    if (last == end)
    {
        GMX_THROW(InvalidInputError(formatString(
                "No line with moleculetype '%s' found the [ molecules ] section of file '%s'",
                moleculeTypeName.c_str(),
                path_.string().c_str())));
    }
    const int available = moleculeCount(lines_[last]);
    if (available < count)
    {
        GMX_THROW(InvalidInputError(
                formatString("The last entry for moleculetype '%s' in the [ molecules ] section of "
                             "file '%s' has less solvent molecules (%d) than were replaced (%d)",
                             moleculeTypeName.c_str(),
                             path_.string().c_str(),
                             available,
                             count)));
    }
    std::vector<std::string> newLines;
    if (available > count)
    {
        newLines.push_back(moleculeLine(moleculeTypeName, available - count));
    }
    for (const auto& [name, replacementCount] : replacements)
    {
        if (replacementCount > 0)
        {
            newLines.push_back(moleculeLine(name, replacementCount));
        }
    }
    lines_.erase(lines_.begin() + last);
    lines_.insert(lines_.begin() + last, newLines.begin(), newLines.end());
}

FilePointer createTopFile(const std::filesystem::path& path, FileMode /*mode*/)
{
    return std::make_unique<TopFile>(path);
}

} // namespace gmx
