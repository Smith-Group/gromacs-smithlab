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
/*! \libinternal \file
 * \brief
 * Declares gmx::TopFile, the file object of a topology (.top) file.
 *
 * \author Amr Alhossary <amr_alhossary@hotmail.com>
 * \inlibraryapi
 * \ingroup module_fileio
 */
#ifndef GMX_FILEIO_TOPFILE_H
#define GMX_FILEIO_TOPFILE_H

#include <filesystem>
#include <string>
#include <utility>
#include <vector>

#include "gromacs/fileio/topologyfile.h"
#include "gromacs/utility/arrayref.h"

namespace gmx
{

/*! \libinternal \brief
 * A topology file, as far as tools edit it.
 *
 * read() keeps the lines of the file; write() writes them back and keeps
 * a backup of the old file. In between, the [ molecules ] and [ system ]
 * sections can be changed. The rest of the file, including \#include and
 * \#define lines, is kept as it is: grompp interprets it.
 */
class TopFile : public TopologyFile
{
public:
    //! Creates the object for \p path; no I/O is done.
    explicit TopFile(std::filesystem::path path);
    ~TopFile() override;

    void        read() override;
    void        write() override;
    bool        supportsMode(FileMode mode) const override { return mode != FileMode::Append; }
    const char* formatName() const override { return "topology"; }

    TopologyCompleteness getCompleteness() const override { return TopologyCompleteness::Full; }
    //! Returns the molecule types in [ molecules ], in the order they first appear.
    std::vector<std::string> getMoleculeTypes() const override;
    //! Returns the number of molecules of \p moleculeTypeName in [ molecules ].
    int getNumMolecules(const std::string& moleculeTypeName) const override;

    //! Returns the name in [ system ], or an empty string if there is none.
    std::string getSystemName() const;
    //! Replaces the name in [ system ] with \p name.
    void setSystemName(const std::string& name);
    //! Adds a line for \p count molecules of \p moleculeTypeName at the end of [ molecules ].
    void addMolecules(const std::string& moleculeTypeName, int count);
    /*! \brief
     * Replaces \p count molecules of \p moleculeTypeName with \p replacements.
     *
     * The molecules are taken from the last line of \p moleculeTypeName in
     * [ molecules ]. The replacements follow that line; those with a count of
     * zero are left out.
     *
     * \throws InvalidInputError if there is no such line, or it has fewer
     *         than \p count molecules.
     */
    void replaceMolecules(const std::string&                          moleculeTypeName,
                          int                                         count,
                          ArrayRef<const std::pair<std::string, int>> replacements);

private:
    //! Returns the range of lines of section \p name: first line after the header, one past the last.
    std::pair<size_t, size_t> section(const char* name) const;
    //! Returns the index of the line in [ system ] that holds the name, or the size of lines_.
    size_t systemNameLine() const;

    //! The path of the file.
    std::filesystem::path path_;
    //! The lines of the file, each with its line break.
    std::vector<std::string> lines_;
};

//! Creates the file object for a topology file, for FileFormatRegistry.
FilePointer createTopFile(const std::filesystem::path& path, FileMode mode);

} // namespace gmx

#endif
