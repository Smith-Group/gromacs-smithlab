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
 * Defines gmx::FileFormatRegistry and the description of the formats that
 * confio.h and trxio.h read and write.
 *
 * \author Amr Alhossary <amr_alhossary@hotmail.com>
 * \ingroup module_fileio
 */
#include "gmxpre.h"

#include "gromacs/fileio/fileformatregistry.h"

#include <string>
#include <utility>

#include "gromacs/fileio/boxcontainingfile.h"
#include "gromacs/fileio/filetypes.h"
#include "gromacs/fileio/structurefile.h"
#include "gromacs/fileio/topologyfile.h"
#include "gromacs/utility/cstringutil.h"
#include "gromacs/utility/exceptions.h"
#include "gromacs/utility/stringutil.h"

#include "trajectoryfileformats.h"

namespace gmx
{

namespace
{

//! A format that confio.h or trxio.h handles, described by what it can do.
struct LegacyFormat
{
    //! The GromacsFileType of the format.
    int type;
    //! What the format can do.
    FileCapability capabilities;
    //! Whether the format can be read.
    bool canRead;
    //! Whether the format can be written.
    bool canWrite;
    //! Creates the file object, or nullptr if the format has none yet.
    FilePointer (*factory)(const std::filesystem::path&, FileMode) = nullptr;
};

/*! \brief The formats that confio.h and trxio.h read and write.
 *
 * Extensions and descriptions come from filetypes.cpp. What each format can
 * do is stated here; the generic type lists in filetypes.cpp (efTRX, efTRO,
 * efSTX, efSTO) are checked against this table by the tests.
 */
const LegacyFormat c_legacyFormats[] = {
    { efXTC, FileCapability::Trajectory, true, true, createXtcFile },
    { efTRR, FileCapability::Trajectory | FileCapability::Velocities, true, true, createTrrFile },
    { efCPT, FileCapability::Trajectory | FileCapability::Velocities, true, false },
    { efTNG, FileCapability::Trajectory | FileCapability::Velocities, true, true, createTngFile },
    { efH5MD, FileCapability::Trajectory | FileCapability::Velocities, true, true, createH5mdFile },
    { efGRO,
      FileCapability::Structure | FileCapability::Multimodel | FileCapability::Box
              | FileCapability::Trajectory | FileCapability::Velocities,
      true,
      true },
    { efG96,
      FileCapability::Structure | FileCapability::Box | FileCapability::Trajectory | FileCapability::Velocities,
      true,
      true },
    { efPDB,
      FileCapability::Structure | FileCapability::Multimodel | FileCapability::Box
              | FileCapability::Topology | FileCapability::Trajectory,
      true,
      true },
    { efCIF,
      FileCapability::Structure | FileCapability::Multimodel | FileCapability::Box
              | FileCapability::Topology | FileCapability::Trajectory,
      true,
      true },
    { efBRK,
      FileCapability::Structure | FileCapability::Multimodel | FileCapability::Box | FileCapability::Topology,
      true,
      true },
    { efENT,
      FileCapability::Structure | FileCapability::Multimodel | FileCapability::Box | FileCapability::Topology,
      true,
      true },
    { efESP, FileCapability::Structure | FileCapability::Box | FileCapability::Velocities, true, true },
    { efPQR, FileCapability::Structure | FileCapability::Box, false, true },
    { efTPR,
      FileCapability::Structure | FileCapability::Box | FileCapability::Topology | FileCapability::Velocities,
      true,
      false },
    { efTOP, FileCapability::Topology, true, true },
};

//! Returns the name of \p mode for messages.
const char* modeDescription(FileMode mode)
{
    switch (mode)
    {
        case FileMode::Read: return "reading";
        case FileMode::Write: return "writing";
        case FileMode::Append: return "appending";
        default: GMX_THROW(InternalError("Invalid file mode"));
    }
}

//! Returns whether \p file implements the interface for \p capability.
bool implementsInterface(const File& file, FileCapability capability)
{
    switch (capability)
    {
        case FileCapability::Structure: return dynamic_cast<const StructureFile*>(&file) != nullptr;
        case FileCapability::Multimodel:
            return dynamic_cast<const MultimodelStructureFile*>(&file) != nullptr;
        case FileCapability::Box: return dynamic_cast<const BoxContainingFile*>(&file) != nullptr;
        case FileCapability::Topology: return dynamic_cast<const TopologyFile*>(&file) != nullptr;
        default: GMX_THROW(InternalError("Capability has no interface"));
    }
}

//! Throws InternalError unless \p file implements exactly the interfaces that \p format declares.
void checkInterfacesMatch(const File& file, const FileFormatInfo& format)
{
    for (const FileCapability capability :
         { FileCapability::Structure, FileCapability::Multimodel, FileCapability::Box, FileCapability::Topology })
    {
        const bool declared    = hasCapabilities(format.capabilities, capability);
        const bool implemented = implementsInterface(file, capability);
        if (declared != implemented)
        {
            GMX_THROW(InternalError(formatString(
                    "The file object for %s files %s the interface for a %s, but the format %s it",
                    format.extension.c_str(),
                    implemented ? "implements" : "does not implement",
                    capabilityDescription(capability),
                    declared ? "declares" : "does not declare")));
        }
    }
}

} // namespace

const char* capabilityDescription(FileCapability capability)
{
    switch (capability)
    {
        case FileCapability::None: return "file";
        case FileCapability::Structure: return "structure file";
        case FileCapability::Multimodel: return "multi-model structure file";
        case FileCapability::Box: return "file with a simulation box";
        case FileCapability::Topology: return "topology file";
        case FileCapability::Trajectory: return "trajectory file";
        case FileCapability::Velocities: return "file with velocities";
        default: GMX_THROW(InternalError("Not a single capability"));
    }
}

FileFormatRegistry::FileFormatRegistry()
{
    for (const LegacyFormat& format : c_legacyFormats)
    {
        FileFormatInfo info;
        info.extension    = ftp2ext_with_dot(format.type);
        info.description  = ftp2desc(format.type);
        info.capabilities = format.capabilities;
        info.canRead      = format.canRead;
        info.canWrite     = format.canWrite;
        info.legacyType   = format.type;
        if (format.factory != nullptr)
        {
            info.factory = format.factory;
        }
        registerFormat(std::move(info));
    }
}

FileFormatRegistry& FileFormatRegistry::instance()
{
    static FileFormatRegistry registry;
    return registry;
}

void FileFormatRegistry::registerFormat(FileFormatInfo info)
{
    if (formatOf(info.extension) != nullptr)
    {
        GMX_THROW(InternalError(formatString(
                "A file format with extension '%s' is already registered", info.extension.c_str())));
    }
    if (info.legacyType.has_value() && byLegacyType(info.legacyType.value()) != nullptr)
    {
        GMX_THROW(InternalError(formatString("A file format for file type %d is already registered",
                                             info.legacyType.value())));
    }
    formats_.push_back(std::move(info));
}

const FileFormatInfo* FileFormatRegistry::formatOf(const std::filesystem::path& path) const
{
    // As fn2ftp(), a path that is only an extension, such as ".gro", names that extension.
    const std::string extension =
            path.has_extension() ? path.extension().string() : path.filename().string();
    for (const FileFormatInfo& format : formats_)
    {
        if (gmx_strcasecmp(extension.c_str(), format.extension.c_str()) == 0)
        {
            return &format;
        }
    }
    return nullptr;
}

const FileFormatInfo* FileFormatRegistry::byLegacyType(int legacyType) const
{
    for (const FileFormatInfo& format : formats_)
    {
        if (format.legacyType == legacyType)
        {
            return &format;
        }
    }
    return nullptr;
}

std::vector<const FileFormatInfo*> FileFormatRegistry::formatsProviding(FileCapability required) const
{
    std::vector<const FileFormatInfo*> result;
    for (const FileFormatInfo& format : formats_)
    {
        if (hasCapabilities(format.capabilities, required))
        {
            result.push_back(&format);
        }
    }
    return result;
}

std::vector<const FileFormatInfo*> FileFormatRegistry::allFormats() const
{
    return formatsProviding(FileCapability::None);
}

FilePointer FileFormatRegistry::open(const std::filesystem::path& path, FileMode mode) const
{
    const FileFormatInfo* format = formatOf(path);
    if (format == nullptr)
    {
        GMX_THROW(InvalidInputError(formatString(
                "File '%s' does not have the extension of a known file format", path.string().c_str())));
    }
    const bool canUseMode = (mode == FileMode::Read) ? format->canRead : format->canWrite;
    if (!canUseMode)
    {
        GMX_THROW(InvalidInputError(formatString("File '%s' (%s) cannot be opened for %s",
                                                 path.string().c_str(),
                                                 format->description.c_str(),
                                                 modeDescription(mode))));
    }
    if (!format->factory)
    {
        GMX_THROW(NotImplementedError(
                formatString("File '%s' (%s) cannot yet be opened as a file object",
                             path.string().c_str(),
                             format->description.c_str())));
    }
    FilePointer file = format->factory(path, mode);
    if (!file)
    {
        GMX_THROW(InternalError(formatString("The factory for %s files returned no file object",
                                             format->extension.c_str())));
    }
    if (!file->supportsMode(mode))
    {
        GMX_THROW(InvalidInputError(formatString("File '%s' (%s) cannot be opened for %s",
                                                 path.string().c_str(),
                                                 file->formatName(),
                                                 modeDescription(mode))));
    }
    checkInterfacesMatch(*file, *format);
    return file;
}

FileFormatRegistrar::FileFormatRegistrar(FileFormatInfo info)
{
    FileFormatRegistry::instance().registerFormat(std::move(info));
}

namespace detail
{

void requireCapability(const std::filesystem::path& path, FileCapability capability)
{
    const FileFormatInfo* format = FileFormatRegistry::instance().formatOf(path);
    if (format == nullptr)
    {
        GMX_THROW(InvalidInputError(formatString(
                "File '%s' does not have the extension of a known file format", path.string().c_str())));
    }
    if (!hasCapabilities(format->capabilities, capability))
    {
        GMX_THROW(InvalidInputError(formatString("File '%s' (%s) cannot be used as a %s",
                                                 path.string().c_str(),
                                                 format->description.c_str(),
                                                 capabilityDescription(capability))));
    }
}

} // namespace detail

} // namespace gmx
