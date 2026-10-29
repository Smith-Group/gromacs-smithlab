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
 * Declares gmx::FileFormatRegistry, which finds file formats by what they can do.
 *
 * Every format is described by a FileFormatInfo: its extension, what it can do
 * (FileCapability), whether it can be read and written, and optionally a
 * factory that creates a file object for it. Callers ask the registry for a
 * capability, never for a format:
 * \code
   std::unique_ptr<StructureFile> file = openFileAs<StructureFile>(path, FileMode::Read);
   file->read();
   \endcode
 * A request the format cannot satisfy is refused with the name of the format
 * before any I/O is done.
 *
 * The formats that GROMACS reads and writes today are registered as
 * descriptions without a factory; formats move to file objects one at a time.
 * A format defined outside this module registers itself from its own
 * translation unit with a FileFormatRegistrar.
 *
 * \author Amr Alhossary <amr_alhossary@hotmail.com>
 * \inlibraryapi
 * \ingroup module_fileio
 */
#ifndef GMX_FILEIO_FILEFORMATREGISTRY_H
#define GMX_FILEIO_FILEFORMATREGISTRY_H

#include <deque>
#include <filesystem>
#include <functional>
#include <memory>
#include <optional>
#include <string>
#include <vector>

#include "gromacs/fileio/file.h"
#include "gromacs/utility/gmxassert.h"

namespace gmx
{

/*! \libinternal \brief
 * Everything the registry knows about one file format.
 */
struct FileFormatInfo
{
    //! Factory type: creates a file object for a path and mode, without doing I/O.
    using Factory = std::function<FilePointer(const std::filesystem::path&, FileMode)>;

    //! Extension including the leading dot, such as ".gro". Matched without regard to case.
    std::string extension;
    //! Description for messages and help output.
    std::string description;
    //! What the format can do.
    FileCapability capabilities = FileCapability::None;
    //! Whether files of this format can be read.
    bool canRead = false;
    //! Whether files of this format can be written.
    bool canWrite = false;
    //! Creates a file object; empty for a format that is not yet available as one.
    Factory factory;
    //! The corresponding GromacsFileType, if there is one.
    std::optional<int> legacyType;
};

/*! \libinternal \brief
 * Registry of the file formats known to GROMACS.
 *
 * The registry is created on first use and holds every format GROMACS
 * reads or writes through the functions in confio.h and trxio.h, plus any
 * format registered later.
 *
 * Formats must be registered before the registry is used from more than one
 * thread; lookups are const and can then be made concurrently.
 */
class FileFormatRegistry
{
public:
    //! Returns the registry of the process.
    static FileFormatRegistry& instance();

    /*! \brief
     * Adds a format.
     *
     * \throws InternalError if a format with the same extension is already registered.
     */
    void registerFormat(FileFormatInfo info);

    //! Returns the format of \p path, chosen by extension, or nullptr if it is unknown.
    const FileFormatInfo* formatOf(const std::filesystem::path& path) const;

    //! Returns the format registered for \p legacyType, or nullptr if there is none.
    const FileFormatInfo* byLegacyType(int legacyType) const;

    //! Returns every format whose capabilities include all of \p required.
    std::vector<const FileFormatInfo*> formatsProviding(FileCapability required) const;

    //! Returns every registered format, in order of registration.
    std::vector<const FileFormatInfo*> allFormats() const;

    /*! \brief
     * Creates a file object for \p path in \p mode, without doing I/O.
     *
     * \throws InvalidInputError if the extension is unknown, or if the
     *         format cannot be opened in \p mode. The message names the format.
     * \throws NotImplementedError if the format is known but not yet
     *         available as a file object.
     */
    FilePointer open(const std::filesystem::path& path, FileMode mode) const;

private:
    FileFormatRegistry();

    //! Registered formats; a deque so that pointers stay valid as formats are added.
    std::deque<FileFormatInfo> formats_;
};

/*! \libinternal \brief
 * Registers a format with FileFormatRegistry::instance() when constructed.
 *
 * Intended for a constant at namespace scope in the translation unit that
 * implements the format:
 * \code
   const FileFormatRegistrar c_xyzRegistrar({ ".xyz", "XYZ coordinate file",
           FileCapability::Structure, true, true,
           [](const std::filesystem::path& path, FileMode mode) -> FilePointer
           { return std::make_unique<XyzFile>(path, mode); } });
   \endcode
 */
class FileFormatRegistrar
{
public:
    //! Registers \p info.
    explicit FileFormatRegistrar(FileFormatInfo info);
};

/*! \brief
 * Creates a file object for \p path that implements \p Interface.
 *
 * \tparam Interface  A capability interface, such as StructureFile.
 * \throws InvalidInputError if the format of \p path does not provide
 *         \p Interface, before the file object is created. Otherwise as
 *         FileFormatRegistry::open().
 */
template<typename Interface>
std::unique_ptr<Interface> openFileAs(const std::filesystem::path& path, FileMode mode);

//! Returns a description of the single capability \p capability, such as "structure file".
const char* capabilityDescription(FileCapability capability);

//! \cond internal
namespace detail
{
//! Throws InvalidInputError unless the format of \p path provides \p capability.
void requireCapability(const std::filesystem::path& path, FileCapability capability);
} // namespace detail
//! \endcond

template<typename Interface>
std::unique_ptr<Interface> openFileAs(const std::filesystem::path& path, FileMode mode)
{
    detail::requireCapability(path, Interface::sc_capability);
    FilePointer file  = FileFormatRegistry::instance().open(path, mode);
    auto*       typed = dynamic_cast<Interface*>(file.get());
    // open() has checked that the object implements every interface its format declares.
    GMX_RELEASE_ASSERT(typed != nullptr,
                       "File object does not implement an interface its format declares");
    file.release();
    return std::unique_ptr<Interface>(typed);
}

} // namespace gmx

#endif
