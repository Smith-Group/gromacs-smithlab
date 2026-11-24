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
 * Declares gmx::File, the common base of file objects, and the vocabulary
 * shared by the capability interfaces.
 *
 * A file format is one class that implements the capability interfaces it
 * has, such as StructureFile, BoxContainingFile and TopologyFile. Each of
 * those derives virtually from File, so a format that implements several of
 * them is still a single File.
 *
 * \author Amr Alhossary <amr_alhossary@hotmail.com>
 * \inlibraryapi
 * \ingroup module_fileio
 */
#ifndef GMX_FILEIO_FILE_H
#define GMX_FILEIO_FILE_H

#include <memory>

namespace gmx
{

//! Mode in which a file object is opened.
enum class FileMode : int
{
    //! Read an existing file.
    Read,
    //! Create a file, or overwrite an existing one.
    Write,
    //! Add to the end of an existing file.
    Append,
    //! Number of modes.
    Count
};

/*! \brief
 * What a file format can do, as a set of flags.
 *
 * The first flags name capability interfaces: a file object of a format with
 * FileCapability::Structure implements StructureFile, and so on. The remaining
 * flags are properties of a format that do not need an interface of their own.
 */
enum class FileCapability : unsigned int
{
    //! No capability.
    None = 0,
    //! Particle positions; implemented by StructureFile.
    Structure = 1U << 0U,
    //! Several models of one structure; implemented by MultimodelStructureFile.
    Multimodel = 1U << 1U,
    //! A simulation box; implemented by BoxContainingFile.
    Box = 1U << 2U,
    //! Molecule definitions and connectivity; implemented by TopologyFile.
    Topology = 1U << 3U,
    //! A sequence of frames.
    Trajectory = 1U << 4U,
    //! Velocities can be stored.
    Velocities = 1U << 5U,
};

//! Returns the union of \p a and \p b.
constexpr FileCapability operator|(FileCapability a, FileCapability b)
{
    return static_cast<FileCapability>(static_cast<unsigned int>(a) | static_cast<unsigned int>(b));
}

//! Returns whether \p set contains every flag in \p required.
constexpr bool hasCapabilities(FileCapability set, FileCapability required)
{
    return (static_cast<unsigned int>(set) & static_cast<unsigned int>(required))
           == static_cast<unsigned int>(required);
}

/*! \libinternal \brief
 * Common base of every file object.
 *
 * Capability interfaces derive from this class virtually. Without that, a
 * format implementing two interfaces would contain two File subobjects, and
 * converting a pointer to it into a File pointer would be ambiguous.
 *
 * Creating a file object does no I/O. Whether a format can be opened in a
 * given mode, and which interfaces it implements, is known before read() or
 * write() is called.
 */
class File
{
public:
    virtual ~File();

    //! Reads the file into this object.
    virtual void read() = 0;
    //! Writes the contents of this object to the file.
    virtual void write() = 0;
    //! Returns whether this format can be used in \p mode.
    virtual bool supportsMode(FileMode mode) const = 0;
    //! Returns the name of the format, for use in messages.
    virtual const char* formatName() const = 0;

    //! Returns whether writes are flushed as they are made.
    bool autoFlush() const { return autoFlush_; }
    //! Sets whether writes are flushed as they are made.
    void setAutoFlush(bool autoFlush) { autoFlush_ = autoFlush; }

protected:
    File() = default;

private:
    //! Whether writes are flushed as they are made.
    bool autoFlush_ = false;
};

//! Owning pointer to a file object.
using FilePointer = std::unique_ptr<File>;

} // namespace gmx

#endif
