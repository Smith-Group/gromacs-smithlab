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
 * Declares gmx::IOHelper, the open file that a file object owns.
 *
 * \author Amr Alhossary <amr_alhossary@hotmail.com>
 * \ingroup module_fileio
 */
#ifndef GMX_FILEIO_IOHELPER_H
#define GMX_FILEIO_IOHELPER_H

#include <cstdio>

#include <filesystem>

#include "gromacs/fileio/gmxfio.h"
#include "gromacs/utility/classhelpers.h"

namespace gmx
{

/*! \internal \brief
 * An open file that is closed when the helper goes out of scope.
 *
 * A file object keeps its open file in an IOHelper, so the file is closed
 * when the object is destroyed, also when a reader or writer throws.
 */
class IOHelper
{
public:
    IOHelper() = default;
    //! Opens \p path in \p mode, as gmx_fio_open() does.
    IOHelper(const std::filesystem::path& path, const char* mode) { open(path, mode); }
    ~IOHelper() { close(); }

    //! Opens \p path in \p mode, closing the file that was open.
    void open(const std::filesystem::path& path, const char* mode)
    {
        close();
        fileio_ = gmx_fio_open(path, mode);
    }
    //! Closes the file, if one is open.
    void close()
    {
        if (fileio_ != nullptr)
        {
            gmx_fio_close(fileio_);
            fileio_ = nullptr;
        }
    }
    //! Returns the open file, or nullptr.
    t_fileio* fileio() const { return fileio_; }
    //! Returns the stream of the open file.
    FILE* stream() const { return gmx_fio_getfp(fileio_); }

private:
    //! The open file, or nullptr.
    t_fileio* fileio_ = nullptr;

    GMX_DISALLOW_COPY_MOVE_AND_ASSIGN(IOHelper);
};

} // namespace gmx

#endif
