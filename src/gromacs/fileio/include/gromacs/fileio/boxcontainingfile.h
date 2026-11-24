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
 * Declares gmx::BoxContainingFile.
 *
 * \author Amr Alhossary <amr_alhossary@hotmail.com>
 * \inlibraryapi
 * \ingroup module_fileio
 */
#ifndef GMX_FILEIO_BOXCONTAININGFILE_H
#define GMX_FILEIO_BOXCONTAININGFILE_H

#include <optional>

#include "gromacs/fileio/file.h"
#include "gromacs/utility/matrix.h"

enum class PbcType : int;

namespace gmx
{

/*! \libinternal \brief
 * A file that can carry a simulation box.
 *
 * Examples are the CRYST1 record of a PDB file and the last line of a GRO
 * file. A format implements this interface when it can store a box; a given
 * file of that format may still have none, for example a PDB file without
 * CRYST1, in which case getBoundaries() returns std::nullopt.
 */
class BoxContainingFile : public virtual File
{
public:
    //! The capability flag of this interface.
    static constexpr FileCapability sc_capability = FileCapability::Box;

    ~BoxContainingFile() override;

    //! Returns the box of the current model, or std::nullopt if the file has none.
    virtual std::optional<Matrix3x3> getBoundaries() const = 0;

    //! Sets the box to be written.
    virtual void setBoundaries(const Matrix3x3& box) = 0;

    //! Returns the periodic boundary type of the stored box.
    virtual PbcType getPbcType() const = 0;

    //! Sets the periodic boundary type to be written with the box.
    virtual void setPbcType(PbcType pbcType) = 0;
};

} // namespace gmx

#endif
