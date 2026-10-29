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
 * Declares gmx::StructureFile and gmx::MultimodelStructureFile.
 *
 * \author Amr Alhossary <amr_alhossary@hotmail.com>
 * \inlibraryapi
 * \ingroup module_fileio
 */
#ifndef GMX_FILEIO_STRUCTUREFILE_H
#define GMX_FILEIO_STRUCTUREFILE_H

#include "gromacs/fileio/file.h"
#include "gromacs/utility/arrayref.h"
#include "gromacs/utility/vectypes.h"

struct t_atoms;

namespace gmx
{

/*! \libinternal \brief
 * A file that carries particle positions.
 *
 * Implemented by every format that can describe a set of particles with
 * coordinates, such as GRO, G96 and PDB.
 */
class StructureFile : public virtual File
{
public:
    //! The capability flag of this interface.
    static constexpr FileCapability sc_capability = FileCapability::Structure;

    ~StructureFile() override;

    //! Returns the number of particles in the current model.
    virtual int getNumAtoms() const = 0;

    //! Returns the positions of the current model, in nm.
    virtual ArrayRef<const RVec> getPositions() const = 0;

    /*! \brief
     * Returns the particle and residue information of the current model.
     *
     * Returns nullptr when the format carries none. This is the seam with
     * the existing data model: callers that use t_atoms keep working.
     */
    virtual const t_atoms* getAtoms() const = 0;
};

/*! \libinternal \brief
 * A structure file that can hold more than one model.
 *
 * Examples are PDB files with MODEL records and concatenated GRO files.
 * The StructureFile methods refer to the active model.
 */
class MultimodelStructureFile : public virtual StructureFile
{
public:
    //! The capability flag of this interface.
    static constexpr FileCapability sc_capability = FileCapability::Multimodel;

    ~MultimodelStructureFile() override;

    //! Returns the number of models.
    virtual int getNumModels() const = 0;

    //! Appends a model and returns its index.
    virtual int addModel(ArrayRef<const RVec> positions) = 0;

    //! Makes \p modelIndex the model that StructureFile methods refer to.
    virtual void setActiveModel(int modelIndex) = 0;

    //! Returns the index of the active model.
    virtual int getActiveModel() const = 0;
};

} // namespace gmx

#endif
