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
 * Declares gmx::TopologyFile and gmx::TopologyCompleteness.
 *
 * \author Amr Alhossary <amr_alhossary@hotmail.com>
 * \inlibraryapi
 * \ingroup module_fileio
 */
#ifndef GMX_FILEIO_TOPOLOGYFILE_H
#define GMX_FILEIO_TOPOLOGYFILE_H

#include <string>
#include <vector>

#include "gromacs/fileio/file.h"

namespace gmx
{

//! How much of a molecular topology a file states.
enum class TopologyCompleteness : int
{
    //! No connectivity; it has to be derived, for example from residue templates.
    None,
    //! Some explicit connectivity, such as PDB CONECT records.
    Partial,
    //! A complete description that a simulation can run from, such as a .top or .tpr file.
    Full,
    //! Number of entries.
    Count
};

/*! \libinternal \brief
 * A file that carries molecule definitions and connectivity.
 *
 * A format whose files state only part of a topology, such as PDB with its
 * CONECT records, implements this interface and reports
 * TopologyCompleteness::Partial, so that a caller that needs a complete
 * topology can refuse it before reading.
 */
class TopologyFile : public virtual File
{
public:
    //! The capability flag of this interface.
    static constexpr FileCapability sc_capability = FileCapability::Topology;

    ~TopologyFile() override;

    //! Returns how complete the topology in this format is.
    virtual TopologyCompleteness getCompleteness() const = 0;

    //! Returns the names of the molecule types in the file.
    virtual std::vector<std::string> getMoleculeTypes() const = 0;

    //! Returns the number of molecules of type \p moleculeTypeName.
    virtual int getNumMolecules(const std::string& moleculeTypeName) const = 0;
};

} // namespace gmx

#endif
