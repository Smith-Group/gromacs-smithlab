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
 * Tests ForceProviders
 *
 * \author Amr Alhossary <amr_alhossary@hotmail.com>
 * \ingroup module_mdtypes
 */
#include "gmxpre.h"

#include <vector>

#include <gtest/gtest.h>

#include "gromacs/math/paddedvector.h"
#include "gromacs/mdtypes/enerdata.h"
#include "gromacs/mdtypes/forceoutput.h"
#include "gromacs/mdtypes/iforceprovider.h"
#include "gromacs/utility/arrayref.h"
#include "gromacs/utility/mpicomm.h"
#include "gromacs/utility/vectypes.h"

namespace gmx
{
namespace test
{
namespace
{

//! A force provider that records the output it is given
class RecordingForceProvider final : public IForceProvider
{
public:
    void calculateForces(const ForceProviderInput& /* forceProviderInput */,
                         ForceProviderOutput* forceProviderOutput) override
    {
        forceWithShiftForces_ = forceProviderOutput->forceWithShiftForces_;
        numCalls_++;
    }

    //! The halo force buffer of the last call
    ForceWithShiftForces* forceWithShiftForces_ = nullptr;
    //! The number of calls
    int numCalls_ = 0;
};

TEST(ForceProvidersTest, GivesHaloForcesOnlyToHaloForceProviders)
{
    ForceProviders         forceProviders;
    RecordingForceProvider homeProvider;
    RecordingForceProvider haloProvider;
    forceProviders.addForceProvider(&homeProvider);
    EXPECT_EQ(forceProviders.numHaloForceProviders(), 0);
    forceProviders.addForceProvider(&haloProvider, "", ForceProviderAtoms::HomeAndHalo);
    EXPECT_EQ(forceProviders.numHaloForceProviders(), 1);

    matrix             box = { { 1, 0, 0 }, { 0, 1, 0 }, { 0, 0, 1 } };
    MpiComm            mpiComm(MpiComm::SingleRank{});
    ForceProviderInput input({}, 0, {}, {}, 0.0, 0, box, mpiComm, nullptr);

    PaddedVector<RVec>   force(1);
    ForceWithVirial      forceWithVirial(force.arrayRefWithPadding().unpaddedArrayRef(), false);
    ForceWithShiftForces forceWithShiftForces(force.arrayRefWithPadding(), false, {});
    gmx_enerdata_t       enerd(1, nullptr);
    ForceProviderOutput  output(&forceWithVirial, &enerd);
    output.forceWithShiftForces_ = &forceWithShiftForces;

    forceProviders.calculateForces(input, &output);

    EXPECT_EQ(homeProvider.numCalls_, 1);
    EXPECT_EQ(homeProvider.forceWithShiftForces_, nullptr);
    EXPECT_EQ(haloProvider.numCalls_, 1);
    EXPECT_EQ(haloProvider.forceWithShiftForces_, &forceWithShiftForces);
}

} // namespace
} // namespace test
} // namespace gmx
