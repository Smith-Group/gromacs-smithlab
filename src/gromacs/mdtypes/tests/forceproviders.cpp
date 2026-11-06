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

#include <string>
#include <vector>

#include <gmock/gmock.h>
#include <gtest/gtest.h>

#include "gromacs/math/paddedvector.h"
#include "gromacs/mdtypes/enerdata.h"
#include "gromacs/mdtypes/forceoutput.h"
#include "gromacs/mdtypes/iforceprovider.h"
#include "gromacs/utility/arrayref.h"
#include "gromacs/utility/exceptions.h"
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

//! A hook that appends its name to a shared record when it is called
class NamedHook final : public IForceProvider, public IPostEnergyProvider
{
public:
    NamedHook(std::string name, std::vector<std::string>* calls) :
        name_(std::move(name)), calls_(calls)
    {
    }
    void calculateForces(const ForceProviderInput& /* forceProviderInput */,
                         ForceProviderOutput* /* forceProviderOutput */) override
    {
        calls_->push_back(name_);
    }
    void applyAfterPotentialEnergy(const PostEnergyInput& /* input */, ForceProviderOutput* /* output */) override
    {
        calls_->push_back(name_);
    }

private:
    std::string               name_;
    std::vector<std::string>* calls_;
};

//! Input and output for calling the hooks, with no atoms
class HookCall
{
public:
    HookCall() :
        mpiComm_(MpiComm::SingleRank{}),
        input_({}, 0, {}, {}, 0.0, 0, box_, mpiComm_, nullptr),
        forceWithVirial_({}, true),
        enerd_(1, nullptr),
        output_(&forceWithVirial_, &enerd_)
    {
    }

    matrix              box_ = { { 1, 0, 0 }, { 0, 1, 0 }, { 0, 0, 1 } };
    MpiComm             mpiComm_;
    ForceProviderInput  input_;
    ForceWithVirial     forceWithVirial_;
    gmx_enerdata_t      enerd_;
    ForceProviderOutput output_;
};

TEST(ForceProvidersTest, CallsForceProvidersByPhaseAndPriority)
{
    std::vector<std::string> calls;
    NamedHook                biasLate("bias-late", &calls);
    NamedHook                biasDefault("bias-default", &calls);
    NamedHook                biasDefaultSecond("bias-default-second", &calls);
    NamedHook                collectiveVariables("collective-variables", &calls);
    NamedHook                geometry("geometry", &calls);

    ForceProviders forceProviders;
    forceProviders.addForceProvider(&biasLate, HookOrder{ StepPhase::Biases, 10 });
    forceProviders.addForceProvider(&biasDefault);
    forceProviders.addForceProvider(&collectiveVariables, HookOrder{ StepPhase::CollectiveVariables, 0 });
    forceProviders.addForceProvider(&biasDefaultSecond);
    forceProviders.addForceProvider(&geometry, HookOrder{ StepPhase::PrepareGeometry, 5 });

    HookCall call;
    forceProviders.calculateForces(call.input_, &call.output_);

    EXPECT_THAT(calls,
                ::testing::ElementsAre("geometry",
                                       "collective-variables",
                                       "bias-default",
                                       "bias-default-second",
                                       "bias-late"));
}

TEST(ForceProvidersTest, CallsPostEnergyProvidersByPhaseAndPriority)
{
    std::vector<std::string> calls;
    NamedHook                observer("observer", &calls);
    NamedHook                energyFirst("energy-first", &calls);
    NamedHook                energySecond("energy-second", &calls);

    ForceProviders forceProviders;
    EXPECT_FALSE(forceProviders.hasPostEnergyProvider());
    forceProviders.addPostEnergyProvider(&observer, HookOrder{ StepPhase::Observers, -5 });
    forceProviders.addPostEnergyProvider(&energySecond, HookOrder{ StepPhase::EnergyDependent, 1 });
    forceProviders.addPostEnergyProvider(&energyFirst, HookOrder{ StepPhase::EnergyDependent, 0 });
    EXPECT_TRUE(forceProviders.hasPostEnergyProvider());

    HookCall              call;
    matrix                virial = { { 0 } };
    const PostEnergyInput input{ 0, true, 0, virial };
    forceProviders.applyAfterPotentialEnergy(input, &call.output_);

    EXPECT_THAT(calls, ::testing::ElementsAre("energy-first", "energy-second", "observer"));
}

TEST(ForceProvidersTest, RejectsHooksInTheWrongPhase)
{
    std::vector<std::string> calls;
    NamedHook                hook("hook", &calls);
    ForceProviders           forceProviders;
    EXPECT_THROW(forceProviders.addForceProvider(&hook, HookOrder{ StepPhase::EnergyDependent, 0 }),
                 APIError);
    EXPECT_THROW(forceProviders.addPostEnergyProvider(&hook, HookOrder{ StepPhase::Biases, 0 }), APIError);
}

//! A step preparer that asks for what the test sets
class RequestingPreparer final : public IStepPreparer
{
public:
    RequestingPreparer(bool needPotentialEnergy, bool needVirial) :
        needPotentialEnergy_(needPotentialEnergy), needVirial_(needVirial)
    {
    }
    void prepareStep(const StepInfo& stepInfo, StepRequest* request) override
    {
        preparedSteps_.push_back(stepInfo.step);
        requestWasClearOnEntry_ =
                requestWasClearOnEntry_ && !request->needPotentialEnergy && !request->needVirial;
        request->needPotentialEnergy = needPotentialEnergy_;
        request->needVirial          = needVirial_;
    }

    //! The steps this preparer was called for
    std::vector<int64_t> preparedSteps_;
    //! Whether the request was clear at every call
    bool requestWasClearOnEntry_ = true;

private:
    bool needPotentialEnergy_;
    bool needVirial_;
};

TEST(ForceProvidersTest, CombinesTheRequestsOfAllStepPreparers)
{
    RequestingPreparer energy(true, false);
    RequestingPreparer virial(false, true);
    RequestingPreparer nothing(false, false);

    ForceProviders forceProviders;
    forceProviders.addStepPreparer(&energy);
    forceProviders.addStepPreparer(&virial);
    forceProviders.addStepPreparer(&nothing);

    StepInfo stepInfo;
    stepInfo.step = 7;
    StepRequest request;
    forceProviders.prepareStep(stepInfo, &request);

    EXPECT_TRUE(request.needPotentialEnergy);
    EXPECT_TRUE(request.needVirial);
    for (const RequestingPreparer* preparer : { &energy, &virial, &nothing })
    {
        EXPECT_THAT(preparer->preparedSteps_, ::testing::ElementsAre(7));
        EXPECT_TRUE(preparer->requestWasClearOnEntry_);
    }
}

} // namespace
} // namespace test
} // namespace gmx
