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
 * Tests the hooks before the step and after the potential energy
 *
 * A test module asks for the potential energy at one step where mdrun
 * would not compute it, records what its hooks get, and adds a force
 * after the potential energy is known.
 *
 * \author Amr Alhossary <amr_alhossary@hotmail.com>
 * \ingroup module_mdrun_integration_tests
 */
#include "gmxpre.h"

#include "config.h"

#include <memory>
#include <mutex>
#include <string>
#include <vector>

#include <gmock/gmock.h>
#include <gtest/gtest.h>

#include "gromacs/mdrun/mdmodules.h"
#include "gromacs/mdtypes/forceoutput.h"
#include "gromacs/mdtypes/iforceprovider.h"
#include "gromacs/mdtypes/imdmodule.h"
#include "gromacs/trajectory/energyframe.h"
#include "gromacs/trajectory/trajectoryframe.h"
#include "gromacs/utility/arrayref.h"
#include "gromacs/utility/vectypes.h"

#include "testutils/mpitest.h"
#include "testutils/simulationdatabase.h"
#include "testutils/testasserts.h"
#include "testutils/trajectoryreader.h"

#include "energyreader.h"
#include "moduletest.h"

namespace gmx
{
namespace test
{
namespace
{

//! The step at which the module asks for the potential energy
constexpr int64_t c_energyStep = 3;
//! The force the module adds to each home atom after the potential energy
const RVec c_addedForce = { 0, 0, 1 };

//! What the post-energy hook got at one call
struct PostEnergyCall
{
    int64_t step;
    bool    energyWasComputed;
    real    potentialEnergy;
};

//! A module that asks for the potential energy at c_energyStep and records its hooks
class StepHooksModule final : public IMDModule, public IStepPreparer, public IPostEnergyProvider
{
public:
    void describe(MDModuleDescriptor* descriptor) const override
    {
        descriptor->addsForcesAfterPotentialEnergy = true;
    }
    IMdpOptionProvider* mdpOptionProvider() override { return nullptr; }
    IMDOutputProvider*  outputProvider() override { return nullptr; }
    void                initForceProviders(ForceProviders* forceProviders) override
    {
        forceProviders->addStepPreparer(this);
        forceProviders->addPostEnergyProvider(this, HookOrder{ StepPhase::EnergyDependent, 0 });
    }
    void subscribeToPreProcessingNotifications(MDModulesNotifiers* /* notifiers */) override {}
    void subscribeToSimulationSetupNotifications(MDModulesNotifiers* /* notifiers */) override {}
    void subscribeToSimulationRunNotifications(MDModulesNotifiers* /* notifiers */) override {}

    void prepareStep(const StepInfo& stepInfo, StepRequest* request) override
    {
        preparedSteps_.push_back(stepInfo.step);
        if (stepInfo.isLastStep)
        {
            lastSteps_.push_back(stepInfo.step);
        }
        request->needPotentialEnergy = (stepInfo.step == c_energyStep);
    }

    void applyAfterPotentialEnergy(const PostEnergyInput& input, ForceProviderOutput* output) override
    {
        calls_.push_back({ input.step_, input.energyWasComputed_, input.potentialEnergy_ });
        for (RVec& force : output->forceWithVirial_.force_)
        {
            force += c_addedForce;
        }
    }

    //! The steps prepareStep() was called for
    std::vector<int64_t> preparedSteps_;
    //! The steps prepareStep() was told are the last
    std::vector<int64_t> lastSteps_;
    //! The calls of the post-energy hook
    std::vector<PostEnergyCall> calls_;
};

//! The test modules made on the ranks of this process
struct MadeModules
{
    std::mutex                                    mutex;
    std::vector<std::shared_ptr<StepHooksModule>> modules;
};

//! Returns the forces of the first frame of \p trajectoryName
std::vector<RVec> readFirstForces(const std::string& trajectoryName)
{
    TrajectoryFrameReader reader(trajectoryName);
    EXPECT_TRUE(reader.readNextFrame());
    const auto frame = reader.frame();
    return { frame.f().begin(), frame.f().end() };
}

//! Returns the potential energy of the first frame of \p edrFileName
real readFirstPotentialEnergy(const std::string& edrFileName)
{
    auto reader = openEnergyFileToReadTerms(edrFileName, { "Potential" });
    EXPECT_TRUE(reader->readNextFrame());
    return reader->frame().at("Potential");
}

using StepHooksTest = MdrunTestFixture;

TEST_F(StepHooksTest, ModuleGetsRequestedEnergyAndAddsForces)
{
    const std::string simulationName = "argon5832";
    const int         nsteps         = 4;
    auto              mdpFieldValues = prepareMdpFieldValues(simulationName, "md", "no", "no");
    mdpFieldValues["nsteps"]         = std::to_string(nsteps);
    mdpFieldValues["nstxout"]        = std::to_string(nsteps);
    mdpFieldValues["nstfout"]        = std::to_string(nsteps);
    mdpFieldValues["nstvout"]        = "0";
    mdpFieldValues["nstenergy"]      = std::to_string(nsteps);
    mdpFieldValues["nstcalcenergy"]  = std::to_string(nsteps);
    runner_.useTopGroAndNdxFromDatabase(simulationName);
    runner_.useStringAsMdpFile(prepareMdpFileContents(mdpFieldValues));
    ASSERT_EQ(0, runner_.callGrompp());

    // Without the module
    const std::string trajectoryWithout = fileManager_.getTemporaryFilePath("without.trr").string();
    runner_.fullPrecisionTrajectoryFileName_ = trajectoryWithout;
    ASSERT_EQ(0, runner_.callMdrun());

    // With the module, made on each rank
    auto madeModules = std::make_shared<MadeModules>();
    runner_.mdModuleFactories_.push_back({ "step-hooks",
                                           [madeModules]()
                                           {
                                               auto module = std::make_shared<StepHooksModule>();
                                               std::lock_guard<std::mutex> lock(madeModules->mutex);
                                               madeModules->modules.push_back(module);
                                               return module;
                                           } });
    const std::string trajectoryWith = fileManager_.getTemporaryFilePath("with.trr").string();
    const std::string edrWith        = fileManager_.getTemporaryFilePath("with.edr").string();
    runner_.fullPrecisionTrajectoryFileName_ = trajectoryWith;
    runner_.edrFileName_                     = edrWith;
    ASSERT_EQ(0, runner_.callMdrun());

    const int numRanksInThisProcess = (GMX_THREAD_MPI ? getNumberOfTestMpiRanks() : 1);
    ASSERT_EQ(madeModules->modules.size(), static_cast<size_t>(numRanksInThisProcess));
    real potentialEnergyAtFirstStep = 0;
    for (const auto& module : madeModules->modules)
    {
        EXPECT_THAT(module->preparedSteps_, ::testing::ElementsAre(0, 1, 2, 3, 4));
        EXPECT_THAT(module->lastSteps_, ::testing::ElementsAre(nsteps));
        ASSERT_EQ(module->calls_.size(), static_cast<size_t>(nsteps + 1));
        for (const PostEnergyCall& call : module->calls_)
        {
            // The energy is computed at the energy steps 0 and 4, and where the module asked
            const bool expectEnergy = (call.step == 0 || call.step == nsteps || call.step == c_energyStep);
            EXPECT_EQ(call.energyWasComputed, expectEnergy) << "at step " << call.step;
        }
        potentialEnergyAtFirstStep += module->calls_[0].potentialEnergy;
    }
    // The hook gets the potential energy of its rank
    if (GMX_THREAD_MPI)
    {
        EXPECT_REAL_EQ_TOL(readFirstPotentialEnergy(edrWith),
                           potentialEnergyAtFirstStep,
                           relativeToleranceAsFloatingPoint(potentialEnergyAtFirstStep, 1e-5));
    }

    // Each atom gets the force the module added once, on its home rank
    const std::vector<RVec> without = readFirstForces(trajectoryWithout);
    const std::vector<RVec> with    = readFirstForces(trajectoryWith);
    ASSERT_EQ(without.size(), with.size());
    const FloatingPointTolerance tolerance = absoluteTolerance(1e-3);
    for (size_t i = 0; i < with.size(); i++)
    {
        const RVec added = with[i] - without[i];
        for (int d = 0; d < DIM; d++)
        {
            EXPECT_REAL_EQ_TOL(c_addedForce[d], added[d], tolerance) << "atom " << i << " dimension " << d;
        }
    }
}

} // namespace
} // namespace test
} // namespace gmx
