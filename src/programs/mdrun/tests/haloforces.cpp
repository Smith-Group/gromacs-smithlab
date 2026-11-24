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
 * Tests that a module can write forces on halo atoms
 *
 * A test module adds a pair force between all atoms closer than a
 * short distance, computed on the rank where the first atom is home,
 * with the second atom home or in the halo. The test compares what the
 * module adds to the forces and the virial with the analytical values,
 * so the result is the same for any number of ranks.
 *
 * \author Amr Alhossary <amr_alhossary@hotmail.com>
 * \ingroup module_mdrun_integration_tests
 */
#include "gmxpre.h"

#include "config.h"

#include <array>
#include <atomic>
#include <memory>
#include <mutex>
#include <optional>
#include <string>
#include <tuple>
#include <vector>

#include <gtest/gtest.h>

#include "gromacs/domdec/domdec_struct.h"
#include "gromacs/gpu_utils/capabilities.h"
#include "gromacs/hardware/device_management.h"
#include "gromacs/hardware/hw_info.h"
#include "gromacs/mdrun/mdmodules.h"
#include "gromacs/mdtypes/forceoutput.h"
#include "gromacs/mdtypes/iforceprovider.h"
#include "gromacs/mdtypes/imdmodule.h"
#include "gromacs/pbcutil/ishift.h"
#include "gromacs/pbcutil/pbc.h"
#include "gromacs/trajectory/energyframe.h"
#include "gromacs/trajectory/trajectoryframe.h"
#include "gromacs/utility/arrayref.h"
#include "gromacs/utility/functions.h"
#include "gromacs/utility/message_string_collector.h"
#include "gromacs/utility/stringutil.h"
#include "gromacs/utility/vec.h"
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

//! The distance below which the test module adds a pair force
constexpr real c_pairCutoff = 0.5;
//! The strength of the pair force
constexpr real c_pairForceConstant = 1e4;

/*! \brief Returns the test pair force on the first atom of a pair
 *
 * The force pushes the atoms apart and goes to zero smoothly at
 * c_pairCutoff, so rounding near the cutoff does not change the result.
 */
std::optional<RVec> pairForce(const RVec& dx)
{
    const real distance = norm(dx);
    if (distance >= c_pairCutoff)
    {
        return std::nullopt;
    }
    const real magnitude = c_pairForceConstant * gmx::square(c_pairCutoff - distance);
    return (magnitude / distance) * dx;
}

//! Adds the test pair force between home atoms and home or halo atoms
class HaloPairForceProvider final : public IForceProvider
{
public:
    void calculateForces(const ForceProviderInput& input, ForceProviderOutput* output) override
    {
        numCalls_++;
        ArrayRef<RVec>       force       = output->forceWithShiftForces_->force();
        ArrayRef<RVec>       shiftForces = output->forceWithShiftForces_->shiftForces();
        const bool           virial      = output->forceWithShiftForces_->computeVirial();
        ArrayRef<const RVec> x           = input.x_;

        // With direct halo communication, the local atoms include filler
        // particles, which are not atoms
        const auto isAtom = [&input](int localIndex) {
            return input.dd_ == nullptr || isValidGlobalAtom(input.dd_->globalAtomIndices[localIndex]);
        };

        // PBC only in the dimensions without domain decomposition,
        // the halo has the right images in the others
        t_pbc        pbc;
        const t_pbc* pbcPtr = set_pbc_dd(
                &pbc, PbcType::Xyz, input.dd_ != nullptr ? &input.dd_->numCells : nullptr, true, input.box_);

        for (int i = 0; i < input.homenr_; i++)
        {
            if (!isAtom(i))
            {
                continue;
            }
            for (int j = i + 1; j < input.numAtomsWithHalo_; j++)
            {
                if (!isAtom(j))
                {
                    continue;
                }
                RVec dx;
                int  shift = c_centralShiftIndex;
                if (pbcPtr != nullptr)
                {
                    shift = pbc_dx_aiuc(pbcPtr, x[i], x[j], dx);
                }
                else
                {
                    dx = x[i] - x[j];
                }
                const auto forceOnI = pairForce(dx);
                if (!forceOnI)
                {
                    continue;
                }
                force[i] += *forceOnI;
                force[j] -= *forceOnI;
                if (virial)
                {
                    shiftForces[shift] += *forceOnI;
                    shiftForces[c_centralShiftIndex] -= *forceOnI;
                }
            }
        }
    }

    //! The number of calls on this rank
    std::atomic<int> numCalls_ = 0;
};

//! A module that declares and adds the halo pair force provider
class HaloPairForceModule final : public IMDModule
{
public:
    void describe(MDModuleDescriptor* descriptor) const override
    {
        descriptor->writesHaloForces = true;
    }
    IMdpOptionProvider* mdpOptionProvider() override { return nullptr; }
    IMDOutputProvider*  outputProvider() override { return nullptr; }
    void                initForceProviders(ForceProviders* forceProviders) override
    {
        forceProviders->addForceProvider(&forceProvider_, "", ForceProviderAtoms::HomeAndHalo);
    }
    void subscribeToPreProcessingNotifications(MDModulesNotifiers* /* notifiers */) override {}
    void subscribeToSimulationSetupNotifications(MDModulesNotifiers* /* notifiers */) override {}
    void subscribeToSimulationRunNotifications(MDModulesNotifiers* /* notifiers */) override {}

    //! The provider, which counts its calls
    HaloPairForceProvider forceProvider_;
};

//! The test modules made on the ranks of this process
struct MadeModules
{
    std::mutex                                        mutex;
    std::vector<std::shared_ptr<HaloPairForceModule>> modules;
};

//! Returns the test pair forces on all atoms, computed with full PBC
std::vector<RVec> expectedForces(ArrayRef<const RVec> x, const matrix box, matrix* virial)
{
    t_pbc pbc;
    set_pbc(&pbc, PbcType::Xyz, box);
    std::vector<RVec> force(x.size(), { 0, 0, 0 });
    clear_mat(*virial);
    for (Index i = 0; i < x.ssize(); i++)
    {
        for (Index j = i + 1; j < x.ssize(); j++)
        {
            RVec dx;
            pbc_dx(&pbc, x[i], x[j], dx);
            const auto forceOnI = pairForce(dx);
            if (!forceOnI)
            {
                continue;
            }
            force[i] += *forceOnI;
            force[j] -= *forceOnI;
            for (int d1 = 0; d1 < DIM; d1++)
            {
                for (int d2 = 0; d2 < DIM; d2++)
                {
                    (*virial)[d1][d2] -= 0.5 * dx[d1] * (*forceOnI)[d2];
                }
            }
        }
    }
    return force;
}

//! The forces, positions and box of the first frame of a trajectory
struct FirstFrame
{
    std::vector<RVec> x;
    std::vector<RVec> f;
    matrix            box;
};

//! Reads the first frame of \p trajectoryName
FirstFrame readFirstFrame(const std::string& trajectoryName)
{
    TrajectoryFrameReader reader(trajectoryName);
    EXPECT_TRUE(reader.readNextFrame());
    const auto frame = reader.frame();
    FirstFrame first;
    first.x.assign(frame.x().begin(), frame.x().end());
    first.f.assign(frame.f().begin(), frame.f().end());
    for (int d1 = 0; d1 < DIM; d1++)
    {
        for (int d2 = 0; d2 < DIM; d2++)
        {
            first.box[d1][d2] = frame.box()(d1, d2);
        }
    }
    return first;
}

//! The names of the virial energy terms
const std::array<std::string, 9> c_virialTermNames = { "Vir-XX", "Vir-XY", "Vir-XZ",
                                                       "Vir-YX", "Vir-YY", "Vir-YZ",
                                                       "Vir-ZX", "Vir-ZY", "Vir-ZZ" };

//! Reads the virial of the first frame of \p edrFileName
std::array<real, 9> readFirstVirial(const std::string& edrFileName)
{
    auto reader = openEnergyFileToReadTerms(
            edrFileName, std::vector<std::string>(c_virialTermNames.begin(), c_virialTermNames.end()));
    EXPECT_TRUE(reader->readNextFrame());
    const auto          frame = reader->frame();
    std::array<real, 9> virial;
    for (size_t term = 0; term < c_virialTermNames.size(); term++)
    {
        virial[term] = frame.at(c_virialTermNames[term]);
    }
    return virial;
}

//! Where the nonbonded forces and the update run
using HaloForcesTestParameters = std::tuple<std::string, std::string>;

class HaloForcesTest : public MdrunTestFixture, public ::testing::WithParamInterface<HaloForcesTestParameters>
{
};

TEST_P(HaloForcesTest, ModuleForcesAndVirialMatchAnalyticalValues)
{
    const auto [nonbonded, update] = GetParam();
#if GMX_GPU
    const bool haveCompatibleDevices = !getCompatibleDevices(s_hwinfo->deviceInfoList).empty();
    const bool haveGpuUpdate         = GpuConfigurationCapabilities::UpdateLeapfrog
                               && GpuConfigurationCapabilities::Constraints;
#else
    const bool haveCompatibleDevices = false;
    const bool haveGpuUpdate         = false;
#endif
    if (nonbonded == "gpu" && !haveCompatibleDevices)
    {
        GTEST_SKIP() << "Needs a compatible GPU";
    }
    if (update == "gpu" && !haveGpuUpdate)
    {
        GTEST_SKIP() << "GPU update is not supported by this build";
    }
    if (update == "gpu" && nonbonded == "cpu")
    {
        GTEST_SKIP() << "GPU update needs GPU nonbondeds";
    }

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

    CommandLine mdrunCaller;
    mdrunCaller.addOption("-nb", nonbonded);
    mdrunCaller.addOption("-update", update);

    // Without the module
    const std::string trajectoryWithout = fileManager_.getTemporaryFilePath("without.trr").string();
    const std::string edrWithout        = fileManager_.getTemporaryFilePath("without.edr").string();
    runner_.fullPrecisionTrajectoryFileName_ = trajectoryWithout;
    runner_.edrFileName_                     = edrWithout;
    ASSERT_EQ(0, runner_.callMdrun(mdrunCaller));

    // With the module, made on each rank
    auto madeModules = std::make_shared<MadeModules>();
    runner_.mdModuleFactories_.push_back({ "halo-pair-force",
                                           [madeModules]()
                                           {
                                               auto module = std::make_shared<HaloPairForceModule>();
                                               std::lock_guard<std::mutex> lock(madeModules->mutex);
                                               madeModules->modules.push_back(module);
                                               return module;
                                           } });
    const std::string trajectoryWith = fileManager_.getTemporaryFilePath("with.trr").string();
    const std::string edrWith        = fileManager_.getTemporaryFilePath("with.edr").string();
    runner_.fullPrecisionTrajectoryFileName_ = trajectoryWith;
    runner_.edrFileName_                     = edrWith;
    ASSERT_EQ(0, runner_.callMdrun(mdrunCaller));

    // Every rank calls the module at every step
    const int numRanksInThisProcess = (GMX_THREAD_MPI ? getNumberOfTestMpiRanks() : 1);
    EXPECT_EQ(madeModules->modules.size(), static_cast<size_t>(numRanksInThisProcess));
    for (const auto& module : madeModules->modules)
    {
        EXPECT_EQ(module->forceProvider_.numCalls_, nsteps + 1);
    }

    const FirstFrame without = readFirstFrame(trajectoryWithout);
    const FirstFrame with    = readFirstFrame(trajectoryWith);
    ASSERT_EQ(without.f.size(), with.f.size());

    matrix                  expectedVirial;
    const std::vector<RVec> expected = expectedForces(with.x, with.box, &expectedVirial);

    const FloatingPointTolerance forceTolerance = absoluteTolerance(0.05);
    for (size_t i = 0; i < expected.size(); i++)
    {
        const RVec added = with.f[i] - without.f[i];
        for (int d = 0; d < DIM; d++)
        {
            EXPECT_REAL_EQ_TOL(expected[i][d], added[d], forceTolerance)
                    << formatString("force on atom %zu, dimension %d", i, d);
        }
    }

    const std::array<real, 9> virialWithout = readFirstVirial(edrWithout);
    const std::array<real, 9> virialWith    = readFirstVirial(edrWith);
    for (int d1 = 0; d1 < DIM; d1++)
    {
        for (int d2 = 0; d2 < DIM; d2++)
        {
            const int                    term  = (d1 * DIM) + d2;
            const real                   added = virialWith[term] - virialWithout[term];
            const FloatingPointTolerance virialTolerance =
                    absoluteTolerance((1e-4 * std::abs(virialWith[term])) + 0.05);
            EXPECT_REAL_EQ_TOL(expectedVirial[d1][d2], added, virialTolerance) << c_virialTermNames[term];
        }
    }
}

INSTANTIATE_TEST_SUITE_P(WithCpuAndGpu,
                         HaloForcesTest,
                         ::testing::Combine(::testing::Values("cpu", "gpu"),
                                            ::testing::Values("cpu", "gpu")));

} // namespace
} // namespace test
} // namespace gmx
