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
 * Tests that grompp and mdrun load an MDModules plugin named in the mdp file
 *
 * \author Amr Alhossary <amr_alhossary@hotmail.com>
 * \ingroup module_mdrun_integration_tests
 */
#include "gmxpre.h"

#include "config.h"

#include <cstdlib>

#include <string>
#include <vector>

#include <gtest/gtest.h>

#include "gromacs/trajectory/trajectoryframe.h"
#include "gromacs/utility/arrayref.h"
#include "gromacs/utility/vectypes.h"

#include "testutils/simulationdatabase.h"
#include "testutils/testasserts.h"
#include "testutils/trajectoryreader.h"

#include "moduletest.h"

namespace gmx
{
namespace test
{
namespace
{

//! Returns the forces of the first frame of \p trajectoryName
std::vector<RVec> readFirstForces(const std::string& trajectoryName)
{
    TrajectoryFrameReader reader(trajectoryName);
    EXPECT_TRUE(reader.readNextFrame());
    const auto frame = reader.frame();
    return { frame.f().begin(), frame.f().end() };
}

/*! \brief Sets GMX_MDMODULE_PLUGIN_PATH for the lifetime of the object
 *
 * The test plugins are in the directory TEST_PLUGIN_DIRECTORY.
 */
class PluginPathSetter
{
public:
    explicit PluginPathSetter(const char* path)
    {
        if (const char* previous = std::getenv("GMX_MDMODULE_PLUGIN_PATH"))
        {
            previous_    = previous;
            hadPrevious_ = true;
        }
        setenv("GMX_MDMODULE_PLUGIN_PATH", path, 1);
    }
    ~PluginPathSetter()
    {
        if (hadPrevious_)
        {
            setenv("GMX_MDMODULE_PLUGIN_PATH", previous_.c_str(), 1);
        }
        else
        {
            unsetenv("GMX_MDMODULE_PLUGIN_PATH");
        }
    }
    GMX_DISALLOW_COPY_MOVE_AND_ASSIGN(PluginPathSetter);

private:
    std::string previous_;
    bool        hadPrevious_ = false;
};

using MDModulePluginTest = MdrunTestFixture;

//! The force along z that the plugin adds to every atom
const real c_pluginForceZ = 2.5;

TEST_F(MDModulePluginTest, PluginNamedInMdpAddsItsForce)
{
    if (!GMX_MDMODULE_PLUGINS)
    {
        GTEST_SKIP() << "Needs a build with GMX_MDMODULE_PLUGINS=ON";
    }
    const PluginPathSetter pluginPath(TEST_PLUGIN_DIRECTORY);

    const std::string simulationName = "argon5832";
    const int         nsteps         = 2;
    auto              mdpFieldValues = prepareMdpFieldValues(simulationName, "md", "no", "no");
    mdpFieldValues["nsteps"]         = std::to_string(nsteps);
    mdpFieldValues["nstxout"]        = std::to_string(nsteps);
    mdpFieldValues["nstfout"]        = std::to_string(nsteps);
    mdpFieldValues["nstvout"]        = "0";
    runner_.useTopGroAndNdxFromDatabase(simulationName);

    // Without the plugin
    runner_.useStringAsMdpFile(prepareMdpFileContents(mdpFieldValues));
    ASSERT_EQ(0, runner_.callGrompp());
    const std::string trajectoryWithout = fileManager_.getTemporaryFilePath("without.trr").string();
    runner_.fullPrecisionTrajectoryFileName_ = trajectoryWithout;
    ASSERT_EQ(0, runner_.callMdrun());

    // With the plugin, named in the mdp file
    runner_.useStringAsMdpFile(prepareMdpFileContents(mdpFieldValues)
                               + "\nmdmodule-plugins = test-force-plugin\n"
                               + "test-force-plugin-z = " + std::to_string(c_pluginForceZ) + "\n");
    ASSERT_EQ(0, runner_.callGrompp());
    const std::string trajectoryWith = fileManager_.getTemporaryFilePath("with.trr").string();
    runner_.fullPrecisionTrajectoryFileName_ = trajectoryWith;
    ASSERT_EQ(0, runner_.callMdrun());

    const std::vector<RVec> without = readFirstForces(trajectoryWithout);
    const std::vector<RVec> with    = readFirstForces(trajectoryWith);
    ASSERT_EQ(without.size(), with.size());
    const FloatingPointTolerance tolerance = relativeToleranceAsFloatingPoint(c_pluginForceZ, 1e-4);
    for (size_t atom = 0; atom < with.size(); atom++)
    {
        EXPECT_REAL_EQ_TOL(with[atom][XX], without[atom][XX], tolerance) << "atom " << atom;
        EXPECT_REAL_EQ_TOL(with[atom][YY], without[atom][YY], tolerance) << "atom " << atom;
        EXPECT_REAL_EQ_TOL(with[atom][ZZ] - without[atom][ZZ], c_pluginForceZ, tolerance)
                << "atom " << atom;
    }
}

TEST_F(MDModulePluginTest, ElectricFieldPluginGivesTheForcesOfTheBuiltInModule)
{
    if (!GMX_MDMODULE_PLUGINS)
    {
        GTEST_SKIP() << "Needs a build with GMX_MDMODULE_PLUGINS=ON";
    }
    const PluginPathSetter pluginPath(TEST_PLUGIN_DIRECTORY);

    const std::string simulationName = "spc216";
    const int         nsteps         = 2;
    auto              mdpFieldValues = prepareMdpFieldValues(simulationName, "md", "no", "no");
    mdpFieldValues["nsteps"]         = std::to_string(nsteps);
    mdpFieldValues["nstxout"]        = std::to_string(nsteps);
    mdpFieldValues["nstfout"]        = std::to_string(nsteps);
    mdpFieldValues["nstvout"]        = "0";
    runner_.useTopGroAndNdxFromDatabase(simulationName);

    // A static field of 1 V/nm along z, from the module named in extraMdp
    const auto runWith = [&](const std::string& extraMdp, const std::string& trajectoryName)
    {
        runner_.useStringAsMdpFile(prepareMdpFileContents(mdpFieldValues) + extraMdp);
        EXPECT_EQ(0, runner_.callGrompp());
        const std::string trajectory = fileManager_.getTemporaryFilePath(trajectoryName).string();
        runner_.fullPrecisionTrajectoryFileName_ = trajectory;
        EXPECT_EQ(0, runner_.callMdrun());
        return readFirstForces(trajectory);
    };
    const std::vector<RVec> withoutField = runWith("", "without.trr");
    const std::vector<RVec> builtIn      = runWith("\nelectric-field-z = 1 0 0 0\n", "builtin.trr");
    const std::vector<RVec> plugin       = runWith(
            "\nmdmodule-plugins = electric-field-plugin\nelectric-field-plugin-z = 1 0 0 0\n",
            "plugin.trr");

    ASSERT_EQ(builtIn.size(), withoutField.size());
    ASSERT_EQ(plugin.size(), builtIn.size());
    EXPECT_NE(builtIn[0][ZZ], withoutField[0][ZZ]) << "the field should move the charged atoms";
    const FloatingPointTolerance tolerance = defaultRealTolerance();
    for (size_t atom = 0; atom < plugin.size(); atom++)
    {
        EXPECT_REAL_EQ_TOL(plugin[atom][XX], builtIn[atom][XX], tolerance) << "atom " << atom;
        EXPECT_REAL_EQ_TOL(plugin[atom][YY], builtIn[atom][YY], tolerance) << "atom " << atom;
        EXPECT_REAL_EQ_TOL(plugin[atom][ZZ], builtIn[atom][ZZ], tolerance) << "atom " << atom;
    }
}

TEST_F(MDModulePluginTest, GromppRefusesPluginThatIsNotFound)
{
    if (!GMX_MDMODULE_PLUGINS)
    {
        GTEST_SKIP() << "Needs a build with GMX_MDMODULE_PLUGINS=ON";
    }
    const PluginPathSetter pluginPath(fileManager_.getOutputTempDirectory().string().c_str());

    const std::string simulationName = "argon5832";
    auto              mdpFieldValues = prepareMdpFieldValues(simulationName, "md", "no", "no");
    runner_.useTopGroAndNdxFromDatabase(simulationName);
    runner_.useStringAsMdpFile(prepareMdpFileContents(mdpFieldValues)
                               + "\nmdmodule-plugins = test-force-plugin\n");
    GMX_EXPECT_DEATH_IF_SUPPORTED(
            runner_.callGrompp(),
            ".*test-force-plugin.so is.*not in any of the plugin directories.*");
}

} // namespace
} // namespace test
} // namespace gmx
