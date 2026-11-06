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
 * Tests the loading of MDModules plugins
 *
 * \author Amr Alhossary <amr_alhossary@hotmail.com>
 * \ingroup module_mdrunutility
 */
#include "gmxpre.h"

#include "config.h"

#include <filesystem>
#include <string>
#include <vector>

#include <gmock/gmock.h>
#include <gtest/gtest.h>

#include "gromacs/mdrun/mdmodulepluginloader.h"
#include "gromacs/mdrun/mdmodules.h"
#include "gromacs/mdtypes/inputrec.h"
#include "gromacs/mdtypes/mdmoduleplugin.h"
#include "gromacs/utility/exceptions.h"
#include "gromacs/utility/keyvaluetree.h"
#include "gromacs/utility/keyvaluetreebuilder.h"
#include "gromacs/utility/keyvaluetreetransform.h"
#include "gromacs/utility/stringcompare.h"
#include "gromacs/utility/stringutil.h"

namespace gmx
{
namespace test
{
namespace
{

#if GMX_MDMODULE_PLUGINS

using ::testing::HasSubstr;

//! The directory with the test plugins
const std::vector<std::filesystem::path> c_testPluginSearchPath = { TEST_PLUGIN_DIRECTORY };

//! Returns the message of the InvalidInputError that loading plugin \p name throws
std::string loadingError(const std::string& name)
{
    MDModules mdModules;
    try
    {
        const std::vector<std::string> names = { name };
        mdModules.addPlugins(names, c_testPluginSearchPath);
    }
    catch (const InvalidInputError& error)
    {
        return error.what();
    }
    return "no error";
}

TEST(MDModulePluginsTest, LoadsPluginWithItsMdpOptions)
{
    MDModules                      mdModules;
    const std::vector<std::string> names = { "test-plugin" };
    mdModules.addPlugins(names, c_testPluginSearchPath);
    ASSERT_EQ(mdModules.plugins().size(), 1);
    EXPECT_EQ(mdModules.plugins()[0].name, "test-plugin");
    EXPECT_EQ(mdModules.plugins()[0].version, "1.2");

    KeyValueTreeBuilder mdpValues;
    mdpValues.rootObject().addValue<std::string>("test-plugin-value", "4");
    KeyValueTreeTransformer transform;
    transform.rules()->addRule().keyMatchType("/", StringCompareType::CaseAndDashInsensitive);
    mdModules.initMdpTransform(transform.rules());
    const KeyValueTreeObject params = transform.transform(mdpValues.build(), nullptr).object();
    t_inputrec               ir;
    mdModules.assignOptionsToModules(params, nullptr, &ir);

    KeyValueTreeBuilder builder;
    auto                builderObject = builder.rootObject();
    mdModules.buildMdpOutput(&builderObject);
    const KeyValueTreeObject output = builder.build();
    ASSERT_TRUE(output.keyExists("test-plugin-value"));
    EXPECT_EQ(output["test-plugin-value"].cast<std::string>(), "4");
}

TEST(MDModulePluginsTest, LoadsAPluginNamedTwiceOnce)
{
    MDModules                      mdModules;
    const std::vector<std::string> names = { "test-plugin", "test-plugin" };
    mdModules.addPlugins(names, c_testPluginSearchPath);
    EXPECT_EQ(mdModules.plugins().size(), 1);
}

TEST(MDModulePluginsTest, RecordsLoadedPlugins)
{
    MDModules                      mdModules;
    const std::vector<std::string> names = { "test-plugin" };
    mdModules.addPlugins(names, c_testPluginSearchPath);
    KeyValueTreeBuilder builder;
    mdModules.writePluginRecord(builder.rootObject());
    const KeyValueTreeObject record = builder.build();
    EXPECT_EQ(MDModules::recordedPluginNames(record), std::vector<std::string>{ "test-plugin" });
    // The test plugin is built with the headers of this host
    const GmxHostApiInfo host = hostApiInfo();
    const KeyValueTreeObject& entry = record["mdmodule-plugins"].asObject()["test-plugin"].asObject();
    EXPECT_EQ(entry["version"].cast<std::string>(), "1.2");
    EXPECT_EQ(entry["api-major"].cast<int64_t>(), static_cast<int64_t>(host.apiMajor));
    EXPECT_EQ(entry["api-hash"].cast<std::string>(),
              formatString("%016llx", static_cast<unsigned long long>(host.apiHash)));
}

TEST(MDModulePluginsTest, RefusesMissingPlugin)
{
    EXPECT_THAT(loadingError("missing-plugin"), HasSubstr("is not in any of the plugin directories"));
}

TEST(MDModulePluginsTest, RefusesPluginThatRefusesTheHost)
{
    const std::string error = loadingError("incompatible-plugin");
    EXPECT_THAT(error, HasSubstr("it refuses every host"));
    EXPECT_THAT(error, HasSubstr("This GROMACS:"));
}

TEST(MDModulePluginsTest, RefusesPluginWithAnotherName)
{
    EXPECT_THAT(loadingError("misnamed-plugin"), HasSubstr("is the plugin of module 'another-name'"));
}

TEST(MDModulePluginsTest, RefusesLibraryWithoutEntryPoints)
{
    EXPECT_THAT(loadingError("not-a-plugin"), HasSubstr("is not an MDModules plugin"));
}

TEST(MDModulePluginsTest, RefusesPath)
{
    EXPECT_THAT(loadingError("../test-plugin"), HasSubstr("not by a path"));
}

#else

TEST(MDModulePluginsTest, RefusesPluginsWithoutPluginSupport)
{
    MDModules                      mdModules;
    const std::vector<std::string> names = { "test-plugin" };
    EXPECT_THROW(mdModules.addPlugins(names, {}), InvalidInputError);
}

#endif

TEST(MDModulePluginsTest, RecordsNothingWithoutPlugins)
{
    MDModules           mdModules;
    KeyValueTreeBuilder builder;
    mdModules.writePluginRecord(builder.rootObject());
    const KeyValueTreeObject record = builder.build();
    EXPECT_TRUE(record.properties().empty());
    EXPECT_TRUE(MDModules::recordedPluginNames(record).empty());
}

TEST(MDModulePluginsTest, SearchPathComesFromTheEnvironment)
{
    const char*       previous      = std::getenv("GMX_MDMODULE_PLUGIN_PATH");
    const std::string previousValue = previous != nullptr ? previous : "";
    setenv("GMX_MDMODULE_PLUGIN_PATH", "/first:/second", 1);
    const std::vector<std::filesystem::path> searchPath = mdModulePluginSearchPath();
    if (previous != nullptr)
    {
        setenv("GMX_MDMODULE_PLUGIN_PATH", previousValue.c_str(), 1);
    }
    else
    {
        unsetenv("GMX_MDMODULE_PLUGIN_PATH");
    }
    EXPECT_EQ(searchPath, (std::vector<std::filesystem::path>{ "/first", "/second" }));
}

} // namespace
} // namespace test
} // namespace gmx
