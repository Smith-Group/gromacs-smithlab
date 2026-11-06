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
 * An MDModules plugin with one mdp option, for the tests of the loader
 *
 * \author Amr Alhossary <amr_alhossary@hotmail.com>
 * \ingroup module_mdrunutility
 */
#include "gmxpre.h"

#include <memory>
#include <string>
#include <string_view>

#include "gromacs/mdtypes/imdmodule.h"
#include "gromacs/mdtypes/imdpoptionprovider.h"
#include "gromacs/mdtypes/mdmodulepluginentrypoints.h"
#include "gromacs/options/basicoptions.h"
#include "gromacs/options/ioptionscontainerwithsections.h"
#include "gromacs/options/optionsection.h"
#include "gromacs/utility/keyvaluetreebuilder.h"
#include "gromacs/utility/keyvaluetreetransform.h"
#include "gromacs/utility/strconvert.h"

#ifndef TEST_PLUGIN_NAME
#    error "TEST_PLUGIN_NAME must name the module"
#endif

namespace gmx
{
namespace test
{
namespace
{

//! A module with one integer mdp option, test-plugin-value, in its own mdp section
class TestPluginModule final : public IMDModule, public IMdpOptionProvider
{
public:
    // From IMDModule
    IMdpOptionProvider* mdpOptionProvider() override { return this; }
    IMDOutputProvider*  outputProvider() override { return nullptr; }
    void                initForceProviders(ForceProviders* /* forceProviders */) override {}
    void subscribeToPreProcessingNotifications(MDModulesNotifiers* /* notifiers */) override {}
    void subscribeToSimulationSetupNotifications(MDModulesNotifiers* /* notifiers */) override {}
    void subscribeToSimulationRunNotifications(MDModulesNotifiers* /* notifiers */) override {}

    // From IMdpOptionProvider
    std::string_view mdpSectionName() const override { return "test-plugin"; }
    void             initMdpTransform(IKeyValueTreeTransformRules* rules) override
    {
        rules->addRule()
                .from<std::string>("/test-plugin-value")
                .to<int>("/test-plugin/value")
                .transformWith(&fromStdString<int>);
    }
    void initMdpOptions(IOptionsContainerWithSections* options) override
    {
        auto section = options->addSection(OptionSection("test-plugin"));
        section.addOption(IntegerOption("value").store(&value_));
    }
    void buildMdpOutput(KeyValueTreeObjectBuilder* builder) const override
    {
        builder->addValue<std::string>("test-plugin-value", toString(value_));
    }

private:
    int value_ = 0;
};

//! Creates the module of the plugin
std::unique_ptr<IMDModule> createTestPluginModule()
{
    return std::make_unique<TestPluginModule>();
}

} // namespace
} // namespace test
} // namespace gmx

GMX_MDMODULE_PLUGIN(TEST_PLUGIN_NAME, "1.2", ::gmx::test::createTestPluginModule)
