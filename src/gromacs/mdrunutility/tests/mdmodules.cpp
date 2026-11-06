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
 * Tests registration of modules with MDModules
 *
 * \author Amr Alhossary <amr_alhossary@hotmail.com>
 * \ingroup module_mdrunutility
 */
#include "gmxpre.h"

#include "gromacs/mdrun/mdmodules.h"

#include <memory>
#include <string>
#include <string_view>

#include <gtest/gtest.h>

#include "gromacs/applied_forces/electricfield.h"
#include "gromacs/mdrunutility/mdmodulesnotifiers.h"
#include "gromacs/mdtypes/iforceprovider.h"
#include "gromacs/mdtypes/imdmodule.h"
#include "gromacs/mdtypes/imdoutputprovider.h"
#include "gromacs/mdtypes/imdpoptionprovider.h"
#include "gromacs/mdtypes/inputrec.h"
#include "gromacs/options/basicoptions.h"
#include "gromacs/options/ioptionscontainerwithsections.h"
#include "gromacs/options/optionsection.h"
#include "gromacs/topology/ifunc.h"
#include "gromacs/utility/exceptions.h"
#include "gromacs/utility/keyvaluetree.h"
#include "gromacs/utility/keyvaluetreebuilder.h"
#include "gromacs/utility/keyvaluetreetransform.h"
#include "gromacs/utility/strconvert.h"
#include "gromacs/utility/stringcompare.h"

namespace gmx
{
namespace test
{
namespace
{

//! A module with one integer mdp option, placed in a given mdp section
class SectionedTestModule final : public IMDModule, public IMdpOptionProvider
{
public:
    explicit SectionedTestModule(std::string_view sectionName) : sectionName_(sectionName) {}

    // From IMDModule
    IMdpOptionProvider* mdpOptionProvider() override { return this; }
    IMDOutputProvider*  outputProvider() override { return nullptr; }
    void                initForceProviders(ForceProviders* /* forceProviders */) override {}
    void subscribeToPreProcessingNotifications(MDModulesNotifiers* /* notifiers */) override {}
    void subscribeToSimulationSetupNotifications(MDModulesNotifiers* /* notifiers */) override {}
    void subscribeToSimulationRunNotifications(MDModulesNotifiers* /* notifiers */) override {}

    // From IMdpOptionProvider
    std::string_view mdpSectionName() const override { return sectionName_; }
    void             initMdpTransform(IKeyValueTreeTransformRules* rules) override
    {
        rules->addRule()
                .from<std::string>("/test-module-value")
                .to<int>("/test-module/value")
                .transformWith(&fromStdString<int>);
    }
    void initMdpOptions(IOptionsContainerWithSections* options) override
    {
        auto section = options->addSection(OptionSection("test-module"));
        section.addOption(IntegerOption("value").store(&value_));
    }
    void buildMdpOutput(KeyValueTreeObjectBuilder* builder) const override
    {
        builder->addValue<std::string>("test-module-value", toString(value_));
    }

    //! The value assigned from mdp input
    int value_ = 0;

private:
    std::string sectionName_;
};

//! A module that implements none of the optional interfaces
class EmptyTestModule final : public IMDModule
{
public:
    IMdpOptionProvider* mdpOptionProvider() override { return nullptr; }
    IMDOutputProvider*  outputProvider() override { return nullptr; }
    void                initForceProviders(ForceProviders* /* forceProviders */) override {}
    void subscribeToPreProcessingNotifications(MDModulesNotifiers* /* notifiers */) override {}
    void subscribeToSimulationSetupNotifications(MDModulesNotifiers* /* notifiers */) override {}
    void subscribeToSimulationRunNotifications(MDModulesNotifiers* /* notifiers */) override {}
};

//! A module that declares what the test puts in \c descriptor_, and may subscribe to the multi-simulation notification
class DeclaringTestModule final : public IMDModule
{
public:
    explicit DeclaringTestModule(bool subscribesToMultisim) :
        subscribesToMultisim_(subscribesToMultisim)
    {
    }

    MDModuleDescriptor  descriptor() const override { return descriptor_; }
    IMdpOptionProvider* mdpOptionProvider() override { return nullptr; }
    IMDOutputProvider*  outputProvider() override { return nullptr; }
    void                initForceProviders(ForceProviders* /* forceProviders */) override {}
    void subscribeToPreProcessingNotifications(MDModulesNotifiers* /* notifiers */) override {}
    void subscribeToSimulationSetupNotifications(MDModulesNotifiers* notifiers) override
    {
        if (subscribesToMultisim_)
        {
            notifiers->simulationSetupNotifier_.subscribe([](const gmx_multisim_t* /* ms */) {});
        }
    }
    void subscribeToSimulationRunNotifications(MDModulesNotifiers* /* notifiers */) override {}

    //! What the module declares
    MDModuleDescriptor descriptor_;

private:
    bool subscribesToMultisim_;
};

//! A force provider that adds no forces
class NullForceProvider final : public IForceProvider
{
public:
    void calculateForces(const ForceProviderInput& /* forceProviderInput */,
                         ForceProviderOutput* /* forceProviderOutput */) override
    {
    }
};

//! A module that adds a force provider for halo atoms and declares what the test puts in \c descriptor_
class HaloForceTestModule final : public IMDModule
{
public:
    MDModuleDescriptor  descriptor() const override { return descriptor_; }
    IMdpOptionProvider* mdpOptionProvider() override { return nullptr; }
    IMDOutputProvider*  outputProvider() override { return nullptr; }
    void                initForceProviders(ForceProviders* forceProviders) override
    {
        forceProviders->addForceProvider(&forceProvider_, "", ForceProviderAtoms::HomeAndHalo);
    }
    void subscribeToPreProcessingNotifications(MDModulesNotifiers* /* notifiers */) override {}
    void subscribeToSimulationSetupNotifications(MDModulesNotifiers* /* notifiers */) override {}
    void subscribeToSimulationRunNotifications(MDModulesNotifiers* /* notifiers */) override {}

    //! What the module declares
    MDModuleDescriptor descriptor_;

private:
    NullForceProvider forceProvider_;
};

//! Transforms flat mdp-style \c input with the rules of all modules in \c mdModules
KeyValueTreeObject transformMdpInput(MDModules* mdModules, const KeyValueTreeObject& input)
{
    KeyValueTreeTransformer transform;
    transform.rules()->addRule().keyMatchType("/", StringCompareType::CaseAndDashInsensitive);
    mdModules->initMdpTransform(transform.rules());
    return transform.transform(input, nullptr).object();
}

//! Returns flat mdp-style input that sets the test module value to \c value
KeyValueTreeObject testModuleMdpInput(const std::string& value)
{
    KeyValueTreeBuilder mdpValues;
    mdpValues.rootObject().addValue<std::string>("test-module-value", value);
    return mdpValues.build();
}

TEST(MDModulesTest, AddedModuleGetsItsOwnMdpSection)
{
    MDModules mdModules;
    auto      module = std::make_shared<SectionedTestModule>("test-section");
    mdModules.add("test", module);

    const KeyValueTreeObject params = transformMdpInput(&mdModules, testModuleMdpInput("7"));
    ASSERT_TRUE(params.keyExists("test-section"));
    EXPECT_TRUE(params["test-section"].asObject().keyExists("test-module"));

    t_inputrec ir;
    mdModules.assignOptionsToModules(params, nullptr, &ir);
    EXPECT_EQ(7, module->value_);
}

TEST(MDModulesTest, AddedModuleJoinsAnExistingMdpSection)
{
    MDModules mdModules;
    auto      module = std::make_shared<SectionedTestModule>("applied-forces");
    mdModules.add("test", module);

    const KeyValueTreeObject params = transformMdpInput(&mdModules, testModuleMdpInput("3"));
    ASSERT_TRUE(params.keyExists("applied-forces"));
    EXPECT_TRUE(params["applied-forces"].asObject().keyExists("test-module"));

    t_inputrec ir;
    mdModules.assignOptionsToModules(params, nullptr, &ir);
    EXPECT_EQ(3, module->value_);
}

TEST(MDModulesTest, AddedModuleWritesMdpOutputAfterBuiltInSections)
{
    MDModules mdModules;
    auto      module = std::make_shared<SectionedTestModule>("test-section");
    module->value_   = 5;
    mdModules.add("test", module);

    KeyValueTreeBuilder builder;
    auto                builderObject = builder.rootObject();
    mdModules.buildMdpOutput(&builderObject);
    const KeyValueTreeObject output = builder.build();

    ASSERT_TRUE(output.keyExists("test-module-value"));
    EXPECT_EQ("5", output["test-module-value"].cast<std::string>());
    EXPECT_EQ("test-module-value", output.properties().back().key());
}

TEST(MDModulesTest, AcceptsModuleWithoutOptionalInterfaces)
{
    MDModules mdModules;
    mdModules.add("empty", std::make_shared<EmptyTestModule>());

    const KeyValueTreeObject params = transformMdpInput(&mdModules, KeyValueTreeBuilder().build());
    t_inputrec               ir;
    EXPECT_NO_THROW(mdModules.assignOptionsToModules(params, nullptr, &ir));

    KeyValueTreeBuilder builder;
    auto                builderObject = builder.rootObject();
    EXPECT_NO_THROW(mdModules.buildMdpOutput(&builderObject));
    EXPECT_NO_THROW(mdModules.outputProvider()->finishOutput());
}

TEST(MDModulesTest, RejectsDuplicateModuleName)
{
    MDModules mdModules;
    mdModules.add("test", std::make_shared<EmptyTestModule>());
    EXPECT_THROW(mdModules.add("test", std::make_shared<EmptyTestModule>()), APIError);
    EXPECT_THROW(mdModules.add(ElectricFieldModuleInfo::sc_name, std::make_shared<EmptyTestModule>()),
                 APIError);
}

TEST(MDModulesTest, RejectsMultisimSubscriptionWithoutDeclaration)
{
    MDModules mdModules;
    mdModules.add("test", std::make_shared<DeclaringTestModule>(true));

    const KeyValueTreeObject params = transformMdpInput(&mdModules, KeyValueTreeBuilder().build());
    t_inputrec               ir;
    mdModules.assignOptionsToModules(params, nullptr, &ir);
    EXPECT_FALSE(mdModules.anyModuleSharesStateAcrossSimulations());
    EXPECT_THROW(mdModules.subscribeToSimulationSetupNotifications(), APIError);
}

TEST(MDModulesTest, AcceptsDeclaredMultisimSubscription)
{
    MDModules mdModules;
    auto      module                                 = std::make_shared<DeclaringTestModule>(true);
    module->descriptor_.sharesStateAcrossSimulations = true;
    mdModules.add("test", module);

    const KeyValueTreeObject params = transformMdpInput(&mdModules, KeyValueTreeBuilder().build());
    t_inputrec               ir;
    mdModules.assignOptionsToModules(params, nullptr, &ir);
    EXPECT_NO_THROW(mdModules.subscribeToSimulationSetupNotifications());
    EXPECT_TRUE(mdModules.anyModuleSharesStateAcrossSimulations());
}

TEST(MDModulesTest, RejectsHaloForceProviderWithoutDeclaration)
{
    MDModules mdModules;
    mdModules.add("test", std::make_shared<HaloForceTestModule>());
    EXPECT_THROW(mdModules.initForceProviders(nullptr), APIError);
}

TEST(MDModulesTest, AcceptsDeclaredHaloForceProvider)
{
    MDModules mdModules;
    auto      module                     = std::make_shared<HaloForceTestModule>();
    module->descriptor_.writesHaloForces = true;
    mdModules.add("test", module);
    const ForceProviders* forceProviders = mdModules.initForceProviders(nullptr);
    EXPECT_EQ(forceProviders->numHaloForceProviders(), 1);
}

TEST(MDModulesTest, RejectsTwoDirectCoulombProviders)
{
    MDModules mdModules;
    EXPECT_FALSE(mdModules.anyModuleProvidesDirectCoulomb());

    auto first                               = std::make_shared<DeclaringTestModule>(false);
    first->descriptor_.providesDirectCoulomb = true;
    mdModules.add("first", first);
    EXPECT_TRUE(mdModules.anyModuleProvidesDirectCoulomb());

    auto second                               = std::make_shared<DeclaringTestModule>(false);
    second->descriptor_.providesDirectCoulomb = true;
    mdModules.add("second", second);
    EXPECT_THROW(mdModules.anyModuleProvidesDirectCoulomb(), APIError);
}

TEST(MDModulesTest, ReportsDeclaredEnergyTerms)
{
    MDModules mdModules;
    auto      module                = std::make_shared<DeclaringTestModule>(false);
    module->descriptor_.energyTerms = { InteractionFunction::DensityFitting,
                                        InteractionFunction::NeuralNetworkPotentialEnergy };
    mdModules.add("test", module);

    const KeyValueTreeObject params = transformMdpInput(&mdModules, KeyValueTreeBuilder().build());
    t_inputrec               ir;
    mdModules.assignOptionsToModules(params, nullptr, &ir);
    mdModules.subscribeToSimulationSetupNotifications();

    MDModulesEnergyTerms energyTerms;
    mdModules.notifiers().simulationSetupNotifier_.notify(&energyTerms);
    EXPECT_EQ(module->descriptor_.energyTerms, energyTerms.terms_);
}

} // namespace
} // namespace test
} // namespace gmx
