/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright 2016- The GROMACS Authors
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
#include "gmxpre.h"

#include "mdmodules.h"

#include <cstdio>

#include <algorithm>
#include <array>
#include <iterator>
#include <memory>
#include <string>
#include <string_view>
#include <utility>
#include <vector>

#include "gromacs/applied_forces/colvars/colvarsMDModule.h"
#include "gromacs/applied_forces/densityfitting/densityfitting.h"
#include "gromacs/applied_forces/electricfield.h"
#include "gromacs/applied_forces/nnpot/nnpot.h"
#include "gromacs/applied_forces/plumed/plumedMDModule.h"
#include "gromacs/applied_forces/qmmm/qmmm.h"
#include "gromacs/fmm/fmm_mdmodule.h"
#include "gromacs/imd/imd.h"
#include "gromacs/mdrunutility/mdmodulesnotifiers.h"
#include "gromacs/mdtypes/iforceprovider.h"
#include "gromacs/mdtypes/imdmodule.h"
#include "gromacs/mdtypes/imdoutputprovider.h"
#include "gromacs/mdtypes/imdpoptionprovider.h"
#include "gromacs/mdtypes/inputrec.h"
#include "gromacs/mdtypes/output_control.h"
#include "gromacs/options/options.h"
#include "gromacs/options/optionsection.h"
#include "gromacs/options/treesupport.h"
#include "gromacs/swap/swapcoords.h"
#include "gromacs/utility/exceptions.h"
#include "gromacs/utility/gmxassert.h"
#include "gromacs/utility/keyvaluetree.h"
#include "gromacs/utility/keyvaluetreebuilder.h"
#include "gromacs/utility/keyvaluetreetransform.h"
#include "gromacs/utility/stringutil.h"

struct gmx_output_env_t;

namespace gmx
{

namespace
{

//! Creates a built-in module
using ModuleFactory = std::unique_ptr<IMDModule> (*)();

/*! \brief The built-in modules, in registration order.
 *
 * Modules are called in this order. Their mdp options are grouped by
 * section in the same order, so reordering the table changes the layout
 * of mdp output and of the parameters in tpr files. The order reproduces
 * the mdp layout from before the table existed: the applied-forces
 * modules, fast multipole method, output control, then the modules
 * without mdp options.
 */
const std::array<std::pair<std::string_view, ModuleFactory>, 10> c_builtInModules = { {
        { ElectricFieldModuleInfo::sc_name, &ElectricFieldModuleInfo::create },
        { DensityFittingModuleInfo::sc_name, &DensityFittingModuleInfo::create },
        { QMMMModuleInfo::sc_name, &QMMMModuleInfo::create },
        { ColvarsModuleInfo::sc_name, &ColvarsModuleInfo::create },
        { NNPotModuleInfo::sc_name, &NNPotModuleInfo::create },
        { FmmModuleInfo::sc_name, &FmmModuleInfo::create },
        { OutputControlModuleInfo::sc_name, &OutputControlModuleInfo::create },
        { InteractiveMolecularDynamicsModuleInfo::sc_name, &InteractiveMolecularDynamicsModuleInfo::create },
        { SwapCoordinatesModuleInfo::sc_name, &SwapCoordinatesModuleInfo::create },
        { PlumedModuleInfo::sc_name, &PlumedModuleInfo::create },
} };

//! A module and the name it was registered under
struct RegisteredModule
{
    //! Unique name of the module
    std::string name;
    //! The module
    std::shared_ptr<IMDModule> module;
};

//! The mdp option providers that share one mdp section
struct MdpSection
{
    //! Name of the section, empty for the top level
    std::string name;
    //! Providers in the section, in registration order
    std::vector<IMdpOptionProvider*> providers;
};

} // namespace

class MDModules::Impl : public IMDOutputProvider
{
public:
    Impl()
    {
        for (const auto& [name, create] : c_builtInModules)
        {
            add(name, create());
        }
    }

    //! Registers \p module under \p name, after all modules registered so far
    void add(std::string_view name, std::shared_ptr<IMDModule> module)
    {
        if (find(name) != nullptr)
        {
            GMX_THROW(APIError(
                    formatString("Module with name %s already added", std::string(name).c_str())));
        }
        modules_.push_back({ std::string(name), std::move(module) });
    }

    //! Returns the module registered under \p name, or nullptr
    IMDModule* find(std::string_view name) const
    {
        const auto found =
                std::find_if(modules_.begin(),
                             modules_.end(),
                             [name](const RegisteredModule& entry) { return entry.name == name; });
        return found != modules_.end() ? found->module.get() : nullptr;
    }

    //! Returns the output control module, which every MDModules has
    IMDModule* outputControlModule() const
    {
        IMDModule* module = find(OutputControlModuleInfo::sc_name);
        GMX_RELEASE_ASSERT(module, "The output control module is always registered");
        return module;
    }

    /*! \brief Returns the mdp sections of all modules.
     *
     * Sections are in the order their first module was registered, and
     * the providers within a section in registration order. Modules
     * without an mdp option provider are skipped.
     */
    std::vector<MdpSection> mdpSections() const
    {
        std::vector<MdpSection> sections;
        for (const auto& entry : modules_)
        {
            IMdpOptionProvider* provider = entry.module->mdpOptionProvider();
            if (provider == nullptr)
            {
                continue;
            }
            const std::string_view sectionName = provider->mdpSectionName();

            auto section = std::find_if(sections.begin(),
                                        sections.end(),
                                        [sectionName](const MdpSection& existing)
                                        { return existing.name == sectionName; });
            if (section == sections.end())
            {
                sections.push_back({ std::string(sectionName), {} });
                section = std::prev(sections.end());
            }
            section->providers.push_back(provider);
        }
        return sections;
    }

    void makeModuleOptions(Options* options) const
    {
        for (const auto& section : mdpSections())
        {
            if (section.name.empty())
            {
                for (IMdpOptionProvider* provider : section.providers)
                {
                    provider->initMdpOptions(options);
                }
            }
            else
            {
                auto sectionOptions = options->addSection(OptionSection(section.name.c_str()));
                for (IMdpOptionProvider* provider : section.providers)
                {
                    provider->initMdpOptions(&sectionOptions);
                }
            }
        }
    }

    // From IMDOutputProvider
    void initOutput(FILE* fplog, int nfile, const t_filenm fnm[], bool bAppendFiles, const gmx_output_env_t* oenv) override
    {
        for (const auto& entry : modules_)
        {
            if (entry.module->outputProvider())
            {
                entry.module->outputProvider()->initOutput(fplog, nfile, fnm, bAppendFiles, oenv);
            }
        }
    }
    void finishOutput() override
    {
        for (const auto& entry : modules_)
        {
            if (entry.module->outputProvider())
            {
                entry.module->outputProvider()->finishOutput();
            }
        }
    }

    /*! \brief Manages callbacks and notifies the MD modules.
     *
     * \note The notifier must be constructed before the modules and shall
     *       not be destructed before the modules are destructed.
     */
    MDModulesNotifiers notifiers_;

    std::unique_ptr<ForceProviders> forceProviders_;

    /*! \brief Registered modules, in registration order */
    std::vector<RegisteredModule> modules_;
};

MDModules::MDModules() : impl_(new Impl) {}

MDModules::~MDModules() {}

void MDModules::initMdpTransform(IKeyValueTreeTransformRules* rules)
{
    for (const auto& section : impl_->mdpSections())
    {
        if (section.name.empty())
        {
            for (IMdpOptionProvider* provider : section.providers)
            {
                provider->initMdpTransform(rules);
            }
        }
        else
        {
            auto scope = rules->scopedTransform("/" + section.name);
            for (IMdpOptionProvider* provider : section.providers)
            {
                provider->initMdpTransform(scope.rules());
            }
        }
    }
}

void MDModules::buildMdpOutput(KeyValueTreeObjectBuilder* builder)
{
    // The order here is the order of the fields in the mdp output file,
    // so it follows the sections, as the options do.
    for (const auto& section : impl_->mdpSections())
    {
        for (const IMdpOptionProvider* provider : section.providers)
        {
            provider->buildMdpOutput(builder);
        }
    }
}

void MDModules::assignOptionsToModules(const KeyValueTreeObject&  params,
                                       IKeyValueTreeErrorHandler* errorHandler,
                                       t_inputrec*                inputrec,
                                       inputrec_strings*          preprocessingStrings)
{
    // Configure module targets if inputrec provided
    if (inputrec)
    {
        // Route OutputControl sub-object to output-control module
        IMDModule* outputControlModule = impl_->outputControlModule();
        setOutputControlTarget(outputControlModule, &inputrec->outputControl);

        // Set preprocessing strings if provided
        // For grompp, both inputrec and preprocessingStrings must be valid
        // For mdrun/tools reading TPR, only inputrec is needed
        if (preprocessingStrings)
        {
            setOutputControlPreprocessingStrings(outputControlModule, preprocessingStrings);
        }

        // Future: route other sub-objects here
    }

    Options moduleOptions;
    impl_->makeModuleOptions(&moduleOptions);
    // The actual output is in the data fields of the modules that
    // were set up in the module options.
    assignOptionsFromKeyValueTree(&moduleOptions, params, errorHandler);
}

void MDModules::adjustInputrecBasedOnModules(t_inputrec* ir, bool routeSubObjects, inputrec_strings* preprocessingStrings)
{
    GMX_RELEASE_ASSERT(ir, "adjustInputrecBasedOnModules() called with null inputrec");
    GMX_RELEASE_ASSERT(ir->params, "adjustInputrecBasedOnModules() called but ir->params is null");

    // Configure module targets if requested
    if (routeSubObjects)
    {
        // Route OutputControl sub-object to output-control module
        IMDModule* outputControlModule = impl_->outputControlModule();
        setOutputControlTarget(outputControlModule, &ir->outputControl);

        // Set preprocessing strings if provided so preprocessing-only options are registered
        if (preprocessingStrings)
        {
            setOutputControlPreprocessingStrings(outputControlModule, preprocessingStrings);
        }

        // Future: route other sub-objects here
    }

    Options moduleOptions;
    impl_->makeModuleOptions(&moduleOptions);

    checkForUnknownOptionsInKeyValueTree(*ir->params, moduleOptions);

    std::unique_ptr<KeyValueTreeObject> adjustedParams(
            new KeyValueTreeObject(adjustKeyValueTreeFromOptions(*ir->params, moduleOptions)));

    delete ir->params;
    ir->params = adjustedParams.release();
}

IMDOutputProvider* MDModules::outputProvider()
{
    return impl_.get();
}

ForceProviders* MDModules::initForceProviders(gmx_wallcycle* wallCycle)
{
    GMX_RELEASE_ASSERT(impl_->forceProviders_ == nullptr,
                       "Force providers initialized multiple times");
    impl_->forceProviders_ = std::make_unique<ForceProviders>(wallCycle);
    for (const auto& entry : impl_->modules_)
    {
        entry.module->initForceProviders(impl_->forceProviders_.get());
    }
    return impl_->forceProviders_.get();
}

void MDModules::subscribeToPreProcessingNotifications()
{
    for (const auto& entry : impl_->modules_)
    {
        entry.module->subscribeToPreProcessingNotifications(&impl_->notifiers_);
    }
}

void MDModules::subscribeToSimulationSetupNotifications()
{
    const auto& setupNotifier          = impl_->notifiers_.simulationSetupNotifier_;
    int         numMultisimSubscribers = setupNotifier.numSubscribers<const gmx_multisim_t*>();
    for (const auto& entry : impl_->modules_)
    {
        entry.module->subscribeToSimulationSetupNotifications(&impl_->notifiers_);
        const int  numSubscribersAfter  = setupNotifier.numSubscribers<const gmx_multisim_t*>();
        const bool subscribedToMultisim = numSubscribersAfter > numMultisimSubscribers;
        numMultisimSubscribers          = numSubscribersAfter;
        if (subscribedToMultisim && !entry.module->descriptor().sharesStateAcrossSimulations)
        {
            GMX_THROW(APIError(
                    formatString("Module %s subscribed to the multi-simulation notification "
                                 "without declaring that it shares state across simulations",
                                 entry.name.c_str())));
        }
    }

    std::vector<InteractionFunction> energyTerms;
    for (const auto& entry : impl_->modules_)
    {
        const auto terms = entry.module->descriptor().energyTerms;
        energyTerms.insert(energyTerms.end(), terms.begin(), terms.end());
    }
    impl_->notifiers_.simulationSetupNotifier_.subscribe([energyTerms](MDModulesEnergyTerms* request)
                                                         { request->terms_ = energyTerms; });
}

void MDModules::subscribeToSimulationRunNotifications()
{
    for (const auto& entry : impl_->modules_)
    {
        entry.module->subscribeToSimulationRunNotifications(&impl_->notifiers_);
    }
}

void MDModules::add(std::string_view name, std::shared_ptr<IMDModule> module)
{
    impl_->add(name, std::move(module));
}

bool MDModules::anyModuleSharesStateAcrossSimulations() const
{
    return std::any_of(impl_->modules_.begin(),
                       impl_->modules_.end(),
                       [](const RegisteredModule& entry)
                       { return entry.module->descriptor().sharesStateAcrossSimulations; });
}

bool MDModules::anyModuleProvidesDirectCoulomb() const
{
    const RegisteredModule* provider = nullptr;
    for (const auto& entry : impl_->modules_)
    {
        if (entry.module->descriptor().providesDirectCoulomb)
        {
            if (provider != nullptr)
            {
                GMX_THROW(APIError(
                        formatString("Modules %s and %s both declare that they provide the "
                                     "direct Coulomb interactions",
                                     provider->name.c_str(),
                                     entry.name.c_str())));
            }
            provider = &entry;
        }
    }
    return provider != nullptr;
}

const MDModulesNotifiers& MDModules::notifiers()
{
    return impl_->notifiers_;
}

} // namespace gmx
