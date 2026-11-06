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
/*! \libinternal \file
 * \brief
 * Declares gmx::MDModules.
 *
 * \author Teemu Murtola <teemu.murtola@gmail.com>
 * \inlibraryapi
 * \ingroup module_mdrun
 */
#ifndef GMX_MDRUN_MDMODULES_H
#define GMX_MDRUN_MDMODULES_H

#include <filesystem>
#include <functional>
#include <memory>
#include <string>
#include <string_view>
#include <vector>

#include "gromacs/mdrun/mdmodulepluginloader.h"
#include "gromacs/utility/arrayref.h"

struct gmx_wallcycle;
struct t_inputrec;

namespace gmx
{

struct inputrec_strings;

class ForceProviders;
class KeyValueTreeObjectBuilder;
class KeyValueTreeObject;
class IKeyValueTreeErrorHandler;
class IKeyValueTreeTransformRules;
class IMDModule;
class IMDOutputProvider;
struct MDModulesNotifiers;
struct OutputControl;

/*! \libinternal \brief
 * Makes a module that the caller of mdrun adds to the built-in ones
 *
 * With thread-MPI each rank has its own MDModules, so each rank makes
 * its own module.
 *
 * \inlibraryapi
 * \ingroup module_mdrunutility
 */
struct MDModuleFactory
{
    //! The name the module is added with, see MDModules::add()
    std::string name;
    //! Makes the module
    std::function<std::shared_ptr<IMDModule>()> make;
};

/*! \libinternal \brief
 * Manages the collection of all modules used for mdrun.
 *
 * This class acts as a central place for constructing modules for mdrun
 * and wiring up dependencies between them.  This class should be the only
 * place that contains the full list of modules, although in the future, some
 * code (e.g., in tools) may benefit from the ability to only create one or a
 * few modules and use them.
 *
 * The general idea is that each module takes care of its own data rather than
 * (for example) mdrun having to know about all the details of each type of force calculation.
 * Initially this is applied for simple things like electric field calculations
 * but later more complex forces will be supported too.
 *
 * Currently, where the set of modules needs to be accessed, either a pointer
 * to MDModules is passed around, or an instance of IMDOutputProvider or
 * ForceProviders returned from MDModules.  These objects returned from
 * MDModules call the corresponding methods in the relevant modules.
 * In the future, some additional logic may need to be introduced at
 * the call sites that can also influence the signature of the methods,
 * similar to what ForceProviders already does for force computation.
 *
 * The assignOptionsToModules() and adjustInputrecBasedOnModules() methods of
 * this class also take responsibility for wiring up the options (and their
 * defaults) for each module.
 *
 * \inlibraryapi
 * \ingroup module_mdrunutility
 */
class MDModules
{
public:
    MDModules();
    ~MDModules();

    /*! \brief
     * Initializes a transform from mdp values to sectioned options.
     *
     * \see IMdpOptionProvider::initMdpTransform()
     *
     * Initializes the combined transform from all modules.
     */
    void initMdpTransform(IKeyValueTreeTransformRules* rules);

    /*! \brief Initializes a builder of flat mdp-style key-value pairs
     * suitable for output.
     *
     * If used as input to initMdpTransform(), the key-value pairs
     * resulting from this function would leave the module
     * settings unchanged.
     *
     * Once the transition from mdp to key-value input is
     * complete, this method will probably not exist.
     */
    void buildMdpOutput(KeyValueTreeObjectBuilder* builder);

    /*! \brief
     * Sets input parameters from `params` for each module.
     *
     * If \c inputrec is provided, also configures modules to write to
     * the appropriate sub-objects within it before assigning values. Once
     * all options are implemented with modules, they can be changed to
     * own their own storage and this parameter will no longer be required.
     *
     * If \c preprocessingStrings is provided, configures modules to write
     * preprocessing-only string parameters (like group specifications) to
     * the provided storage. These are not serialized to TPR.
     *
     * \note REQUIRED FOR GROMPP: Both \c inputrec and \c preprocessingStrings
     *       must be provided when parsing MDP files in grompp. This ensures
     *       modules can register options for preprocessing-only parameters.
     * \note REQUIRED FOR MDRUN: Only \c inputrec is needed when reading TPR
     *       files in mdrun/tools (pass nullptr for \c preprocessingStrings).
     * \note OPTIONAL FOR TESTS: Can be called without inputrec for unit tests
     *       that don't need module sub-object routing.
     *
     * \param[in]  params               Contains keys and values from user
     *     input (and defaults) to configure modules that have
     *     registered options with those keys.
     * \param[out] errorHandler          Called to report errors.
     * \param[in]  inputrec              Optional inputrec for module routing (default nullptr)
     * \param[in]  preprocessingStrings  Optional preprocessing string storage (default nullptr)
     */
    void assignOptionsToModules(const KeyValueTreeObject&  params,
                                IKeyValueTreeErrorHandler* errorHandler,
                                t_inputrec*                inputrec             = nullptr,
                                inputrec_strings*          preprocessingStrings = nullptr);

    /*! \brief
     * Normalizes inputrec parameters to match current code version.
     *
     * This orders the parameters in `ir->param` to match the current code
     * and adds any missing defaults.  It also throws an error if the
     * inputrec contains parameters that are not recognized by any module.
     *
     * If routeSubObjects is true, also configures modules to write to the
     * appropriate sub-objects within inputrec before processing.
     *
     * If preprocessingStrings is provided, also configures modules to register
     * preprocessing-only options so they are recognized during validation.
     *
     * \note REQUIRED: \c ir must not be null - this is a required parameter
     *       for all code paths.
     * \note FOR GROMPP: Pass both routeSubObjects=true and non-null
     *       preprocessingStrings when processing MDP files to ensure all
     *       parameters (including preprocessing-only) are recognized.
     * \note FOR MDRUN: Pass routeSubObjects=true but nullptr for
     *       preprocessingStrings when reading TPR files.
     *
     * \param[in] ir                  The inputrec to adjust (must not be null)
     * \param[in] routeSubObjects     Whether to route sub-objects to modules (default true)
     * \param[in] preprocessingStrings Optional preprocessing string storage (default nullptr)
     */
    void adjustInputrecBasedOnModules(t_inputrec*       ir,
                                      bool              routeSubObjects      = true,
                                      inputrec_strings* preprocessingStrings = nullptr);

    /*! \brief
     * Returns an interface for initializing and finalizing output for modules.
     */
    IMDOutputProvider* outputProvider();
    /*! \brief
     * Returns an object for computing forces from the modules.
     *
     * \param[in] wallCycle  Pointer to a wall cycle object, can be nullptr, in which case
     *                       no cycle counting will be performed for force providers.
     * \returns an object for computing forces from the modules.
     *
     * \throws APIError if a module adds a force provider for halo atoms
     *         without declaring MDModuleDescriptor::writesHaloForces, or a hook
     *         after the potential energy without declaring
     *         MDModuleDescriptor::addsForcesAfterPotentialEnergy.
     */
    ForceProviders* initForceProviders(gmx_wallcycle* wallCycle);

    /*! \brief Subscribe MDModules to notifications during pre-processing.
     *
     * Allows MDModules to subscribe to notifications that are called back
     * during pre processing an MD simulation, after the options were
     * assigned to the modules.
     */
    void subscribeToPreProcessingNotifications();

    /*! \brief Subscribe MDModules to simulation setup notifications.
     *
     * Allows MDModules to subscribe to notifications that are called back
     * during the set up of an MD simulation, after the options were
     * assigned to the modules. MDModules itself answers MDModulesEnergyTerms
     * with the energy terms the modules declare.
     *
     * \throws APIError When a module subscribes to the multi-simulation
     *     notification without declaring that it shares state across
     *     simulations, see MDModuleDescriptor.
     */
    void subscribeToSimulationSetupNotifications();

    /*! \brief Subscribe MDModules to simulation run notifications.
     *
     * Allows MDModules to subscribe to notifications that are called back
     * during an MD simulation, after the options were assigned to the modules
     * and setup notifications have been emitted.
     */
    void subscribeToSimulationRunNotifications();

    /*!
     * \brief Add a module to the container associated with the given name.
     *
     * An object may be added by a client to the bound MD Modules at run time.
     * Both the client and the MDModules object may need to extend the life
     * of the provided object. However, the MDModules container guarantees
     * to extend the life of a provided object for as long as its consumers
     * may attempt to use its the interfaces accessible through IMDModule
     * methods.
     *
     * The module is registered after all modules registered so far,
     * and is treated like a built-in module: if it provides mdp
     * options, they are placed in the section of the structured input
     * that it names with IMdpOptionProvider::mdpSectionName().
     *
     * \param name   Unique name for this module in the collection
     * \param module Implements some sort of modular functionality for MD.
     *
     * \internal
     * Adding a module at an arbitrary point in the MDModules life breaks
     * some assumptions in the protocol of the other member functions. If
     * MDModules should not change after some point, we should move this
     * to a builder class.
     *
     * \throws APIError When a module with the given name is already added
     */
    void add(std::string_view name, std::shared_ptr<IMDModule> module);

    /*! \brief Loads the plugins of the modules \p names and adds their modules.
     *
     * The plugins are looked for in \p searchPath, see
     * mdModulePluginSearchPath(). Their modules are added after all
     * modules added so far, sorted by name, so the order does not depend
     * on the order of \p names. A name given twice is loaded once.
     *
     * \throws InvalidInputError when a plugin cannot be loaded, see
     *         loadMDModulePlugin()
     * \throws APIError when a module with the name of a plugin is already added
     */
    void addPlugins(ArrayRef<const std::string> names, ArrayRef<const std::filesystem::path> searchPath);

    //! Returns the plugins loaded by addPlugins(), in the order their modules were added
    ArrayRef<const MDModulePlugin> plugins() const;

    /*! \brief Writes the name and version of each loaded plugin to \p builder
     *
     * grompp writes them to the internal parameters in the tpr, from
     * which mdrun reads the names with recordedPluginNames(). Nothing is
     * written when no plugin is loaded.
     */
    void writePluginRecord(KeyValueTreeObjectBuilder builder) const;

    //! Returns the names of the plugins that writePluginRecord() wrote to \p internalParameters
    static std::vector<std::string> recordedPluginNames(const KeyValueTreeObject& internalParameters);

    /*! \brief Returns whether a module declares that it shares state
     * across the simulations of a multi-simulation run.
     *
     * Valid after the options were assigned to the modules.
     */
    bool anyModuleSharesStateAcrossSimulations() const;

    /*! \brief Returns whether a module declares that it adds forces
     * after the potential energy is known.
     *
     * Valid after the options were assigned to the modules.
     */
    bool anyModuleAddsForcesAfterPotentialEnergy() const;

    /*! \brief Returns whether a module declares that it provides the
     * direct Coulomb interactions.
     *
     * Valid after the options were assigned to the modules.
     *
     * \throws APIError When more than one module declares it
     */
    bool anyModuleProvidesDirectCoulomb() const;

    /*! \brief Return a handle to the notifiers used for callbacks between modules.
     */
    const MDModulesNotifiers& notifiers();

private:
    class Impl;

    std::unique_ptr<Impl> impl_;
};

} // namespace gmx

#endif
