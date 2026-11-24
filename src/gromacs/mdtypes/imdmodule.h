/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright 2017- The GROMACS Authors
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
 * Declares gmx::IMDModule.
 *
 * See \ref page_mdmodules for an overview of this and associated interfaces.
 *
 * \inlibraryapi
 * \ingroup module_mdtypes
 */
#ifndef GMX_MDTYPES_IMDMODULE_H
#define GMX_MDTYPES_IMDMODULE_H

#include <vector>

enum class InteractionFunction : int;

namespace gmx
{

class ForceProviders;
class IMDOutputProvider;
class IMdpOptionProvider;
struct MDModulesNotifiers;

/*! \libinternal \brief
 * What a module declares about itself for this run.
 *
 * The host reads the declarations after the mdp options have been
 * assigned to the modules, so a module can derive them from its
 * options. Fields are only ever appended to this struct, and the host
 * allocates it, so a module sets only the fields it knows.
 *
 * \inlibraryapi
 * \ingroup module_mdtypes
 */
struct MDModuleDescriptor
{
    /*! \brief Whether the module couples the simulations of a multi-simulation run.
     *
     * The simulations then act in lockstep: checkpoints and stop
     * conditions are communicated between them. A module that
     * subscribes to the gmx_multisim_t notification has to declare this.
     */
    bool sharesStateAcrossSimulations = false;
    /*! \brief Whether the module computes the direct (short-range) Coulomb
     * interactions instead of the nonbonded kernels. At most one module may.
     */
    bool providesDirectCoulomb = false;
    /*! \brief Whether the module writes forces on halo atoms
     *
     * A module that adds a force provider with ForceProviderAtoms::HomeAndHalo
     * has to declare this.
     */
    bool writesHaloForces = false;
    /*! \brief Whether the module adds forces after the potential energy is known
     *
     * A module that adds an IPostEnergyProvider has to declare this. Such
     * forces are added on the CPU, so the update then runs on the CPU.
     */
    bool addsForcesAfterPotentialEnergy = false;
    //! The energy terms the module writes, which the host writes to the energy file
    std::vector<InteractionFunction> energyTerms;
};

/*! \libinternal \brief
 * Extension module for \Gromacs simulations.
 *
 * The methods that return other interfaces return null for those
 * interfaces that the module does not implement.
 *
 *  See also:
 * - imdpprovider_test_helper.h: Provides a helper to populate IMdpOptionProviders from key-value trees for testing.
 * - imdpoptionprovider_helpers.h for helpers to transform and output mdp options.
 *
 * \inlibraryapi
 * \ingroup module_mdtypes
 */
class IMDModule
{
public:
    virtual ~IMDModule() {}

    /*! \brief Fills in what the module declares about itself for this run.
     *
     * Called after the mdp options have been assigned to the module.
     * The host passes a descriptor that declares nothing, and the module
     * sets what it declares; the default sets nothing. A module built
     * against an older header therefore leaves the fields it does not
     * know at their defaults.
     */
    virtual void describe(MDModuleDescriptor* /* descriptor */) const {}
    //! Returns an interface for handling mdp input (and tpr I/O).
    virtual IMdpOptionProvider* mdpOptionProvider() = 0;
    //! Returns an interface for handling output files during simulation.
    virtual IMDOutputProvider* outputProvider() = 0;
    //! Initializes force providers from this module.
    virtual void initForceProviders(ForceProviders* forceProviders) = 0;
    //! Subscribe to pre processing notifications
    virtual void subscribeToPreProcessingNotifications(MDModulesNotifiers* notifiers) = 0;
    //! Subscribe to simulation setup notifications
    virtual void subscribeToSimulationSetupNotifications(MDModulesNotifiers* notifiers) = 0;
    //! Subscribe to simulation run notifications
    virtual void subscribeToSimulationRunNotifications(MDModulesNotifiers* notifiers) = 0;
};

} // namespace gmx

#endif
