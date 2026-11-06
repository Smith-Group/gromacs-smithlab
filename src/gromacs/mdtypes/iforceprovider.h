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
 * Declares gmx::IForceProvider and ForceProviders.
 *
 * See \ref page_mdmodules for an overview of this and associated interfaces.
 *
 * \author Teemu Murtola <teemu.murtola@gmail.com>
 * \author Carsten Kutzner <ckutzne@gwdg.de>
 * \inlibraryapi
 * \ingroup module_mdtypes
 */
#ifndef GMX_MDTYPES_IFORCEPROVIDER_H
#define GMX_MDTYPES_IFORCEPROVIDER_H

#include <cstdint>

#include <memory>
#include <string>

#include "gromacs/utility/arrayref.h"
#include "gromacs/utility/gmxassert.h"
#include "gromacs/utility/real.h"
#include "gromacs/utility/vec.h"
#include "gromacs/utility/vectypes.h"

class DDBalanceRegionHandler;
struct gmx_domdec_t;
struct gmx_enerdata_t;
struct gmx_wallcycle;
struct t_forcerec;

namespace gmx
{

template<typename T>
class ArrayRef;
class ForceWithShiftForces;
class ForceWithVirial;
class MpiComm;

/*! \brief The atoms a force provider writes forces on
 *
 * Declared when the provider is added to ForceProviders.
 */
enum class ForceProviderAtoms
{
    //! Home atoms only, through ForceProviderOutput::forceWithVirial_
    Home,
    /*! \brief Home and halo atoms, through ForceProviderOutput::forceWithShiftForces_
     *
     * For pair forces between home atoms and the halo: the forces on the
     * halo atoms are communicated to their home ranks and enter the virial
     * through the shift forces, as the forces of the nonbonded kernels do.
     */
    HomeAndHalo
};

/*! \brief Where in the MD step a hook runs
 *
 * The values leave room for phases in between.
 */
enum class StepPhase : int
{
    //! Before anything reads the coordinates, e.g. making molecules whole
    PrepareGeometry = 100,
    //! Reads the coordinates to compute collective variables for the phases below
    CollectiveVariables = 200,
    //! Adds forces; where all force providers ran before phases existed
    Biases = 300,
    //! Needs the potential energy, see IPostEnergyProvider
    EnergyDependent = 400,
    //! Reads the final forces and energies and does not change them
    Observers = 500
};

/*! \brief When a hook runs within the MD step
 *
 * Hooks run by phase, then by ascending priority. Hooks with the same
 * phase and priority run in the order they were added.
 */
struct HookOrder
{
    //! The phase the hook runs in
    StepPhase phase = StepPhase::Biases;
    //! The order within the phase, ascending
    int priority = 0;
};

//! What the host knows about a step before its work is set up
struct StepInfo
{
    //! The step
    int64_t step = 0;
    //! The time at the step
    double time = 0;
    //! Whether the pair list is constructed at the step
    bool isNeighborSearchStep = false;
    //! Whether a checkpoint is written at the step
    bool isCheckpointingStep = false;
    //! Whether this is the last step
    bool isLastStep = false;
};

//! What a module needs to be computed at a step
struct StepRequest
{
    //! The potential energy, for an IPostEnergyProvider
    bool needPotentialEnergy = false;
    //! The virial
    bool needVirial = false;
};


/*! \libinternal \brief
 * Helper struct that bundles data for passing it over to the force providers
 *
 * This is a short-lived container that bundles up all necessary input data for the
 * force providers. Its only purpose is to allow calling forceProviders->calculateForces()
 * with just two arguments, one being the container for the input data,
 * the other the container for the output data.
 *
 * Both ForceProviderInput as well as ForceProviderOutput only package existing
 * data structs together for handing it over to calculateForces(). Apart from the
 * POD entries they own nothing.
 */
class ForceProviderInput
{
public:
    /*! \brief Constructor assembles all necessary force provider input data
     *
     * \param[in]  x        Atomic positions.
     * \param[in]  homenr   Number of atoms on the domain.
     * \param[in]  chargeA  Atomic charges for atoms on the domain.
     * \param[in]  massT    Atomic masses for atoms on the domain.
     * \param[in]  time     The current time in the simulation.
     * \param[in]  step     The current step in the simulation
     * \param[in]  box      The simulation box.
     * \param[in]  mpiComm  Communication object for my group.
     * \param[in]  dd       Domain decomposition object, pass nullptr when DD is not in use.
     */
    ForceProviderInput(ArrayRef<const RVec> x,
                       int                  homenr,
                       ArrayRef<const real> chargeA,
                       ArrayRef<const real> massT,
                       double               time,
                       int64_t              step,
                       const matrix         box,
                       const MpiComm&       mpiComm,
                       const gmx_domdec_t*  dd) :
        x_(x), homenr_(homenr), chargeA_(chargeA), massT_(massT), t_(time), step_(step), mpiComm_(mpiComm), dd_(dd)
    {
        copy_mat(box, box_);
    }

    ArrayRef<const RVec> x_; //!< The atomic positions
    int                  homenr_;
    ArrayRef<const real> chargeA_;
    ArrayRef<const real> massT_;
    double               t_;    //!< The current time in the simulation
    int64_t              step_; //!< The current step in the simulation
    matrix               box_ = { { 0, 0, 0 }, { 0, 0, 0 }, { 0, 0, 0 } }; //!< The simulation box
    const MpiComm&       mpiComm_; //!< Communication object for my group
    const gmx_domdec_t*  dd_;      //!< Domain decomposition object, deprecated
    /*! \brief The number of local atoms including the halo
     *
     * \p x_ holds them: the home atoms first, then the halo.
     * Without domain decomposition, this is \p homenr_.
     * With direct halo communication, the home and the halo atoms include
     * filler particles, whose global atom index in \p dd_ is -1.
     */
    int numAtomsWithHalo_ = 0;
    /*! \brief The dynamic load balancing region handler, can be nullptr
     *
     * The balancing region is closed when the force providers are called.
     * A provider that writes halo forces can add the cycles of its
     * computation to the load with
     * DDBalanceRegionHandler::addCyclesAfterForceComputationCpu().
     */
    const DDBalanceRegionHandler* ddBalanceRegionHandler_ = nullptr;
};

/*! \brief Take pointer, check if valid, return reference
 */
template<class T>
T& makeRefFromPointer(T* ptr)
{
    GMX_ASSERT(ptr != nullptr, "got null pointer");
    return *ptr;
}

/*! \libinternal \brief
 * Helper struct bundling the output data of a force provider
 *
 * Same as for the ForceProviderInput class, but these variables can be written as well.
 */
class ForceProviderOutput
{
public:
    /*! \brief Constructor assembles all necessary force provider output data
     *
     * \param[in,out]  forceWithVirial  Container for force and virial
     * \param[in,out]  enerd            Structure containing energy data
     */
    ForceProviderOutput(ForceWithVirial* forceWithVirial, gmx_enerdata_t* enerd) :
        forceWithVirial_(makeRefFromPointer(forceWithVirial)), enerd_(makeRefFromPointer(enerd))
    {
    }

    ForceWithVirial& forceWithVirial_; //!< Container for force and virial
    gmx_enerdata_t&  enerd_;           //!< Structure containing energy data
    /*! \brief The forces on home and halo atoms, and the shift forces
     *
     * Only set for a provider added with ForceProviderAtoms::HomeAndHalo,
     * nullptr for the others.
     */
    ForceWithShiftForces* forceWithShiftForces_ = nullptr;
};


/*! \libinternal \brief
 * Interface for a component that provides forces during MD.
 *
 * Modules implementing IMDModule generally implement this internally, and use
 * IMDModule::initForceProviders() to register their implementation in
 * ForceProviders.
 *
 * The interface most likely requires additional generalization for use in
 * other modules than the current electric field implementation.
 *
 * The forces that are produced by force providers are not taken into account
 * in the calculation of the virial. When applicable, the provider should
 * compute its own virial contribution. Forces written to
 * ForceProviderOutput::forceWithShiftForces_ (ForceProviderAtoms::HomeAndHalo)
 * do enter the virial: a pair across a periodic boundary in a dimension without
 * domain decomposition adds its shift force, as the bonded forces do.
 *
 * The providers are called on every PP rank at every step that computes
 * forces, so a provider may communicate over ForceProviderInput::mpiComm_.
 *
 * \inlibraryapi
 * \ingroup module_mdtypes
 */
class IForceProvider
{
public:
    /*! \brief
     * Computes forces.
     *
     * \param[in]    forceProviderInput    struct that collects input data for the force providers
     * \param[in,out] forceProviderOutput   struct that collects output data of the force providers
     */
    virtual void calculateForces(const ForceProviderInput& forceProviderInput,
                                 ForceProviderOutput*      forceProviderOutput) = 0;

protected:
    ~IForceProvider() {}
};

/*! \libinternal \brief
 * Interface for a component that prepares for each MD step.
 *
 * Called on every PP rank before the work of the step is set up, which is
 * the last point at which the module can still ask for the potential
 * energy or the virial at that step. The MD integrators of the legacy
 * simulator and rerun call it; the modular simulator, energy minimization
 * and the other integrators do not.
 *
 * \inlibraryapi
 * \ingroup module_mdtypes
 */
class IStepPreparer
{
public:
    /*! \brief
     * Prepares for a step.
     *
     * \param[in]  stepInfo  What the host knows about the step
     * \param[out] request   What the module needs at the step, all false on entry
     */
    virtual void prepareStep(const StepInfo& stepInfo, StepRequest* request) = 0;

protected:
    ~IStepPreparer() = default;
};

/*! \libinternal \brief
 * Input of a hook that runs after the potential energy is known
 */
struct PostEnergyInput
{
    //! The step
    int64_t step_;
    //! Whether \p potentialEnergy_ and \p virial_ were computed at this step
    bool energyWasComputed_;
    //! The potential energy of this rank, not summed over the ranks
    real potentialEnergy_;
    //! The force virial of this rank, before the contributions of the hooks
    const tensor& virial_;
};

/*! \libinternal \brief
 * Interface for a component that computes forces from the potential energy.
 *
 * Called on every PP rank at every force computation, after the potential
 * energy of the rank is accumulated and the force on the home atoms is
 * complete, e.g. to bias the potential energy. ForceProviderOutput::forceWithVirial_
 * holds that complete force (with multiple time stepping, the combined
 * force); the hook adds its virial contribution with
 * ForceWithVirial::addVirialContribution(). Forces on virtual sites are
 * spread after the hooks.
 *
 * At steps that compute no energy, \p energyWasComputed_ is false; a
 * module that needs the energy asks for it with IStepPreparer.
 *
 * \inlibraryapi
 * \ingroup module_mdtypes
 */
class IPostEnergyProvider
{
public:
    /*! \brief
     * Computes forces from the potential energy.
     *
     * \param[in]     input   The potential energy and the virial
     * \param[in,out] output  The complete force on the home atoms, and the energies
     */
    virtual void applyAfterPotentialEnergy(const PostEnergyInput& input, ForceProviderOutput* output) = 0;

protected:
    ~IPostEnergyProvider() = default;
};

/*! \libinternal \brief
 * Evaluates forces from a collection of gmx::IForceProvider.
 *
 * \inlibraryapi
 * \ingroup module_mdtypes
 */
class ForceProviders
{
public:
    /*! \brief
     * Constructor.
     *
     * \param[in] wallCycle  Pointer to a wallcycle counter struct, can be nullptr
     */
    ForceProviders(gmx_wallcycle* wallCycle = nullptr);
    ~ForceProviders();

    /*! \brief
     * Adds a provider.
     *
     * \param[in] provider          The force provider callback function
     * \param[in] cycleCounterName  A non-empty string will add a cycle counter with the given name
     *                              that registers the time spent in the force provider function
     * \param[in] atoms             The atoms the provider writes forces on
     */
    void addForceProvider(gmx::IForceProvider* provider,
                          const std::string&   cycleCounterName = "",
                          ForceProviderAtoms   atoms            = ForceProviderAtoms::Home);

    /*! \brief
     * Adds a provider that runs in the given order.
     *
     * \param[in] provider          The force provider callback function
     * \param[in] order             When the provider runs, in a phase up to StepPhase::Biases
     * \param[in] cycleCounterName  A non-empty string will add a cycle counter with the given name
     *                              that registers the time spent in the force provider function
     * \param[in] atoms             The atoms the provider writes forces on
     *
     * \throws APIError if the phase is after StepPhase::Biases
     */
    void addForceProvider(gmx::IForceProvider* provider,
                          HookOrder            order,
                          const std::string&   cycleCounterName = "",
                          ForceProviderAtoms   atoms            = ForceProviderAtoms::Home);

    //! Adds a component that prepares for each step; all are called, in the order added
    void addStepPreparer(IStepPreparer* preparer);

    /*! \brief
     * Adds a hook that runs after the potential energy is known.
     *
     * \param[in] provider          The hook
     * \param[in] order             When the hook runs, in StepPhase::EnergyDependent or later
     * \param[in] cycleCounterName  A non-empty string will add a cycle counter with the given name
     *
     * \throws APIError if the phase is before StepPhase::EnergyDependent
     */
    void addPostEnergyProvider(IPostEnergyProvider* provider,
                               HookOrder            order,
                               const std::string&   cycleCounterName = "");

    //! Whether there are modules added.
    bool hasForceProvider() const;

    //! The number of providers that write forces on halo atoms
    int numHaloForceProviders() const;

    //! Whether there are hooks that run after the potential energy is known
    bool hasPostEnergyProvider() const;

    /*! \brief
     * Prepares all components for a step.
     *
     * \param[in]  stepInfo  What the host knows about the step
     * \param[out] request   What any of the components needs at the step
     */
    void prepareStep(const StepInfo& stepInfo, StepRequest* request) const;

    //! Computes forces, by phase and priority of the providers.
    void calculateForces(const gmx::ForceProviderInput& forceProviderInput,
                         gmx::ForceProviderOutput*      forceProviderOutput) const;

    //! Runs the hooks that need the potential energy, by phase and priority.
    void applyAfterPotentialEnergy(const PostEnergyInput& input, ForceProviderOutput* output) const;

private:
    class Impl;

    std::unique_ptr<Impl> impl_;
};

} // namespace gmx

#endif
