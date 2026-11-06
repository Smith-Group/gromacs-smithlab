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
/*! \internal \file
 * \brief
 * Implements classes from iforceprovider.h.
 *
 * \author Teemu Murtola <teemu.murtola@gmail.com>
 * \ingroup module_mdtypes
 */
#include "gmxpre.h"

#include "iforceprovider.h"

#include <algorithm>
#include <tuple>
#include <utility>
#include <vector>

#include "gromacs/timing/wallcycle.h"
#include "gromacs/utility/exceptions.h"

using namespace gmx;

namespace
{

//! Whether a hook in order \p a runs before one in order \p b
bool runsBefore(const HookOrder& a, const HookOrder& b)
{
    return std::tie(a.phase, a.priority) < std::tie(b.phase, b.priority);
}

//! Inserts \p entry after the entries that run before it or together with it
template<typename Entry>
void insertInOrder(std::vector<Entry>* entries, Entry entry)
{
    const auto position = std::upper_bound(entries->begin(),
                                           entries->end(),
                                           entry,
                                           [](const Entry& a, const Entry& b)
                                           { return runsBefore(a.order, b.order); });
    entries->insert(position, std::move(entry));
}

} // namespace

class ForceProviders::Impl
{
public:
    Impl(gmx_wallcycle* wallCycle) : wallCycle_(wallCycle) {}

    //! A provider with its order, its cycle counter and the atoms it writes forces on
    struct Entry
    {
        IForceProvider*                 provider;
        HookOrder                       order;
        std::optional<WallCycleCounter> wallCycleCounter;
        ForceProviderAtoms              atoms;
    };

    //! A hook after the potential energy, with its order and its cycle counter
    struct PostEnergyEntry
    {
        IPostEnergyProvider*            provider;
        HookOrder                       order;
        std::optional<WallCycleCounter> wallCycleCounter;
    };

    //! Calls \p call within the cycle counter \p counter, if there is one
    template<typename Call>
    void timed(const std::optional<WallCycleCounter>& counter, Call call) const
    {
        if (counter.has_value())
        {
            wallcycle_start(wallCycle_, counter.value());
        }
        call();
        if (counter.has_value())
        {
            wallcycle_stop(wallCycle_, counter.value());
        }
    }

    //! The force providers, in the order they run
    std::vector<Entry> providers_;
    //! The components that prepare for each step
    std::vector<IStepPreparer*> stepPreparers_;
    //! The hooks after the potential energy, in the order they run
    std::vector<PostEnergyEntry> postEnergyProviders_;
    gmx_wallcycle*               wallCycle_;
};

ForceProviders::ForceProviders(gmx_wallcycle* wallCycle) : impl_(new Impl(wallCycle)) {}

ForceProviders::~ForceProviders() {}

void ForceProviders::addForceProvider(gmx::IForceProvider* provider,
                                      const std::string&   cycleCounterName,
                                      ForceProviderAtoms   atoms)
{
    addForceProvider(provider, HookOrder{}, cycleCounterName, atoms);
}

void ForceProviders::addForceProvider(gmx::IForceProvider* provider,
                                      HookOrder            order,
                                      const std::string&   cycleCounterName,
                                      ForceProviderAtoms   atoms)
{
    if (order.phase > StepPhase::Biases)
    {
        GMX_THROW(APIError(
                "A force provider runs before the potential energy is known, in phase Biases "
                "or earlier; use an IPostEnergyProvider for a later phase"));
    }
    std::optional<WallCycleCounter> counter;
    if (impl_->wallCycle_)
    {
        counter = impl_->wallCycle_->registerCycleCounter(cycleCounterName);
    }

    insertInOrder(&impl_->providers_, Impl::Entry{ provider, order, counter, atoms });
}

void ForceProviders::addStepPreparer(IStepPreparer* preparer)
{
    impl_->stepPreparers_.push_back(preparer);
}

void ForceProviders::addPostEnergyProvider(IPostEnergyProvider* provider,
                                           HookOrder            order,
                                           const std::string&   cycleCounterName)
{
    if (order.phase < StepPhase::EnergyDependent)
    {
        GMX_THROW(APIError(
                "A hook after the potential energy runs in phase EnergyDependent or later; use "
                "a force provider for an earlier phase"));
    }
    std::optional<WallCycleCounter> counter;
    if (impl_->wallCycle_ && !cycleCounterName.empty())
    {
        counter = impl_->wallCycle_->registerCycleCounter(cycleCounterName);
    }

    insertInOrder(&impl_->postEnergyProviders_, Impl::PostEnergyEntry{ provider, order, counter });
}

bool ForceProviders::hasForceProvider() const
{
    return !impl_->providers_.empty();
}

int ForceProviders::numHaloForceProviders() const
{
    return std::count_if(impl_->providers_.begin(),
                         impl_->providers_.end(),
                         [](const Impl::Entry& entry)
                         { return entry.atoms == ForceProviderAtoms::HomeAndHalo; });
}

bool ForceProviders::hasPostEnergyProvider() const
{
    return !impl_->postEnergyProviders_.empty();
}

void ForceProviders::prepareStep(const StepInfo& stepInfo, StepRequest* request) const
{
    for (IStepPreparer* preparer : impl_->stepPreparers_)
    {
        StepRequest preparerRequest;
        preparer->prepareStep(stepInfo, &preparerRequest);
        request->needPotentialEnergy = request->needPotentialEnergy || preparerRequest.needPotentialEnergy;
        request->needVirial = request->needVirial || preparerRequest.needVirial;
    }
}

void ForceProviders::calculateForces(const ForceProviderInput& forceProviderInput,
                                     ForceProviderOutput*      forceProviderOutput) const
{
    // The providers of home-atom forces do not get the halo force buffer
    ForceProviderOutput homeAtomsOutput   = *forceProviderOutput;
    homeAtomsOutput.forceWithShiftForces_ = nullptr;

    for (const auto& entry : impl_->providers_)
    {
        impl_->timed(entry.wallCycleCounter,
                     [&]()
                     {
                         entry.provider->calculateForces(forceProviderInput,
                                                         entry.atoms == ForceProviderAtoms::HomeAndHalo
                                                                 ? forceProviderOutput
                                                                 : &homeAtomsOutput);
                     });
    }
}

void ForceProviders::applyAfterPotentialEnergy(const PostEnergyInput& input, ForceProviderOutput* output) const
{
    for (const auto& entry : impl_->postEnergyProviders_)
    {
        impl_->timed(entry.wallCycleCounter,
                     [&]() { entry.provider->applyAfterPotentialEnergy(input, output); });
    }
}
