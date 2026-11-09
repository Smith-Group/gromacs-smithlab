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
 * The electric field built as an MDModules plugin from its own source,
 * for the tests of plugins in mdrun
 *
 * The plugin is called electric-field-plugin, and so are its mdp
 * options, so that it can be loaded next to the built-in module.
 *
 * \author Amr Alhossary <amr_alhossary@hotmail.com>
 * \ingroup module_mdrun_integration_tests
 */
#include "gmxpre.h"

#include <memory>

#include "gromacs/applied_forces/electricfield.h"
#include "gromacs/mdtypes/imdmodule.h"
#include "gromacs/mdtypes/mdmodulepluginentrypoints.h"

namespace gmx
{
namespace test
{
namespace
{

//! Creates the module of the plugin
std::unique_ptr<IMDModule> createElectricFieldPluginModule()
{
    return ElectricFieldModuleInfo::create("electric-field-plugin");
}

} // namespace
} // namespace test
} // namespace gmx

GMX_MDMODULE_PLUGIN("electric-field-plugin", "1.0", ::gmx::test::createElectricFieldPluginModule)
