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
/*! \libinternal \file
 * \brief
 * Declares the loader of MDModules plugins.
 *
 * \author Amr Alhossary <amr_alhossary@hotmail.com>
 * \inlibraryapi
 * \ingroup module_mdrun
 */
#ifndef GMX_MDRUN_MDMODULEPLUGINLOADER_H
#define GMX_MDRUN_MDMODULEPLUGINLOADER_H

#include <filesystem>
#include <memory>
#include <string>
#include <vector>

#include "gromacs/utility/arrayref.h"

namespace gmx
{

class IMDModule;

/*! \libinternal \brief
 * A module created by a plugin, and what the plugin said about itself.
 */
struct MDModulePlugin
{
    //! The name of the module
    std::string name;
    //! The version of the plugin
    std::string version;
    //! The file the plugin was loaded from
    std::filesystem::path path;
    //! The module, destroyed by the plugin
    std::shared_ptr<IMDModule> module;
};

/*! \brief Returns the directories in which plugins are looked for
 *
 * These are the directories in the environment variable
 * GMX_MDMODULE_PLUGIN_PATH, separated by colons, or, when it is not set,
 * the plugin directory of the installation.
 */
std::vector<std::filesystem::path> mdModulePluginSearchPath();

/*! \brief Loads the plugin of module \p name and creates its module
 *
 * The plugin is the file named \p name with the shared module suffix
 * in the first directory of \p searchPath that has it. It is loaded only
 * when it says it can run in this host, and when it calls itself
 * \p name. A plugin stays loaded until the program exits, as code in it
 * can still be referenced after its module is destroyed.
 *
 * \throws InvalidInputError when the plugin is not found, cannot be
 *         loaded, is not a plugin, cannot run in this host, or has
 *         another name; or when \Gromacs was built without plugin support.
 */
MDModulePlugin loadMDModulePlugin(const std::string& name, ArrayRef<const std::filesystem::path> searchPath);

} // namespace gmx

#endif
