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
 * Defines the entry points of an MDModules plugin.
 *
 * A plugin includes this header in one source file and uses
 * GMX_MDMODULE_PLUGIN once:
 * \code
 * GMX_MDMODULE_PLUGIN("my-module", "1.0", MyModuleInfo::create)
 * \endcode
 * where \c MyModuleInfo::create returns a std::unique_ptr<gmx::IMDModule>.
 *
 * \author Amr Alhossary <amr_alhossary@hotmail.com>
 * \inlibraryapi
 * \ingroup module_mdtypes
 */
#ifndef GMX_MDTYPES_MDMODULEPLUGINENTRYPOINTS_H
#define GMX_MDTYPES_MDMODULEPLUGINENTRYPOINTS_H

#include <cstddef>

#include "gromacs/mdtypes/imdmodule.h"
#include "gromacs/mdtypes/mdmoduleapiversion.h"
#include "gromacs/mdtypes/mdmoduleplugin.h"
#include "gromacs/utility/real.h"

namespace gmx
{

/*! \brief Whether a plugin built with these headers can run in \p host
 *
 * The host has to have the same major version of the interface, the
 * same precision and the same MPI as the \Gromacs that the plugin was
 * built against. When it does not, \p reason is set to why.
 */
inline bool isHostCompatibleWithPluginBuild(const GmxHostApiInfo& host, const char** reason)
{
    if (host.structSize < offsetof(GmxHostApiInfo, mpi) + sizeof(host.mpi))
    {
        *reason = "the host is older than the plugin";
        return false;
    }
    if (host.apiMajor != GMX_MDMODULE_API_MAJOR)
    {
        *reason = "the host has another major version of the MDModules plugin interface";
        return false;
    }
    if (host.sizeofReal != sizeof(real))
    {
        *reason = "the host was built in another precision";
        return false;
    }
    if (host.mpi != GMX_MDMODULE_API_MPI)
    {
        *reason = "the host was built with another MPI";
        return false;
    }
    return true;
}

} // namespace gmx

//! Gives a plugin entry point C linkage and makes it visible to the host
#define GMX_MDMODULE_PLUGIN_EXPORT extern "C" __attribute__((visibility("default")))

/*! \brief Defines the four entry points of a plugin
 *
 * \param[in] moduleName     The name of the module, a string literal
 * \param[in] pluginVersion  The version of the plugin, a string literal
 * \param[in] createModule   A function that returns a std::unique_ptr<gmx::IMDModule>
 */
#define GMX_MDMODULE_PLUGIN(moduleName, pluginVersion, createModule)                                       \
    GMX_MDMODULE_PLUGIN_EXPORT bool gmxPluginIsCompatible(const GmxHostApiInfo* host, const char** reason) \
    {                                                                                                      \
        return ::gmx::isHostCompatibleWithPluginBuild(*host, reason);                                      \
    }                                                                                                      \
    GMX_MDMODULE_PLUGIN_EXPORT const GmxPluginInfo* gmxPluginInfo()                                        \
    {                                                                                                      \
        static const GmxPluginInfo info = { sizeof(GmxPluginInfo), moduleName, pluginVersion };            \
        return &info;                                                                                      \
    }                                                                                                      \
    GMX_MDMODULE_PLUGIN_EXPORT ::gmx::IMDModule* gmxPluginCreate()                                         \
    {                                                                                                      \
        return (createModule)().release();                                                                 \
    }                                                                                                      \
    GMX_MDMODULE_PLUGIN_EXPORT void gmxPluginDestroy(::gmx::IMDModule* module)                             \
    {                                                                                                      \
        delete module;                                                                                     \
    }

#endif
