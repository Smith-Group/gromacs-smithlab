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
 * Declares the C boundary between the host and an MDModules plugin.
 *
 * The structs here have C layout, and fields are only ever appended to
 * them, so a plugin and a host built at different times agree on the
 * fields both know.
 *
 * \author Amr Alhossary <amr_alhossary@hotmail.com>
 * \inlibraryapi
 * \ingroup module_mdtypes
 */
#ifndef GMX_MDTYPES_MDMODULEPLUGIN_H
#define GMX_MDTYPES_MDMODULEPLUGIN_H

#include <cstdint>

//! Value of GmxHostApiInfo::mpi for a host without MPI
constexpr uint32_t c_gmxHostMpiNone = 0;
//! Value of GmxHostApiInfo::mpi for a host with thread-MPI
constexpr uint32_t c_gmxHostMpiThread = 1;
//! Value of GmxHostApiInfo::mpi for a host with an MPI library
constexpr uint32_t c_gmxHostMpiLibrary = 2;

/*! \libinternal \brief
 * How the host was built, as far as a plugin needs to know.
 *
 * The host fills it before it asks a plugin anything. A plugin built
 * for an older host reads only the fields within its own
 * sizeof(GmxHostApiInfo), which \c structSize lets it check.
 */
struct GmxHostApiInfo
{
    //! sizeof(GmxHostApiInfo) in the host
    uint32_t structSize;
    //! sizeof(real) in the host: 4 in mixed precision, 8 in double precision
    uint32_t sizeofReal;
    //! The major version of the interface (GMX_MDMODULE_API_MAJOR)
    uint64_t apiMajor;
    //! The hash of the interface headers (GMX_MDMODULE_API_HASH)
    uint64_t apiHash;
    //! The \Gromacs version, for messages
    const char* gromacsVersion;
    //! The MPI of the host: c_gmxHostMpiNone, c_gmxHostMpiThread or c_gmxHostMpiLibrary
    uint32_t mpi;
};

/*! \libinternal \brief
 * Who a plugin is.
 *
 * Returned by the plugin's gmxPluginInfo() entry point, before the
 * module is created. \c name is the name of the module; mdp files name
 * the plugin by it, and the plugin file is named after it.
 */
struct GmxPluginInfo
{
    //! sizeof(GmxPluginInfo) in the plugin
    uint32_t structSize;
    //! The name of the module
    const char* name;
    //! The version of the plugin, for messages and the tpr record
    const char* version;
};

namespace gmx
{

class IMDModule;

//! Returns how this host was built
GmxHostApiInfo hostApiInfo();

} // namespace gmx

extern "C"
{
    /*! \brief Entry point of a plugin: whether it can run in \p host
     *
     * When it cannot, it sets \p reason to a message that it owns.
     */
    using GmxPluginIsCompatibleFunction = bool (*)(const GmxHostApiInfo* host, const char** reason);
    //! Entry point of a plugin: who it is
    using GmxPluginInfoFunction = const GmxPluginInfo* (*)();
    //! Entry point of a plugin: creates its module
    using GmxPluginCreateFunction = gmx::IMDModule* (*)();
    //! Entry point of a plugin: destroys a module that gmxPluginCreate() returned
    using GmxPluginDestroyFunction = void (*)(gmx::IMDModule* module);
}

//! The names of the entry points of a plugin
constexpr const char* c_gmxPluginIsCompatibleName = "gmxPluginIsCompatible";
//! \copydoc c_gmxPluginIsCompatibleName
constexpr const char* c_gmxPluginInfoName = "gmxPluginInfo";
//! \copydoc c_gmxPluginIsCompatibleName
constexpr const char* c_gmxPluginCreateName = "gmxPluginCreate";
//! \copydoc c_gmxPluginIsCompatibleName
constexpr const char* c_gmxPluginDestroyName = "gmxPluginDestroy";

#endif
