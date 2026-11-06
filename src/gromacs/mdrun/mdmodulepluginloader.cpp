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
 * Defines the loader of MDModules plugins.
 *
 * \author Amr Alhossary <amr_alhossary@hotmail.com>
 * \ingroup module_mdrun
 */
#include "gmxpre.h"

#include "mdmodulepluginloader.h"

#include "config.h"

#include <cstddef>
#include <cstdlib>

#include <algorithm>
#include <string>

#if GMX_MDMODULE_PLUGINS
#    include <dlfcn.h>
#endif

#include "gromacs/mdtypes/imdmodule.h"
#include "gromacs/mdtypes/mdmoduleplugin.h"
#include "gromacs/utility/exceptions.h"
#include "gromacs/utility/programcontext.h"
#include "gromacs/utility/stringutil.h"

#include "buildinfo.h"

namespace gmx
{

std::vector<std::filesystem::path> mdModulePluginSearchPath()
{
    std::vector<std::filesystem::path> searchPath;
    if (const char* pathVariable = std::getenv("GMX_MDMODULE_PLUGIN_PATH"))
    {
        for (const std::string& directory : splitDelimitedString(pathVariable, ':'))
        {
            if (!directory.empty())
            {
                searchPath.emplace_back(directory);
            }
        }
        return searchPath;
    }
    const InstallationPrefixInfo installationPrefix = getProgramContext().installationPrefix();
    if (!installationPrefix.path_.empty() && !installationPrefix.sourceLayoutTreeLike_)
    {
        searchPath.push_back(installationPrefix.path_ / GMX_INSTALL_MDMODULE_PLUGINDIR);
    }
    return searchPath;
}

namespace
{

//! Describes \p host for messages
std::string describeHost(const GmxHostApiInfo& host)
{
    const char* mpi = host.mpi == c_gmxHostMpiLibrary
                              ? "an MPI library"
                              : (host.mpi == c_gmxHostMpiThread ? "thread-MPI" : "no MPI");
    return formatString(
            "%s, MDModules plugin interface version %llu (hash %016llx), %s precision, %s",
            host.gromacsVersion,
            static_cast<unsigned long long>(host.apiMajor),
            static_cast<unsigned long long>(host.apiHash),
            host.sizeofReal == sizeof(double) ? "double" : "mixed",
            mpi);
}

//! Throws that plugin \p name cannot be loaded because of \p reason
[[noreturn]] void throwCannotLoad(const std::string& name, const std::string& reason)
{
    GMX_THROW(InvalidInputError(formatString(
            "Cannot load the MDModules plugin '%s': %s", name.c_str(), reason.c_str())));
}

} // namespace

MDModulePlugin loadMDModulePlugin(const std::string& name, ArrayRef<const std::filesystem::path> searchPath)
{
    if (name.empty() || name.find_first_of("/\\") != std::string::npos || name == "." || name == "..")
    {
        throwCannotLoad(name, "a plugin is named by its module name, not by a path.");
    }
#if GMX_MDMODULE_PLUGINS
    const std::string fileName = name + GMX_MDMODULE_PLUGIN_SUFFIX;
    const auto        directory =
            std::find_if(searchPath.begin(),
                         searchPath.end(),
                         [&fileName](const std::filesystem::path& candidate)
                         { return std::filesystem::is_regular_file(candidate / fileName); });
    if (directory == searchPath.end())
    {
        std::vector<std::string> directories;
        for (const auto& candidate : searchPath)
        {
            directories.push_back(candidate.string());
        }
        throwCannotLoad(
                name,
                formatString("%s is not in any of the plugin directories (%s). Set "
                             "GMX_MDMODULE_PLUGIN_PATH to the directories that contain it.",
                             fileName.c_str(),
                             directories.empty() ? "none" : joinStrings(directories, ", ").c_str()));
    }
    const std::filesystem::path path = *directory / fileName;

    // The library is never closed: callbacks of its module can outlive the module.
    void* library = dlopen(path.c_str(), RTLD_NOW | RTLD_LOCAL);
    if (library == nullptr)
    {
        throwCannotLoad(name, formatString("%s: %s", path.c_str(), dlerror()));
    }
    const auto entryPoint = [&](const char* entryPointName)
    {
        void* symbol = dlsym(library, entryPointName);
        if (symbol == nullptr)
        {
            throwCannotLoad(name,
                            formatString("%s is not an MDModules plugin: it has no entry point %s.",
                                         path.c_str(),
                                         entryPointName));
        }
        return symbol;
    };
    // Converting the symbols to function pointers is how dlsym() is used
    auto isCompatible =
            reinterpret_cast<GmxPluginIsCompatibleFunction>(entryPoint(c_gmxPluginIsCompatibleName));
    auto pluginInfo = reinterpret_cast<GmxPluginInfoFunction>(entryPoint(c_gmxPluginInfoName));
    auto create     = reinterpret_cast<GmxPluginCreateFunction>(entryPoint(c_gmxPluginCreateName));
    auto destroy = reinterpret_cast<GmxPluginDestroyFunction>(entryPoint(c_gmxPluginDestroyName));

    const GmxHostApiInfo host   = hostApiInfo();
    const char*          reason = nullptr;
    if (!isCompatible(&host, &reason))
    {
        throwCannotLoad(name,
                        formatString("%s cannot run in this GROMACS: %s. This GROMACS: %s.",
                                     path.c_str(),
                                     reason != nullptr ? reason : "no reason given",
                                     describeHost(host).c_str()));
    }
    // A plugin built against an older header fills a shorter struct; only the
    // fields within its structSize are read.
    const GmxPluginInfo* info = pluginInfo();
    if (info == nullptr || info->structSize < offsetof(GmxPluginInfo, version) + sizeof(info->version)
        || info->name == nullptr || info->version == nullptr)
    {
        throwCannotLoad(name, formatString("%s does not say what it is.", path.c_str()));
    }
    if (name != info->name)
    {
        throwCannotLoad(name, formatString("%s is the plugin of module '%s'.", path.c_str(), info->name));
    }
    const bool saysApiVersion =
            info->structSize >= offsetof(GmxPluginInfo, apiHash) + sizeof(info->apiHash);
    IMDModule* module = create();
    if (module == nullptr)
    {
        throwCannotLoad(name, formatString("%s did not create its module.", path.c_str()));
    }
    return { name,
             info->version,
             path,
             saysApiVersion ? info->apiMajor : host.apiMajor,
             saysApiVersion ? info->apiHash : 0,
             std::shared_ptr<IMDModule>(module, destroy) };
#else
    GMX_UNUSED_VALUE(searchPath);
    throwCannotLoad(name,
                    "this GROMACS was built without support for MDModules plugins. Build it "
                    "with -DGMX_MDMODULE_PLUGINS=ON.");
#endif
}

} // namespace gmx
