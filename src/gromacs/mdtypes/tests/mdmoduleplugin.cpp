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
 * Tests the host side of the MDModules plugin boundary
 *
 * \author Amr Alhossary <amr_alhossary@hotmail.com>
 * \ingroup module_mdtypes
 */
#include "gmxpre.h"

#include "gromacs/mdtypes/mdmoduleplugin.h"

#include "config.h"

#include <cstddef>

#include <string>
#include <utility>

#include <gmock/gmock.h>
#include <gtest/gtest.h>

#include "gromacs/mdtypes/mdmoduleapiversion.h"
#include "gromacs/mdtypes/mdmodulepluginentrypoints.h"
#include "gromacs/utility/baseversion.h"
#include "gromacs/utility/real.h"

namespace gmx
{
namespace test
{
namespace
{

TEST(HostApiInfoTest, DescribesThisBuild)
{
    const GmxHostApiInfo info = hostApiInfo();
    EXPECT_EQ(info.structSize, sizeof(GmxHostApiInfo));
    EXPECT_EQ(info.sizeofReal, sizeof(real));
    EXPECT_EQ(info.apiMajor, GMX_MDMODULE_API_MAJOR);
    EXPECT_EQ(info.apiHash, GMX_MDMODULE_API_HASH);
    EXPECT_EQ(std::string(info.gromacsVersion), std::string(gmx_version()));
    const uint32_t expectedMpi = GMX_LIB_MPI ? c_gmxHostMpiLibrary
                                             : (GMX_THREAD_MPI ? c_gmxHostMpiThread : c_gmxHostMpiNone);
    EXPECT_EQ(info.mpi, expectedMpi);
}

//! Returns whether a plugin built with these headers accepts \p host, and why not
std::pair<bool, std::string> pluginAccepts(const GmxHostApiInfo& host)
{
    const char* reason   = nullptr;
    const bool  accepted = isHostCompatibleWithPluginBuild(host, &reason);
    return { accepted, reason != nullptr ? reason : "" };
}

TEST(PluginCompatibilityCheckTest, AcceptsThisBuild)
{
    EXPECT_TRUE(pluginAccepts(hostApiInfo()).first);
}

TEST(PluginCompatibilityCheckTest, RefusesHostOlderThanThePlugin)
{
    GmxHostApiInfo host = hostApiInfo();
    host.structSize     = offsetof(GmxHostApiInfo, mpi);
    const auto result   = pluginAccepts(host);
    EXPECT_FALSE(result.first);
    EXPECT_THAT(result.second, ::testing::HasSubstr("older"));
}

TEST(PluginCompatibilityCheckTest, RefusesAnotherMajor)
{
    GmxHostApiInfo host = hostApiInfo();
    host.apiMajor += 1;
    const auto result = pluginAccepts(host);
    EXPECT_FALSE(result.first);
    EXPECT_THAT(result.second, ::testing::HasSubstr("major version"));
}

TEST(PluginCompatibilityCheckTest, AcceptsAnotherHashOfTheSameMajor)
{
    GmxHostApiInfo host = hostApiInfo();
    host.apiHash ^= 1;
    EXPECT_TRUE(pluginAccepts(host).first);
}

TEST(PluginCompatibilityCheckTest, RefusesAnotherPrecision)
{
    GmxHostApiInfo host = hostApiInfo();
    host.sizeofReal     = sizeof(real) == sizeof(double) ? sizeof(float) : sizeof(double);
    const auto result   = pluginAccepts(host);
    EXPECT_FALSE(result.first);
    EXPECT_THAT(result.second, ::testing::HasSubstr("precision"));
}

TEST(PluginCompatibilityCheckTest, RefusesAnotherMpi)
{
    GmxHostApiInfo host = hostApiInfo();
    host.mpi          = host.mpi == c_gmxHostMpiLibrary ? c_gmxHostMpiThread : c_gmxHostMpiLibrary;
    const auto result = pluginAccepts(host);
    EXPECT_FALSE(result.first);
    EXPECT_THAT(result.second, ::testing::HasSubstr("MPI"));
}

} // namespace
} // namespace test
} // namespace gmx
