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
 * Tests for gmx::TopFile.
 *
 * \author Amr Alhossary <amr_alhossary@hotmail.com>
 * \ingroup module_fileio
 */
#include "gmxpre.h"

#include "gromacs/fileio/topfile.h"

#include <filesystem>
#include <string>
#include <utility>
#include <vector>

#include <gmock/gmock.h>
#include <gtest/gtest.h>

#include "gromacs/utility/exceptions.h"
#include "gromacs/utility/textreader.h"
#include "gromacs/utility/textwriter.h"

#include "testutils/testasserts.h"
#include "testutils/testfilemanager.h"

namespace gmx
{
namespace test
{
namespace
{

//! A topology with two lines of solvent.
const char* const c_topology =
        "#include \"oplsaa.ff/forcefield.itp\"\n"
        "\n"
        "[ system ]\n"
        "; Name\n"
        "Protein\n"
        "\n"
        "[ molecules ]\n"
        "; Compound  #mols\n"
        "Protein_A  1\n"
        "SOL        100\n"
        "NA         2\n"
        "SOL        50\n";

class TopFileTest : public ::testing::Test
{
public:
    TopFileTest() : path_(fileManager_.getTemporaryFilePath("topol.top"))
    {
        TextWriter::writeFileFromString(path_, c_topology);
        file_.read();
    }

    TestFileManager       fileManager_;
    std::filesystem::path path_;
    TopFile               file_{ path_ };
};

TEST_F(TopFileTest, ReadsMoleculesAndSystemName)
{
    EXPECT_EQ(file_.getCompleteness(), TopologyCompleteness::Full);
    EXPECT_THAT(file_.getMoleculeTypes(), ::testing::ElementsAre("Protein_A", "SOL", "NA"));
    EXPECT_EQ(file_.getNumMolecules("SOL"), 150);
    EXPECT_EQ(file_.getNumMolecules("sol"), 150);
    EXPECT_EQ(file_.getNumMolecules("CL"), 0);
    EXPECT_EQ(file_.getSystemName(), "Protein");
}

TEST_F(TopFileTest, AddsMoleculesAndRenamesTheSystem)
{
    file_.setSystemName("Protein in water");
    file_.addMolecules("SOL", 3);
    file_.write();
    EXPECT_EQ(TextReader::readFileToString(path_),
              "#include \"oplsaa.ff/forcefield.itp\"\n"
              "\n"
              "[ system ]\n"
              "; Name\n"
              "Protein in water\n"
              "\n"
              "[ molecules ]\n"
              "; Compound  #mols\n"
              "Protein_A  1\n"
              "SOL        100\n"
              "NA         2\n"
              "SOL        50\n"
              "SOL                 3\n");
}

TEST_F(TopFileTest, ReplacesMoleculesOfTheLastLine)
{
    const std::pair<std::string, int> ions[] = { { "NA", 3 }, { "CL", 0 } };
    file_.replaceMolecules("SOL", 3, ions);
    file_.write();
    EXPECT_EQ(TextReader::readFileToString(path_),
              "#include \"oplsaa.ff/forcefield.itp\"\n"
              "\n"
              "[ system ]\n"
              "; Name\n"
              "Protein\n"
              "\n"
              "[ molecules ]\n"
              "; Compound  #mols\n"
              "Protein_A  1\n"
              "SOL        100\n"
              "NA         2\n"
              "SOL                47\n"
              "NA                  3\n");
}

TEST_F(TopFileTest, RefusesToReplaceMissingMolecules)
{
    const std::pair<std::string, int> ions[] = { { "NA", 60 } };
    EXPECT_THROW_GMX(file_.replaceMolecules("SOL", 60, ions), InvalidInputError);
    EXPECT_THROW_GMX(file_.replaceMolecules("MG", 1, ions), InvalidInputError);
}

} // namespace
} // namespace test
} // namespace gmx
