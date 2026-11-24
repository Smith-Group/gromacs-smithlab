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
 * Tests for the file objects of the structure formats that confio.h reads
 * and writes.
 *
 * \author Amr Alhossary <amr_alhossary@hotmail.com>
 * \ingroup module_fileio
 */
#include "gmxpre.h"

#include <cstdio>

#include <filesystem>
#include <memory>
#include <string>
#include <vector>

#include <gtest/gtest.h>

#include "gromacs/fileio/boxcontainingfile.h"
#include "gromacs/fileio/cifio.h"
#include "gromacs/fileio/confio.h"
#include "gromacs/fileio/espio.h"
#include "gromacs/fileio/fileformatregistry.h"
#include "gromacs/fileio/filetypes.h"
#include "gromacs/fileio/g96io.h"
#include "gromacs/fileio/groio.h"
#include "gromacs/fileio/pdbio.h"
#include "gromacs/fileio/structurefile.h"
#include "gromacs/fileio/trxio.h"
#include "gromacs/mdtypes/md_enums.h"
#include "gromacs/topology/atoms.h"
#include "gromacs/topology/symtab.h"
#include "gromacs/trajectory/trajectoryframe.h"
#include "gromacs/utility/exceptions.h"
#include "gromacs/utility/smalloc.h"
#include "gromacs/utility/stringutil.h"
#include "gromacs/utility/textreader.h"
#include "gromacs/utility/vec.h"

#include "testutils/testasserts.h"
#include "testutils/testfilemanager.h"

namespace gmx
{
namespace test
{
namespace
{

//! A structure read from the test simulation database with readConfAndAtoms().
class StructureFormatsTest : public ::testing::Test
{
public:
    StructureFormatsTest()
    {
        open_symtab(&symtab_);
        readConfAndAtoms(
                TestFileManager::getTestSimulationDatabaseDirectory().append("lysozyme.pdb"),
                &symtab_,
                &title_,
                &atoms_,
                &pbcType_,
                &x_,
                nullptr,
                box_);
        // Any velocities will do to check that they are written.
        snew(v_, atoms_.nr);
        for (int i = 0; i < atoms_.nr; i++)
        {
            copy_rvec(x_[i], v_[i]);
        }
    }
    ~StructureFormatsTest() override
    {
        sfree(v_);
        sfree(x_);
        sfree(title_);
        done_atom(&atoms_);
        done_symtab(&symtab_);
    }

    //! Returns the contents of a file written by \p writer to \p fileName.
    template<typename Writer>
    std::string writtenBy(const std::string& fileName, Writer writer)
    {
        const std::filesystem::path path = fileManager_.getTemporaryFilePath(fileName);
        writer(path);
        return TextReader::readFileToString(path);
    }
    //! Writes with the PDB writer of pdbio.h, as PQR if \p usePqrFormat.
    void writeWithPdbWriter(const std::filesystem::path& path, bool usePqrFormat)
    {
        std::vector<int> all(atoms_.nr);
        for (int i = 0; i < atoms_.nr; i++)
        {
            all[i] = i;
        }
        FILE* out = std::fopen(path.string().c_str(), "w");
        write_pdbfile_indexed(
                out, title_, &atoms_, x_, pbcType_, box_, ' ', -1, atoms_.nr, all.data(), nullptr, usePqrFormat);
        std::fclose(out);
    }

    TestFileManager fileManager_;
    t_symtab        symtab_;
    char*           title_ = nullptr;
    t_atoms         atoms_;
    PbcType         pbcType_ = PbcType::Unset;
    rvec*           x_       = nullptr;
    rvec*           v_       = nullptr;
    matrix          box_;
};

TEST_F(StructureFormatsTest, ReadsAsThePdbReader)
{
    t_symtab symtab;
    open_symtab(&symtab);
    t_atoms                     atoms;
    int                         numAtoms = 0;
    const std::filesystem::path path =
            TestFileManager::getTestSimulationDatabaseDirectory().append("lysozyme.pdb");
    FILE* in = std::fopen(path.string().c_str(), "r");
    get_pdb_coordnum(in, &numAtoms);
    std::fclose(in);
    init_t_atoms(&atoms, numAtoms, true);
    std::vector<RVec> x(numAtoms);
    char*             title = nullptr;
    matrix            box;
    PbcType           pbcType = PbcType::Unset;
    gmx_pdb_read_conf(path, &symtab, &title, &atoms, as_rvec_array(x.data()), &pbcType, box);

    EXPECT_STREQ(title, title_);
    EXPECT_EQ(pbcType, pbcType_);
    ASSERT_EQ(atoms.nr, atoms_.nr);
    ASSERT_EQ(atoms.nres, atoms_.nres);
    for (int i = 0; i < atoms.nr; i++)
    {
        EXPECT_STREQ(*atoms.atomname[i], *atoms_.atomname[i]);
        EXPECT_EQ(atoms.atom[i].resind, atoms_.atom[i].resind);
        EXPECT_EQ(atoms.pdbinfo[i].occup, atoms_.pdbinfo[i].occup);
        EXPECT_REAL_EQ(x[i][XX], x_[i][XX]);
    }
    for (int r = 0; r < atoms.nres; r++)
    {
        EXPECT_STREQ(*atoms.resinfo[r].name, *atoms_.resinfo[r].name);
        EXPECT_EQ(atoms.resinfo[r].chainid, atoms_.resinfo[r].chainid);
    }
    for (int d = 0; d < DIM; d++)
    {
        EXPECT_REAL_EQ(box[d][d], box_[d][d]);
    }
    sfree(title);
    done_atom(&atoms);
    done_symtab(&symtab);
}

TEST_F(StructureFormatsTest, WritesGroAsTheGroWriter)
{
    EXPECT_EQ(writtenBy("new.gro",
                        [&](const std::filesystem::path& path)
                        { write_sto_conf(path, title_, &atoms_, x_, v_, pbcType_, box_); }),
              writtenBy("old.gro",
                        [&](const std::filesystem::path& path)
                        { write_conf_p(path, title_, &atoms_, x_, v_, box_); }));
}

TEST_F(StructureFormatsTest, WritesG96AsTheG96Writer)
{
    EXPECT_EQ(writtenBy("new.g96",
                        [&](const std::filesystem::path& path)
                        { write_sto_conf(path, title_, &atoms_, x_, v_, pbcType_, box_); }),
              writtenBy("old.g96",
                        [&](const std::filesystem::path& path)
                        {
                            t_trxframe frame;
                            clear_trxframe(&frame, TRUE);
                            frame.natoms = atoms_.nr;
                            frame.bAtoms = TRUE;
                            frame.atoms  = &atoms_;
                            frame.bX     = TRUE;
                            frame.x      = x_;
                            frame.bV     = TRUE;
                            frame.v      = v_;
                            frame.bBox   = TRUE;
                            copy_mat(box_, frame.box);
                            FILE* out = std::fopen(path.string().c_str(), "w");
                            write_g96_conf(out, title_, &frame, -1, nullptr);
                            std::fclose(out);
                        }));
}

TEST_F(StructureFormatsTest, WritesPdbAsThePdbWriter)
{
    EXPECT_EQ(writtenBy("new.pdb",
                        [&](const std::filesystem::path& path)
                        { write_sto_conf(path, title_, &atoms_, x_, nullptr, pbcType_, box_); }),
              writtenBy("old.pdb",
                        [&](const std::filesystem::path& path) { writeWithPdbWriter(path, false); }));
}

TEST_F(StructureFormatsTest, WritesCifAsTheCifWriter)
{
    EXPECT_EQ(writtenBy("new.cif",
                        [&](const std::filesystem::path& path)
                        { write_sto_conf(path, title_, &atoms_, x_, nullptr, pbcType_, box_); }),
              writtenBy("old.cif",
                        [&](const std::filesystem::path& path)
                        { CifIO::write(path, title_, &atoms_, x_, pbcType_, box_); }));
}

TEST_F(StructureFormatsTest, ReadsCifAsTheCifReader)
{
    const std::filesystem::path path = fileManager_.getTemporaryFilePath("conf.cif");
    CifIO::write(path, title_, &atoms_, x_, pbcType_, box_);

    t_symtab symtab;
    open_symtab(&symtab);
    t_atoms atoms;
    int     numAtoms = 0;
    CifIO::countAtoms(path, &numAtoms);
    init_t_atoms(&atoms, numAtoms, true);
    std::vector<RVec> x(numAtoms);
    char*             title = nullptr;
    matrix            box;
    PbcType           pbcType = PbcType::Unset;
    CifIO::read(path, &symtab, &title, &atoms, as_rvec_array(x.data()), &pbcType, box);

    t_symtab newSymtab;
    open_symtab(&newSymtab);
    t_atoms newAtoms;
    char*   newTitle   = nullptr;
    PbcType newPbcType = PbcType::Unset;
    rvec*   newX       = nullptr;
    matrix  newBox;
    readConfAndAtoms(path, &newSymtab, &newTitle, &newAtoms, &newPbcType, &newX, nullptr, newBox);

    EXPECT_STREQ(newTitle, title);
    EXPECT_EQ(newPbcType, pbcType);
    ASSERT_EQ(newAtoms.nr, atoms.nr);
    ASSERT_EQ(newAtoms.nres, atoms.nres);
    for (int i = 0; i < atoms.nr; i++)
    {
        EXPECT_STREQ(*newAtoms.atomname[i], *atoms.atomname[i]);
        EXPECT_EQ(newAtoms.atom[i].resind, atoms.atom[i].resind);
        EXPECT_EQ(newAtoms.pdbinfo[i].occup, atoms.pdbinfo[i].occup);
        EXPECT_REAL_EQ(newX[i][XX], x[i][XX]);
    }
    for (int r = 0; r < atoms.nres; r++)
    {
        EXPECT_STREQ(*newAtoms.resinfo[r].name, *atoms.resinfo[r].name);
        EXPECT_EQ(newAtoms.resinfo[r].chainid, atoms.resinfo[r].chainid);
    }
    for (int d = 0; d < DIM; d++)
    {
        EXPECT_REAL_EQ(newBox[d][d], box[d][d]);
    }
    sfree(newX);
    sfree(newTitle);
    done_atom(&newAtoms);
    done_symtab(&newSymtab);
    sfree(title);
    done_atom(&atoms);
    done_symtab(&symtab);
}

TEST_F(StructureFormatsTest, WritesEspAsTheEspWriter)
{
    EXPECT_EQ(writtenBy("new.esp",
                        [&](const std::filesystem::path& path)
                        { write_sto_conf(path, title_, &atoms_, x_, v_, pbcType_, box_); }),
              writtenBy("old.esp",
                        [&](const std::filesystem::path& path)
                        {
                            FILE* out = std::fopen(path.string().c_str(), "w");
                            write_espresso_conf_indexed(
                                    out, title_, &atoms_, atoms_.nr, nullptr, x_, v_, box_);
                            std::fclose(out);
                        }));
}

TEST_F(StructureFormatsTest, WritesPqrWithoutAnIndex)
{
    // write_sto_conf() used to refuse PQR; only write_sto_conf_indexed() wrote it.
    EXPECT_EQ(writtenBy("new.pqr",
                        [&](const std::filesystem::path& path)
                        { write_sto_conf(path, title_, &atoms_, x_, nullptr, pbcType_, box_); }),
              writtenBy("old.pqr",
                        [&](const std::filesystem::path& path) { writeWithPdbWriter(path, true); }));
}

TEST_F(StructureFormatsTest, WritesSelectedParticlesInOrder)
{
    int index[2] = { 5, 1 };
    EXPECT_EQ(writtenBy("new.gro",
                        [&](const std::filesystem::path& path) {
                            write_sto_conf_indexed(
                                    path, title_, &atoms_, x_, nullptr, pbcType_, box_, 2, index);
                        }),
              writtenBy("old.gro",
                        [&](const std::filesystem::path& path)
                        {
                            FILE* out = std::fopen(path.string().c_str(), "w");
                            write_hconf_indexed_p(out, title_, &atoms_, 2, index, x_, nullptr, box_);
                            std::fclose(out);
                        }));
}

TEST_F(StructureFormatsTest, WritesNoParticlesForAnEmptySelection)
{
    const std::filesystem::path path = fileManager_.getTemporaryFilePath("empty.gro");
    write_sto_conf_indexed(path, title_, &atoms_, x_, nullptr, pbcType_, box_, 0, nullptr);
    std::unique_ptr<StructureFile> file = openFileAs<StructureFile>(path, FileMode::Read);
    file->read();
    EXPECT_EQ(file->getNumAtoms(), 0);
}

TEST_F(StructureFormatsTest, GroReportsWhetherItHasVelocities)
{
    const std::filesystem::path withVelocities    = fileManager_.getTemporaryFilePath("v.gro");
    const std::filesystem::path withoutVelocities = fileManager_.getTemporaryFilePath("x.gro");
    write_sto_conf(withVelocities, title_, &atoms_, x_, v_, pbcType_, box_);
    write_sto_conf(withoutVelocities, title_, &atoms_, x_, nullptr, pbcType_, box_);

    std::unique_ptr<StructureFile> file = openFileAs<StructureFile>(withVelocities, FileMode::Read);
    file->read();
    EXPECT_EQ(file->getVelocities().ssize(), atoms_.nr);
    file = openFileAs<StructureFile>(withoutVelocities, FileMode::Read);
    file->read();
    EXPECT_TRUE(file->getVelocities().empty());
}

TEST_F(StructureFormatsTest, PdbWritesEachModel)
{
    const std::filesystem::path path = fileManager_.getTemporaryFilePath("models.pdb");
    std::vector<RVec>           shifted(x_, x_ + atoms_.nr);
    for (RVec& position : shifted)
    {
        position[XX] += 1;
    }
    {
        std::unique_ptr<MultimodelStructureFile> file =
                openFileAs<MultimodelStructureFile>(path, FileMode::Write);
        file->setStructure(
                title_, atoms_, constArrayRefFromArray(reinterpret_cast<const RVec*>(x_), atoms_.nr), {});
        EXPECT_EQ(file->addModel(shifted), 1);
        EXPECT_EQ(file->getNumModels(), 2);
        file->setActiveModel(1);
        EXPECT_REAL_EQ(shifted[0][XX], file->getPositions()[0][XX]);
        EXPECT_THROW_GMX(file->setActiveModel(2), RangeError);
        file->write();
    }
    const std::string contents = TextReader::readFileToString(path);
    EXPECT_NE(contents.find("MODEL        1"), std::string::npos);
    EXPECT_NE(contents.find("MODEL        2"), std::string::npos);

    std::unique_ptr<StructureFile> file = openFileAs<StructureFile>(path, FileMode::Read);
    file->read();
    EXPECT_REAL_EQ(x_[0][XX], file->getPositions()[0][XX]);
}

TEST_F(StructureFormatsTest, CifWritesEachModel)
{
    const std::filesystem::path path = fileManager_.getTemporaryFilePath("models.cif");
    std::vector<RVec>           shifted(x_, x_ + atoms_.nr);
    for (RVec& position : shifted)
    {
        position[XX] += 1;
    }
    {
        std::unique_ptr<MultimodelStructureFile> file =
                openFileAs<MultimodelStructureFile>(path, FileMode::Write);
        file->setStructure(
                title_, atoms_, constArrayRefFromArray(reinterpret_cast<const RVec*>(x_), atoms_.nr), {});
        file->addModel(shifted);
        file->write();
    }
    const std::unique_ptr<CifTrajectory> trajectory = CifIO::openTrajectory(path);
    ASSERT_NE(trajectory, nullptr);
    ASSERT_EQ(trajectory->numAtoms(), atoms_.nr);
    t_trxframe frame;
    clear_trxframe(&frame, TRUE);
    std::vector<RVec> x(atoms_.nr);
    frame.x = as_rvec_array(x.data());
    ASSERT_TRUE(trajectory->readNextFrame(&frame));
    EXPECT_NEAR(x[0][XX], x_[0][XX], 1e-3);
    ASSERT_TRUE(trajectory->readNextFrame(&frame));
    EXPECT_NEAR(x[0][XX], shifted[0][XX], 1e-3);
    EXPECT_FALSE(trajectory->readNextFrame(&frame));

    std::unique_ptr<StructureFile> file = openFileAs<StructureFile>(path, FileMode::Read);
    file->read();
    EXPECT_STREQ(file->getTitle().c_str(), title_);
}

TEST_F(StructureFormatsTest, PdbHasABoxOnlyWithCryst1)
{
    const std::filesystem::path path = fileManager_.getTemporaryFilePath("nobox.pdb");
    matrix                      noBox;
    clear_mat(noBox);
    write_sto_conf(path, title_, &atoms_, x_, nullptr, PbcType::Unset, noBox);

    std::unique_ptr<BoxContainingFile> file = openFileAs<BoxContainingFile>(path, FileMode::Read);
    file->read();
    EXPECT_FALSE(file->getBoundaries().has_value());
    EXPECT_EQ(file->getPbcType(), PbcType::No);
}

TEST(StructureFormatsWithoutStructureTest, RefusesToWriteBeforeAStructureIsSet)
{
    TestFileManager                fileManager;
    std::unique_ptr<StructureFile> file =
            openFileAs<StructureFile>(fileManager.getTemporaryFilePath("unset.gro"), FileMode::Write);
    EXPECT_THROW_GMX(file->write(), APIError);
}

} // namespace
} // namespace test
} // namespace gmx
