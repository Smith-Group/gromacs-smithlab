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
 * Tests that the trajectory formats read and written through the file
 * objects are the ones xtcio.h and trrio.h read and write.
 *
 * \author Amr Alhossary <amr_alhossary@hotmail.com>
 * \ingroup module_fileio
 */
#include "gmxpre.h"

#include "config.h"

#include <array>
#include <filesystem>
#include <fstream>
#include <iterator>
#include <memory>
#include <string>

#include <gtest/gtest.h>

#include "gromacs/fileio/fileformatregistry.h"
#include "gromacs/fileio/oenv.h"
#include "gromacs/fileio/timecontrol.h"
#include "gromacs/fileio/trajectoryfile.h"
#include "gromacs/fileio/trrio.h"
#include "gromacs/fileio/trxio.h"
#include "gromacs/fileio/xtcio.h"
#include "gromacs/trajectory/trajectoryframe.h"
#include "gromacs/utility/exceptions.h"
#include "gromacs/utility/smalloc.h"
#include "gromacs/utility/vec.h"

#include "testutils/testasserts.h"
#include "testutils/testfilemanager.h"

namespace gmx
{
namespace test
{
namespace
{

//! Returns the bytes of the file at \p path.
std::string contentsOf(const std::filesystem::path& path)
{
    std::ifstream in(path, std::ios::binary);
    return std::string(std::istreambuf_iterator<char>(in), {});
}

class TrajectoryFormatsTest : public ::testing::Test
{
public:
    static constexpr int c_numAtoms  = 6;
    static constexpr int c_numFrames = 3;

    TrajectoryFormatsTest()
    {
        clear_trxframe(&frame_, true);
        frame_.natoms = c_numAtoms;
        frame_.bStep  = true;
        frame_.bTime  = true;
        frame_.bPrec  = true;
        frame_.prec   = 1000;
        frame_.bBox   = true;
        frame_.bX     = true;
        frame_.bV     = true;
        frame_.bF     = true;
        frame_.x      = as_rvec_array(x_.data());
        frame_.v      = as_rvec_array(v_.data());
        frame_.f      = as_rvec_array(f_.data());
        clear_mat(frame_.box);
    }

    //! Makes frame_ the frame numbered \p number.
    void setFrame(int number)
    {
        frame_.step = 10 * number;
        frame_.time = 0.5 * number;
        for (int i = 0; i < c_numAtoms; i++)
        {
            x_[i] = { (1.0F * i) + number, 0.5F * i, 0.25F * i };
            v_[i] = { -1.0F * i, 0.1F * number, 0.0F };
            f_[i] = { 100.0F, 10.0F * i, 1.0F * number };
        }
        frame_.box[XX][XX] = frame_.box[YY][YY] = frame_.box[ZZ][ZZ] = 3.0F + number;
    }

    //! Writes the frames through trxio.h to \p path.
    void writeThroughTrxio(const std::filesystem::path& path)
    {
        t_trxstatus* status = open_trx(path, "w");
        for (int number = 0; number < c_numFrames; number++)
        {
            setFrame(number);
            write_trxframe(status, &frame_, nullptr);
        }
        close_trx(status);
    }

    TestFileManager              fileManager_;
    t_trxframe                   frame_;
    std::array<RVec, c_numAtoms> x_;
    std::array<RVec, c_numAtoms> v_;
    std::array<RVec, c_numAtoms> f_;
};

TEST_F(TrajectoryFormatsTest, WritesXtcAsTheXtcWriter)
{
    const std::filesystem::path newFile = fileManager_.getTemporaryFilePath("new.xtc");
    const std::filesystem::path oldFile = fileManager_.getTemporaryFilePath("old.xtc");
    writeThroughTrxio(newFile);

    t_fileio* out = open_xtc(oldFile, "w");
    for (int number = 0; number < c_numFrames; number++)
    {
        setFrame(number);
        write_xtc(out, c_numAtoms, frame_.step, frame_.time, frame_.box, frame_.x, frame_.prec);
    }
    close_xtc(out);

    EXPECT_EQ(contentsOf(newFile), contentsOf(oldFile));
}

TEST_F(TrajectoryFormatsTest, WritesTrrAsTheTrrWriter)
{
    const std::filesystem::path newFile = fileManager_.getTemporaryFilePath("new.trr");
    const std::filesystem::path oldFile = fileManager_.getTemporaryFilePath("old.trr");
    writeThroughTrxio(newFile);

    t_fileio* out = gmx_trr_open(oldFile, "w");
    for (int number = 0; number < c_numFrames; number++)
    {
        setFrame(number);
        gmx_trr_write_frame(
                out, frame_.step, frame_.time, frame_.lambda, frame_.box, c_numAtoms, frame_.x, frame_.v, frame_.f);
    }
    gmx_trr_close(out);

    EXPECT_EQ(contentsOf(newFile), contentsOf(oldFile));
}

TEST_F(TrajectoryFormatsTest, WritesSelectedParticlesInOrder)
{
    const std::filesystem::path newFile  = fileManager_.getTemporaryFilePath("new.trr");
    const std::filesystem::path oldFile  = fileManager_.getTemporaryFilePath("old.trr");
    const int                   index[2] = { 5, 1 };
    setFrame(1);
    frame_.bF = false;

    t_trxstatus* status = open_trx(newFile, "w");
    write_trxframe_indexed(status, &frame_, 2, index, nullptr);
    close_trx(status);

    const rvec x[2] = { { x_[5][XX], x_[5][YY], x_[5][ZZ] }, { x_[1][XX], x_[1][YY], x_[1][ZZ] } };
    const rvec v[2] = { { v_[5][XX], v_[5][YY], v_[5][ZZ] }, { v_[1][XX], v_[1][YY], v_[1][ZZ] } };
    t_fileio*  out  = gmx_trr_open(oldFile, "w");
    gmx_trr_write_frame(out, frame_.step, frame_.time, frame_.lambda, frame_.box, 2, x, v, nullptr);
    gmx_trr_close(out);

    EXPECT_EQ(contentsOf(newFile), contentsOf(oldFile));
}

TEST_F(TrajectoryFormatsTest, TrrReadsOnlyWhatIsAskedFor)
{
    const std::filesystem::path path = fileManager_.getTemporaryFilePath("frames.trr");
    writeThroughTrxio(path);

    std::unique_ptr<TrajectoryFile> file = openFileAs<TrajectoryFile>(path, FileMode::Read);
    file->read();
    t_trxframe frame;
    clear_trxframe(&frame, true);
    ASSERT_TRUE(file->readNextFrame(&frame, TRX_READ_X));
    EXPECT_EQ(frame.natoms, c_numAtoms);
    EXPECT_EQ(frame.step, 0);
    EXPECT_TRUE(frame.bX);
    EXPECT_FALSE(frame.bV);
    EXPECT_FALSE(frame.bF);
    EXPECT_REAL_EQ(frame.x[c_numAtoms - 1][XX], 1.0F * (c_numAtoms - 1));

    clear_trxframe(&frame, false);
    ASSERT_TRUE(file->readNextFrame(&frame, TRX_READ_X | TRX_READ_V));
    EXPECT_EQ(frame.step, 10);
    EXPECT_TRUE(frame.bV);
    EXPECT_FALSE(frame.bF);
    EXPECT_REAL_EQ(frame.v[0][YY], 0.1F);

    clear_trxframe(&frame, false);
    ASSERT_TRUE(file->readNextFrame(&frame, TRX_READ_X));
    clear_trxframe(&frame, false);
    EXPECT_FALSE(file->readNextFrame(&frame, TRX_READ_X));
    sfree(frame.x);
    sfree(frame.v);
}

TEST_F(TrajectoryFormatsTest, XtcSeeksToTheStartTime)
{
    const std::filesystem::path path = fileManager_.getTemporaryFilePath("frames.xtc");
    writeThroughTrxio(path);

    gmx_output_env_t* oenv;
    output_env_init_default(&oenv);
    TimeControl timeControl;
    timeControl.begin = 1.0;
    t_trxstatus* status;
    t_trxframe   frame;
    ASSERT_TRUE(read_first_frame(oenv, &status, path, &frame, &timeControl, TRX_NEED_X));
    EXPECT_EQ(frame.step, 20);
    EXPECT_REAL_EQ(frame.x[0][XX], 2.0F);
    EXPECT_FALSE(read_next_frame(oenv, status, &frame));
    close_trx(status);
    sfree(frame.x);
    output_env_done(oenv);
}

TEST_F(TrajectoryFormatsTest, XtcRefusesAFrameWithoutPositions)
{
    const std::filesystem::path     path = fileManager_.getTemporaryFilePath("none.xtc");
    std::unique_ptr<TrajectoryFile> file = openFileAs<TrajectoryFile>(path, FileMode::Write);
    file->write();
    setFrame(0);
    frame_.bX = false;
    EXPECT_THROW_GMX(file->writeFrame(frame_, {}), InvalidInputError);
}

TEST_F(TrajectoryFormatsTest, XtcKnowsTheTimeOfItsFinalFrame)
{
    const std::filesystem::path path = fileManager_.getTemporaryFilePath("frames.xtc");
    writeThroughTrxio(path);

    gmx_output_env_t* oenv;
    output_env_init_default(&oenv);
    t_trxstatus* status;
    t_trxframe   frame;
    ASSERT_TRUE(read_first_frame(oenv, &status, path, &frame, nullptr, TRX_NEED_X));
    EXPECT_REAL_EQ(trx_get_time_of_final_frame(status), 1.0);
    close_trx(status);
    sfree(frame.x);
    output_env_done(oenv);
}

TEST_F(TrajectoryFormatsTest, TngReadsThroughTheFileObject)
{
    if (!GMX_USE_TNG)
    {
        GTEST_SKIP() << "Built without TNG";
    }
    const std::filesystem::path path =
            TestFileManager::getTestSimulationDatabaseDirectory().append("spc2-traj.tng");
    std::unique_ptr<TrajectoryFile> file = openFileAs<TrajectoryFile>(path, FileMode::Read);
    file->read();
    t_trxframe frame;
    clear_trxframe(&frame, true);
    ASSERT_TRUE(file->readNextFrame(&frame, TRX_NEED_X));
    EXPECT_EQ(frame.natoms, 6);
    EXPECT_TRUE(frame.bX);
    int numFrames = 1;
    while (file->readNextFrame(&frame, TRX_NEED_X))
    {
        numFrames++;
    }
    EXPECT_GT(numFrames, 1);
    done_frame(&frame);
}

TEST_F(TrajectoryFormatsTest, TngGivesTheTngWriterItsOpenFile)
{
    if (!GMX_USE_TNG)
    {
        GTEST_SKIP() << "Built without TNG";
    }
    const std::filesystem::path path =
            TestFileManager::getTestSimulationDatabaseDirectory().append("spc2-traj.tng");
    std::unique_ptr<TrajectoryFile> file = openFileAs<TrajectoryFile>(path, FileMode::Read);
    EXPECT_EQ(file->legacyTng(), nullptr);
    file->read();
    EXPECT_NE(file->legacyTng(), nullptr);
    EXPECT_EQ(file->legacyH5md(), nullptr);
}

TEST_F(TrajectoryFormatsTest, H5mdGivesTheH5mdWriterItsOpenFile)
{
    if (!GMX_USE_HDF5)
    {
        GTEST_SKIP() << "Built without HDF5";
    }
    const std::filesystem::path path =
            TestFileManager::getTestSimulationDatabaseDirectory().append("spc2-traj.h5md");
    std::unique_ptr<TrajectoryFile> file = openFileAs<TrajectoryFile>(path, FileMode::Read);
    EXPECT_EQ(file->legacyH5md(), nullptr);
    file->read();
    EXPECT_NE(file->legacyH5md(), nullptr);
    EXPECT_EQ(file->legacyTng(), nullptr);
}

TEST_F(TrajectoryFormatsTest, H5mdWritesAndReadsFrames)
{
    if (!GMX_USE_HDF5)
    {
        GTEST_SKIP() << "Built without HDF5";
    }
    const std::filesystem::path path = fileManager_.getTemporaryFilePath("frames.h5md");
    writeThroughTrxio(path);

    gmx_output_env_t* oenv;
    output_env_init_default(&oenv);
    t_trxstatus* status;
    t_trxframe   frame;
    ASSERT_TRUE(read_first_frame(oenv, &status, path, &frame, nullptr, TRX_NEED_X | TRX_READ_V));
    EXPECT_EQ(frame.natoms, c_numAtoms);
    EXPECT_EQ(frame.step, 0);
    EXPECT_TRUE(frame.bV);
    // With SZ3, positions are compressed to the precision of the frame
    EXPECT_NEAR(frame.x[1][XX], 1.0, 1.0 / frame_.prec);
    ASSERT_TRUE(read_next_frame(oenv, status, &frame));
    EXPECT_EQ(frame.step, 10);
    EXPECT_NEAR(frame.x[1][XX], 2.0, 1.0 / frame_.prec);
    EXPECT_REAL_EQ(frame.box[XX][XX], 4.0F);
    ASSERT_TRUE(read_next_frame(oenv, status, &frame));
    EXPECT_FALSE(read_next_frame(oenv, status, &frame));
    close_trx(status);
    done_frame(&frame);
    output_env_done(oenv);
}

} // namespace
} // namespace test
} // namespace gmx
