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
 * Defines the file objects of the trajectory formats that trxio.h reads
 * and writes.
 *
 * Each format reads and writes with the existing functions of xtcio.h,
 * trrio.h, tngio.h and h5md.h, so its files are the same as before.
 *
 * \author Amr Alhossary <amr_alhossary@hotmail.com>
 * \ingroup module_fileio
 */
#include "gmxpre.h"

#include "trajectoryfileformats.h"

#include <memory>
#include <optional>
#include <type_traits>
#include <utility>
#include <vector>

#include "gromacs/fileio/gmxfio.h"
#include "gromacs/fileio/gmxfio_xdr.h"
#include "gromacs/fileio/h5md/binary_information.h"
#include "gromacs/fileio/h5md/h5md.h"
#include "gromacs/fileio/tngio.h"
#include "gromacs/fileio/trajectoryfile.h"
#include "gromacs/fileio/trrio.h"
#include "gromacs/fileio/trxio.h"
#include "gromacs/fileio/xdrf.h"
#include "gromacs/fileio/xtcio.h"
#include "gromacs/mdtypes/inputrec.h"
#include "gromacs/topology/topology.h"
#include "gromacs/trajectory/trajectoryframe.h"
#include "gromacs/utility/exceptions.h"
#include "gromacs/utility/fatalerror.h"
#include "gromacs/utility/gmxassert.h"
#include "gromacs/utility/smalloc.h"
#include "gromacs/utility/vectypes.h"

namespace gmx
{

namespace
{

//! Returns the vectors of \p all at \p index, in that order.
std::vector<RVec> selected(const rvec* all, ArrayRef<const int> index)
{
    std::vector<RVec> result;
    result.reserve(index.size());
    for (const int i : index)
    {
        result.emplace_back(all[i]);
    }
    return result;
}

//! A file that xtcio.h or trrio.h opens on a t_fileio.
class XdrTrajectoryFile : public TrajectoryFile
{
public:
    XdrTrajectoryFile(std::filesystem::path path, FileMode mode) :
        path_(std::move(path)), mode_(mode)
    {
    }
    ~XdrTrajectoryFile() override
    {
        if (fio_ != nullptr)
        {
            gmx_fio_close(fio_);
        }
    }

    void read() override { fio_ = gmx_fio_open(path_, "r"); }
    void write() override { fio_ = gmx_fio_open(path_, mode_ == FileMode::Append ? "a" : "w"); }
    bool supportsMode(FileMode /*mode*/) const override { return true; }

    t_fileio* legacyFileio() const override { return fio_; }

protected:
    //! The path of the file.
    std::filesystem::path path_;
    //! The mode the file is opened in.
    FileMode mode_;
    //! The open file, or nullptr before read() or write().
    t_fileio* fio_ = nullptr;
};

//! An XTC file.
class XtcFile : public XdrTrajectoryFile
{
public:
    using XdrTrajectoryFile::XdrTrajectoryFile;

    const char* formatName() const override { return "XTC"; }

    bool readNextFrame(t_trxframe* frame, int /*flags*/) override
    {
        gmx_bool ok   = FALSE;
        real     time = 0;
        bool     haveFrame;
        if (numAtoms_ < 0)
        {
            haveFrame =
                    (read_first_xtc(
                             fio_, &frame->natoms, &frame->step, &time, frame->box, &frame->x, &frame->prec, &ok)
                     != 0);
            if (!haveFrame)
            {
                GMX_RELEASE_ASSERT(!ok,
                                   "Inconsistent results - OK status from read_first_xtc, but 0 "
                                   "atom coords read");
                frame->not_ok = DATA_NOT_OK;
            }
            numAtoms_ = frame->natoms;
        }
        else
        {
            haveFrame =
                    (read_next_xtc(fio_, numAtoms_, &frame->step, &time, frame->box, frame->x, &frame->prec, &ok)
                     != 0);
            if (!ok)
            {
                /* Actually the header could also be not ok,
                   but from bOK from read_next_xtc this can't be distinguished */
                frame->not_ok = DATA_NOT_OK;
            }
        }
        if (haveFrame)
        {
            frame->bPrec        = (frame->prec > 0);
            frame->bStep        = TRUE;
            frame->bTime        = TRUE;
            frame->timeIsDouble = false;
            frame->bX           = TRUE;
            frame->bBox         = TRUE;
            frame->time         = static_cast<double>(time);
        }
        return haveFrame;
    }

    bool seekToTime(real time) override
    {
        if (xtc_seek_time(fio_, time, numAtoms_, TRUE))
        {
            gmx_fatal(FARGS,
                      "Specified frame (time %f) doesn't exist or file "
                      "corrupt/inconsistent.",
                      time);
        }
        return true;
    }

    std::optional<double> timeOfFinalFrame() override
    {
        gmx_bool     ok = FALSE;
        const double time =
                xdr_xtc_get_last_frame_time(gmx_fio_getfp(fio_), gmx_fio_getxdr(fio_), numAtoms_, &ok);
        if (!ok)
        {
            gmx_fatal(FARGS, "Error reading last frame. Maybe seek not supported.");
        }
        return time;
    }

    void writeFrame(const t_trxframe& frame, ArrayRef<const int> index) override
    {
        if (!frame.bX)
        {
            GMX_THROW(InvalidInputError("Need coordinates to write an XTC trajectory"));
        }
        const real prec = frame.bPrec ? frame.prec : 1000.0;
        if (index.empty())
        {
            write_xtc(fio_, frame.natoms, frame.step, frame.time, frame.box, frame.x, prec);
        }
        else
        {
            const std::vector<RVec> x = selected(frame.x, index);
            write_xtc(fio_, index.ssize(), frame.step, frame.time, frame.box, as_rvec_array(x.data()), prec);
        }
    }

private:
    //! The number of atoms, known after the first frame is read.
    int numAtoms_ = -1;
};

//! A TRR file.
class TrrFile : public XdrTrajectoryFile
{
public:
    using XdrTrajectoryFile::XdrTrajectoryFile;

    const char* formatName() const override { return "TRR"; }

    bool readNextFrame(t_trxframe* frame, int flags) override
    {
        gmx_trr_header_t header;
        gmx_bool         ok = FALSE;
        if (!gmx_trr_read_frame_header(fio_, &header, &ok))
        {
            if (!ok)
            {
                frame->not_ok = HEADER_NOT_OK;
            }
            return false;
        }
        frame->natoms = header.natoms;
        frame->bStep  = TRUE;
        frame->step   = header.step;
        frame->bTime  = TRUE;
        frame->timeIsDouble = (std::is_same_v<decltype(header.t), double> && gmx_fio_is_double(fio_));
        frame->time      = static_cast<double>(header.t);
        frame->bLambda   = TRUE;
        frame->bFepState = TRUE;
        frame->lambda    = header.lambda;
        frame->bBox      = header.box_size > 0;
        if (flags & (TRX_READ_X | TRX_NEED_X))
        {
            if (frame->x == nullptr)
            {
                snew(frame->x, header.natoms);
            }
            frame->bX = header.x_size > 0;
        }
        if (flags & (TRX_READ_V | TRX_NEED_V))
        {
            if (frame->v == nullptr)
            {
                snew(frame->v, header.natoms);
            }
            frame->bV = header.v_size > 0;
        }
        if (flags & (TRX_READ_F | TRX_NEED_F))
        {
            if (frame->f == nullptr)
            {
                snew(frame->f, header.natoms);
            }
            frame->bF = header.f_size > 0;
        }
        if (!gmx_trr_read_frame_data(fio_, &header, frame->box, frame->x, frame->v, frame->f))
        {
            frame->not_ok = DATA_NOT_OK;
            return false;
        }
        return true;
    }

    void writeFrame(const t_trxframe& frame, ArrayRef<const int> index) override
    {
        if (index.empty())
        {
            gmx_trr_write_frame(fio_,
                                frame.step,
                                frame.time,
                                frame.lambda,
                                frame.box,
                                frame.natoms,
                                frame.bX ? frame.x : nullptr,
                                frame.bV ? frame.v : nullptr,
                                frame.bF ? frame.f : nullptr);
            return;
        }
        const std::vector<RVec> x = frame.bX ? selected(frame.x, index) : std::vector<RVec>();
        const std::vector<RVec> v = frame.bV ? selected(frame.v, index) : std::vector<RVec>();
        const std::vector<RVec> f = frame.bF ? selected(frame.f, index) : std::vector<RVec>();
        gmx_trr_write_frame(fio_,
                            frame.step,
                            frame.time,
                            frame.lambda,
                            frame.box,
                            index.ssize(),
                            frame.bX ? as_rvec_array(x.data()) : nullptr,
                            frame.bV ? as_rvec_array(v.data()) : nullptr,
                            frame.bF ? as_rvec_array(f.data()) : nullptr);
    }
};

//! A TNG file. Writing goes through trjtools_gmx_prepare_tng_writing().
class TngFile : public TrajectoryFile
{
public:
    TngFile(std::filesystem::path path, FileMode /*mode*/) : path_(std::move(path)) {}
    ~TngFile() override
    {
        if (tng_ != nullptr)
        {
            gmx_tng_close(&tng_);
        }
    }

    void              read() override { gmx_tng_open(path_, 'r', &tng_); }
    [[noreturn]] void write() override
    {
        GMX_THROW(NotImplementedError(
                "TNG files are written through trjtools_gmx_prepare_tng_writing()"));
    }
    bool        supportsMode(FileMode mode) const override { return mode == FileMode::Read; }
    const char* formatName() const override { return "TNG"; }

    bool readNextFrame(t_trxframe* frame, int /*flags*/) override
    {
        if (!haveReadFrame_)
        {
            // The reader looks for the frame after this step.
            frame->step = -1;
        }
        const bool haveFrame = gmx_read_next_tng_frame(tng_, frame, nullptr, 0);
        if (!haveFrame && !haveReadFrame_)
        {
            frame->not_ok = DATA_NOT_OK;
        }
        haveReadFrame_ = true;
        return haveFrame;
    }
    std::optional<double> timeOfFinalFrame() override
    {
        return gmx_tng_get_time_of_final_frame(tng_);
    }
    [[noreturn]] void writeFrame(const t_trxframe& /*frame*/, ArrayRef<const int> /*index*/) override
    {
        GMX_THROW(NotImplementedError(
                "TNG files are written through trjtools_gmx_prepare_tng_writing()"));
    }

private:
    //! The path of the file.
    std::filesystem::path path_;
    //! The open file, or nullptr before read().
    gmx_tng_trajectory_t tng_ = nullptr;
    //! Whether a frame has been read.
    bool haveReadFrame_ = false;
};

//! An H5MD file.
class H5mdFile : public TrajectoryFile
{
public:
    H5mdFile(std::filesystem::path path, FileMode /*mode*/) : path_(std::move(path)) {}

    void read() override
    {
        h5md_ = std::make_unique<H5md>(path_, H5mdFileMode::Read);
        h5md_->setupFromExistingFile();
    }
    void        write() override { h5md_ = std::make_unique<H5md>(path_, H5mdFileMode::Write); }
    bool        supportsMode(FileMode mode) const override { return mode != FileMode::Append; }
    const char* formatName() const override { return "H5MD"; }

    bool readNextFrame(t_trxframe* frame, int /*flags*/) override
    {
        const bool haveFrame = h5md_->readNextFrame(frame);
        if (!haveFrame && !haveReadFrame_)
        {
            frame->not_ok = DATA_NOT_OK;
        }
        haveReadFrame_ = true;
        return haveFrame;
    }

    void writeFrame(const t_trxframe& frame, ArrayRef<const int> index) override
    {
        const int numAtoms = index.empty() ? frame.natoms : index.ssize();
        if (!isSetUp_)
        {
            // Only the number of particles and the periodicity are known
            // here, so the file gets no topology. A frame without a PBC
            // type but with a box is taken as periodic in all directions.
            gmx_mtop_t topology;
            topology.natoms = numAtoms;
            t_inputrec inputRecord;
            inputRecord.pbcType = frame.bPBC ? frame.pbcType : (frame.bBox ? PbcType::Xyz : PbcType::No);
            h5md_->setupFileFromInput(topology, inputRecord);
            isSetUp_ = true;
        }
        const std::vector<RVec> x =
                (frame.bX && !index.empty()) ? selected(frame.x, index) : std::vector<RVec>();
        const std::vector<RVec> v =
                (frame.bV && !index.empty()) ? selected(frame.v, index) : std::vector<RVec>();
        const std::vector<RVec> f =
                (frame.bF && !index.empty()) ? selected(frame.f, index) : std::vector<RVec>();
        h5md_->writeNextFrame(vectorsToWrite(frame.bX, frame.x, x, numAtoms),
                              vectorsToWrite(frame.bV, frame.v, v, numAtoms),
                              vectorsToWrite(frame.bF, frame.f, f, numAtoms),
                              frame.box,
                              frame.step,
                              frame.time,
                              frame.bPrec && sz3CompressionIsAvailable() ? std::make_optional(frame.prec)
                                                                         : std::nullopt);
    }

private:
    //! Returns \p all, or \p chosen when particles were selected, or nothing when \p have is false.
    static ArrayRef<const RVec> vectorsToWrite(bool have, const rvec* all, const std::vector<RVec>& chosen, int numAtoms)
    {
        if (!have)
        {
            return {};
        }
        if (!chosen.empty())
        {
            return chosen;
        }
        return constArrayRefFromArray(reinterpret_cast<const RVec*>(all), numAtoms);
    }

    //! The path of the file.
    std::filesystem::path path_;
    //! The open file, or nullptr before read() or write().
    std::unique_ptr<H5md> h5md_;
    //! Whether a frame has been read.
    bool haveReadFrame_ = false;
    //! Whether the file has been set up for writing.
    bool isSetUp_ = false;
};

} // namespace

FilePointer createXtcFile(const std::filesystem::path& path, FileMode mode)
{
    return std::make_unique<XtcFile>(path, mode);
}

FilePointer createTrrFile(const std::filesystem::path& path, FileMode mode)
{
    return std::make_unique<TrrFile>(path, mode);
}

FilePointer createTngFile(const std::filesystem::path& path, FileMode mode)
{
    return std::make_unique<TngFile>(path, mode);
}

FilePointer createH5mdFile(const std::filesystem::path& path, FileMode mode)
{
    return std::make_unique<H5mdFile>(path, mode);
}

} // namespace gmx
