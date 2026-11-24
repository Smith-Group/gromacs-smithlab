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
 * Declares gmx::TrajectoryFile.
 *
 * \author Amr Alhossary <amr_alhossary@hotmail.com>
 * \inlibraryapi
 * \ingroup module_fileio
 */
#ifndef GMX_FILEIO_TRAJECTORYFILE_H
#define GMX_FILEIO_TRAJECTORYFILE_H

#include <optional>

#include "gromacs/fileio/file.h"
#include "gromacs/utility/arrayref.h"
#include "gromacs/utility/real.h"

struct gmx_conect_t;
struct gmx_tng_trajectory;
struct t_fileio;
struct t_trxframe;

namespace gmx
{

class H5md;

/*! \libinternal \brief
 * A file that carries a sequence of frames.
 *
 * openForReadingFrames() and openForWritingFrames() open the file; the
 * frames are then read one at a time with readNextFrame() or written with
 * writeFrame(). The frame is
 * a t_trxframe until gmx::TrajectoryFrame owns its data.
 */
class TrajectoryFile : public virtual File
{
public:
    //! The capability flag of this interface.
    static constexpr FileCapability sc_capability = FileCapability::Trajectory;

    ~TrajectoryFile() override;

    /*! \brief
     * Opens the file to read frames.
     *
     * This is read(), except for a format that also holds a structure:
     * there read() reads the structure.
     */
    virtual void openForReadingFrames() { read(); }
    /*! \brief
     * Opens the file to write frames.
     *
     * This is write(), except for a format that also holds a structure:
     * there write() writes the structure.
     */
    virtual void openForWritingFrames() { write(); }

    /*! \brief
     * Reads the next frame into \p frame.
     *
     * \param[in,out] frame  The frame to fill. Arrays it already has are reused.
     * \param[in]     flags  The TRX_READ_* and TRX_NEED_* flags of trxio.h; a
     *                       format that stores more reads only what they ask for.
     * \returns false at the end of the file, or when the frame is not
     *          complete, which \c frame->not_ok then says.
     */
    virtual bool readNextFrame(t_trxframe* frame, int flags) = 0;

    /*! \brief
     * Moves to the first frame at or after \p time, when the format can
     * do that without reading the frames before it.
     *
     * \returns whether the position changed.
     */
    virtual bool seekToTime(real /*time*/) { return false; }

    /*! \brief
     * Returns the time of the last frame, when the format can find it
     * without reading the frames before it.
     */
    virtual std::optional<double> timeOfFinalFrame() { return std::nullopt; }

    /*! \brief
     * Writes \p frame.
     *
     * \param[in] frame        The frame to write.
     * \param[in] index        The particles to write, in that order; all when empty.
     * \param[in] connections  The bonds that PDB and mmCIF write as CONECT
     *                         records, or nullptr; other formats ignore it.
     * \throws InvalidInputError if \p frame lacks what the format needs.
     */
    virtual void writeFrame(const t_trxframe& frame, ArrayRef<const int> index, gmx_conect_t* connections) = 0;

    /*! \brief
     * Returns the t_fileio of the open file, or nullptr if the format
     * does not use one.
     *
     * This is the seam with trxio.h: trx_get_fileio() and rewind_trj()
     * work on it.
     */
    virtual t_fileio* legacyFileio() const { return nullptr; }
    /*! \brief
     * Returns the open TNG file, or nullptr for another format.
     *
     * trjtools_gmx_prepare_tng_writing() copies the setup of an input
     * TNG file from it.
     */
    virtual gmx_tng_trajectory* legacyTng() const { return nullptr; }
    /*! \brief
     * Returns the open H5MD file, or nullptr for another format.
     *
     * trjtools_gmx_prepare_h5md_writing() copies the metadata of an
     * input H5MD file from it.
     */
    virtual const H5md* legacyH5md() const { return nullptr; }
};

} // namespace gmx

#endif
