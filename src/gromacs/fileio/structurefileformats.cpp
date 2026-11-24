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
 * Defines the file objects of the structure formats that confio.h reads
 * and writes. GRO, G96, PDB and mmCIF files also hold the frames that
 * trxio.h reads and writes.
 *
 * Each format reads and writes with the existing functions of groio.h,
 * g96io.h, pdbio.h, cifio.h and espio.h, so its files are the same as before.
 *
 * \author Amr Alhossary <amr_alhossary@hotmail.com>
 * \ingroup module_fileio
 */
#include "gmxpre.h"

#include "structurefileformats.h"

#include <cinttypes>
#include <cstdio>
#include <cstring>

#include <deque>
#include <memory>
#include <numeric>
#include <optional>
#include <string>
#include <utility>
#include <vector>

#include "gromacs/fileio/boxcontainingfile.h"
#include "gromacs/fileio/cifio.h"
#include "gromacs/fileio/espio.h"
#include "gromacs/fileio/filetypes.h"
#include "gromacs/fileio/g96io.h"
#include "gromacs/fileio/gmxfio.h"
#include "gromacs/fileio/groio.h"
#include "gromacs/fileio/iohelper.h"
#include "gromacs/fileio/pdbio.h"
#include "gromacs/fileio/structurefile.h"
#include "gromacs/fileio/topologyfile.h"
#include "gromacs/fileio/trajectoryfile.h"
#include "gromacs/fileio/trxio.h"
#include "gromacs/mdtypes/md_enums.h"
#include "gromacs/topology/atoms.h"
#include "gromacs/topology/symtab.h"
#include "gromacs/trajectory/trajectoryframe.h"
#include "gromacs/utility/classhelpers.h"
#include "gromacs/utility/cstringutil.h"
#include "gromacs/utility/exceptions.h"
#include "gromacs/utility/smalloc.h"
#include "gromacs/utility/stringutil.h"
#include "gromacs/utility/vec.h"

namespace gmx
{

namespace
{

/*! \brief
 * The structure and box that read() fills or that setStructure() points to.
 *
 * \tparam StructureInterface  StructureFile, or MultimodelStructureFile for
 *                             the formats that can hold several models.
 */
template<class StructureInterface>
class BuiltinStructureFile : public StructureInterface, public BoxContainingFile
{
public:
    explicit BuiltinStructureFile(std::filesystem::path path) : path_(std::move(path))
    {
        open_symtab(&symtab_);
        init_t_atoms(&readAtoms_, 0, false);
    }
    ~BuiltinStructureFile() override
    {
        done_atom(&readAtoms_);
        done_symtab(&symtab_);
    }

    int                  getNumAtoms() const override { return gmx::ssize(getPositions()); }
    ArrayRef<const RVec> getPositions() const override
    {
        return models_.empty() ? ArrayRef<const RVec>() : models_[activeModel_];
    }
    const t_atoms*       getAtoms() const override { return atoms_; }
    ArrayRef<const RVec> getVelocities() const override { return velocities_; }
    const std::string&   getTitle() const override { return title_; }
    void                 setStructure(const std::string&   title,
                                      const t_atoms&       atoms,
                                      ArrayRef<const RVec> positions,
                                      ArrayRef<const RVec> velocities) override
    {
        title_ = title;
        atoms_ = &atoms;
        models_.assign(1, positions);
        addedModels_.clear();
        activeModel_ = 0;
        velocities_  = velocities;
        selection_.reset();
    }
    void setSelection(ArrayRef<const int> selection) override { selection_ = selection; }

    std::optional<Matrix3x3> getBoundaries() const override { return box_; }
    void                     setBoundaries(const Matrix3x3& box) override { box_ = box; }
    PbcType                  getPbcType() const override { return pbcType_; }
    void                     setPbcType(PbcType pbcType) override { pbcType_ = pbcType; }

protected:
    //! Prepares the storage that a reader of \p numAtoms particles fills.
    void prepareForReading(int numAtoms, bool withPdbInfo)
    {
        done_atom(&readAtoms_);
        init_t_atoms(&readAtoms_, numAtoms, withPdbInfo);
        readPositions_.assign(numAtoms, RVec(0, 0, 0));
        readVelocities_.assign(numAtoms, RVec(0, 0, 0));
    }
    /*! \brief
     * Makes what a reader filled the current structure.
     *
     * Takes ownership of \p title, which the readers allocate. A box that
     * is all zero means the file has none.
     */
    void finishReading(char* title, bool haveVelocities, const matrix box, PbcType pbcType)
    {
        title_ = (title != nullptr) ? title : "";
        sfree(title);
        atoms_ = &readAtoms_;
        models_.assign(1, readPositions_);
        addedModels_.clear();
        activeModel_ = 0;
        velocities_ = haveVelocities ? ArrayRef<const RVec>(readVelocities_) : ArrayRef<const RVec>();
        selection_.reset();
        if (norm2(box[XX]) != 0 || norm2(box[YY]) != 0 || norm2(box[ZZ]) != 0)
        {
            box_ = createMatrix3x3FromLegacyMatrix(box);
        }
        else
        {
            box_.reset();
        }
        pbcType_ = pbcType;
    }

    //! The symbol table that a reader puts names in.
    t_symtab* symtabForReading() { return &symtab_; }
    //! The particles that a reader fills.
    t_atoms* atomsForReading() { return &readAtoms_; }
    //! The positions that a reader fills.
    rvec* positionsForReading() { return as_rvec_array(readPositions_.data()); }
    //! The velocities that a reader fills.
    rvec* velocitiesForReading() { return as_rvec_array(readVelocities_.data()); }

    //! The particles to be written; throws if setStructure() was not called.
    const t_atoms* atomsForWriting() const
    {
        if (atoms_ == nullptr)
        {
            GMX_THROW(APIError("Set a structure with setStructure() before writing it"));
        }
        return atoms_;
    }
    //! The title, as passed to the writers.
    const char* titleForWriting() const { return title_.c_str(); }
    //! The positions of \p model, as passed to the writers.
    const rvec* positionsForWriting(int model) const
    {
        return as_rvec_array(models_[model].data());
    }
    //! The velocities, as passed to the writers, or nullptr if there are none.
    const rvec* velocitiesForWriting() const
    {
        return velocities_.empty() ? nullptr : as_rvec_array(velocities_.data());
    }
    //! Whether only some of the particles are written.
    bool haveSelection() const { return selection_.has_value(); }
    //! The particles to be written, all of them when there is no selection.
    std::vector<int> particlesToWrite() const
    {
        if (haveSelection())
        {
            return { selection_->begin(), selection_->end() };
        }
        std::vector<int> all(atomsForWriting()->nr);
        std::iota(all.begin(), all.end(), 0);
        return all;
    }
    //! The box, as passed to the writers; all zero if there is none.
    void boxForWriting(matrix box) const
    {
        clear_mat(box);
        if (box_.has_value())
        {
            fillLegacyMatrix(box_.value(), box);
        }
    }
    //! The number of models.
    int numModels() const { return gmx::ssize(models_); }

    //! Writes every model with the PDB writer, as PQR if \p usePqrFormat.
    void writePdbModels(bool usePqrFormat) const
    {
        matrix box;
        boxForWriting(box);
        const std::vector<int> particles = particlesToWrite();
        IOHelper               out(path_, "w");
        for (int model = 0; model < numModels(); model++)
        {
            write_pdbfile_indexed(out.stream(),
                                  titleForWriting(),
                                  atomsForWriting(),
                                  positionsForWriting(model),
                                  pbcType_,
                                  box,
                                  ' ',
                                  numModels() > 1 ? model + 1 : -1,
                                  gmx::ssize(particles),
                                  particles.data(),
                                  nullptr,
                                  usePqrFormat);
        }
    }

    //! The path of the file.
    std::filesystem::path path_;
    //! The positions of each model.
    std::vector<ArrayRef<const RVec>> models_;
    //! Positions of the models added with addModel(), owned here.
    std::deque<std::vector<RVec>> addedModels_;
    //! The model that the StructureFile methods refer to.
    int activeModel_ = 0;

private:
    //! The symbol table that the names in readAtoms_ point into.
    t_symtab symtab_;
    //! The particles read from the file.
    t_atoms readAtoms_;
    //! The positions read from the file.
    std::vector<RVec> readPositions_;
    //! The velocities read from the file.
    std::vector<RVec> readVelocities_;
    //! The current particles: read from the file, or set with setStructure().
    const t_atoms* atoms_ = nullptr;
    //! The title.
    std::string title_;
    //! The current velocities.
    ArrayRef<const RVec> velocities_;
    //! The particles to be written, all of them when there is no selection.
    std::optional<ArrayRef<const int>> selection_;
    //! The box, if there is one.
    std::optional<Matrix3x3> box_;
    //! The periodic boundary type of the box.
    PbcType pbcType_ = PbcType::Unset;

    GMX_DISALLOW_COPY_MOVE_AND_ASSIGN(BuiltinStructureFile);
};

/*! \brief
 * A built-in format that can hold several models.
 *
 * read() reads the first model. Models added with addModel() are written
 * after it, and the velocities belong to the first model.
 */
class BuiltinMultimodelStructureFile : public BuiltinStructureFile<MultimodelStructureFile>
{
public:
    using BuiltinStructureFile::BuiltinStructureFile;

    int getNumModels() const override { return numModels(); }
    int addModel(ArrayRef<const RVec> positions) override
    {
        addedModels_.emplace_back(positions.begin(), positions.end());
        models_.emplace_back(addedModels_.back());
        return numModels() - 1;
    }
    void setActiveModel(int modelIndex) override
    {
        if (modelIndex < 0 || modelIndex >= numModels())
        {
            GMX_THROW(RangeError(formatString(
                    "Model %d does not exist in a file with %d models", modelIndex, numModels())));
        }
        activeModel_ = modelIndex;
    }
    int getActiveModel() const override { return activeModel_; }
};

/*! \brief
 * A built-in structure format whose files can also hold a sequence of frames.
 *
 * read() and write() handle the structure. The frames go through a stream
 * of their own, opened with openForReadingFrames() or
 * openForWritingFrames(). In append mode, frames are added to the end of
 * the file, and write() still replaces it.
 *
 * \tparam Base  The structure file class of the format.
 */
template<class Base>
class BuiltinTrajectoryStructureFile : public Base, public TrajectoryFile
{
public:
    BuiltinTrajectoryStructureFile(std::filesystem::path path, FileMode mode) :
        Base(std::move(path)), mode_(mode)
    {
    }
    void openForReadingFrames() override { frames_.open(this->path_, "r"); }
    void openForWritingFrames() override
    {
        frames_.open(this->path_, mode_ == FileMode::Append ? "a" : "w");
    }
    bool      supportsMode(FileMode /*mode*/) const override { return true; }
    t_fileio* legacyFileio() const override { return frames_.fileio(); }

protected:
    //! The stream of frames.
    FILE* frameStream() const { return frames_.stream(); }
    //! Opens the stream of frames again at the start of the file.
    void reopenFramesForReading() { frames_.open(this->path_, "r"); }
    //! Throws unless \p frame has positions, and particle names if \p needNames.
    void requireToWrite(const t_trxframe& frame, bool needNames) const
    {
        if (!frame.bX)
        {
            GMX_THROW(InvalidInputError(
                    formatString("Need coordinates to write a %s trajectory", this->formatName())));
        }
        if (needNames && !frame.bAtoms)
        {
            GMX_THROW(InvalidInputError(
                    formatString("Can not write a %s file without atom names", this->formatName())));
        }
    }
    //! The title that a frame is written with.
    static std::string frameTitle(const t_trxframe& frame)
    {
        return formatString("frame t= %.3f", frame.time);
    }

    //! The number of frames read or written.
    int numFrames_ = 0;

private:
    //! Whether frames replace the file or are added to it.
    FileMode mode_;
    //! The stream of frames, closed with the object.
    IOHelper frames_;
};

//! A GRO file. Several models or frames are written one after the other.
class GroFile : public BuiltinTrajectoryStructureFile<BuiltinMultimodelStructureFile>
{
public:
    using BuiltinTrajectoryStructureFile::BuiltinTrajectoryStructureFile;

    void read() override
    {
        int numAtoms = 0;
        get_coordnum(path_, &numAtoms);
        prepareForReading(numAtoms, false);
        char*      title = nullptr;
        matrix     box;
        const bool haveVelocities = gmx_gro_read_conf(path_,
                                                      symtabForReading(),
                                                      &title,
                                                      atomsForReading(),
                                                      positionsForReading(),
                                                      velocitiesForReading(),
                                                      box);
        finishReading(title, haveVelocities, box, PbcType::Unset);
    }
    void write() override
    {
        matrix box;
        boxForWriting(box);
        const std::vector<int> particles = particlesToWrite();
        IOHelper               out(path_, "w");
        for (int model = 0; model < numModels(); model++)
        {
            write_hconf_indexed_p(out.stream(),
                                  titleForWriting(),
                                  atomsForWriting(),
                                  gmx::ssize(particles),
                                  particles.data(),
                                  positionsForWriting(model),
                                  model == 0 ? velocitiesForWriting() : nullptr,
                                  box);
        }
    }
    const char* formatName() const override { return "GRO"; }

    bool readNextFrame(t_trxframe* frame, int /*flags*/) override
    {
        const bool isFirst = (numFrames_++ == 0);
        return isFirst ? gro_first_x_or_v(frameStream(), frame) > 0
                       : gro_next_x_or_v(frameStream(), frame);
    }
    void writeFrame(const t_trxframe& frame, ArrayRef<const int> index, gmx_conect_t* /*connections*/) override
    {
        requireToWrite(frame, true);
        const std::string title      = frameTitle(frame);
        const rvec*       velocities = frame.bV ? frame.v : nullptr;
        if (index.empty())
        {
            write_hconf_p(frameStream(), title.c_str(), frame.atoms, frame.x, velocities, frame.box);
        }
        else
        {
            write_hconf_indexed_p(frameStream(),
                                  title.c_str(),
                                  frame.atoms,
                                  index.ssize(),
                                  index.data(),
                                  frame.x,
                                  velocities,
                                  frame.box);
        }
    }
};

//! A G96 file.
class G96File : public BuiltinTrajectoryStructureFile<BuiltinStructureFile<StructureFile>>
{
public:
    using BuiltinTrajectoryStructureFile::BuiltinTrajectoryStructureFile;

    void read() override
    {
        t_trxframe frame;
        char       line[STRLEN + 1];
        frame.natoms = -1;
        frame.atoms  = nullptr;
        frame.x      = nullptr;
        frame.v      = nullptr;
        frame.f      = nullptr;
        IOHelper  in(path_, "r");
        const int numAtoms = read_g96_conf(in.stream(), path_, nullptr, &frame, nullptr, line);

        prepareForReading(numAtoms, false);
        char* title  = nullptr;
        frame.natoms = numAtoms;
        frame.atoms  = atomsForReading();
        frame.x      = positionsForReading();
        frame.v      = velocitiesForReading();
        frame.f      = nullptr;
        in.open(path_, "r");
        read_g96_conf(in.stream(), path_, &title, &frame, symtabForReading(), line);
        in.close();
        if (!frame.bBox)
        {
            clear_mat(frame.box);
        }
        finishReading(title, frame.bV, frame.box, PbcType::Unset);
    }
    void write() override
    {
        t_trxframe frame;
        clear_trxframe(&frame, TRUE);
        frame.natoms = atomsForWriting()->nr;
        frame.bAtoms = TRUE;
        frame.atoms  = const_cast<t_atoms*>(atomsForWriting());
        frame.bX     = TRUE;
        frame.x      = const_cast<rvec*>(positionsForWriting(0));
        if (velocitiesForWriting() != nullptr)
        {
            frame.bV = TRUE;
            frame.v  = const_cast<rvec*>(velocitiesForWriting());
        }
        frame.bBox = TRUE;
        boxForWriting(frame.box);
        const std::vector<int> particles = particlesToWrite();
        IOHelper               out(path_, "w");
        write_g96_conf(out.stream(),
                       titleForWriting(),
                       &frame,
                       haveSelection() ? gmx::ssize(particles) : -1,
                       haveSelection() ? particles.data() : nullptr);
    }
    const char* formatName() const override { return "G96"; }

    bool readNextFrame(t_trxframe* frame, int flags) override
    {
        if (numFrames_ == 0)
        {
            // Count the particles first; a compressed file cannot be
            // rewound, so it is opened again.
            persistentLine_.resize(STRLEN + 1);
            read_g96_conf(frameStream(), path_, nullptr, frame, nullptr, persistentLine_.data());
            reopenFramesForReading();
            clear_trxframe(frame, FALSE);
            if ((flags & (TRX_READ_X | TRX_NEED_X)) != 0)
            {
                snew(frame->x, frame->natoms);
            }
            if ((flags & (TRX_READ_V | TRX_NEED_V)) != 0)
            {
                snew(frame->v, frame->natoms);
            }
        }
        numFrames_++;
        read_g96_conf(frameStream(), {}, nullptr, frame, nullptr, persistentLine_.data());
        return frame->natoms > 0;
    }
    void writeFrame(const t_trxframe& frame, ArrayRef<const int> index, gmx_conect_t* /*connections*/) override
    {
        requireToWrite(frame, false);
        // write_trxframe() wrote G96 frames without a title, and
        // write_trxframe_indexed() with one.
        const std::string title = index.empty() ? std::string() : frameTitle(frame);
        write_g96_conf(frameStream(),
                       title.c_str(),
                       &frame,
                       index.empty() ? -1 : index.ssize(),
                       index.empty() ? nullptr : index.data());
    }

private:
    //! The line that the G96 reader keeps between frames.
    std::vector<char> persistentLine_;
};

/*! \brief
 * A PDB file, also used for BRK and ENT files.
 *
 * Several models are written with MODEL records. CONECT records are not
 * read yet, so the file reports no topology.
 */
class PdbFile : public BuiltinTrajectoryStructureFile<BuiltinMultimodelStructureFile>, public TopologyFile
{
public:
    using BuiltinTrajectoryStructureFile::BuiltinTrajectoryStructureFile;

    void read() override
    {
        int numAtoms = 0;
        {
            IOHelper in(path_, "r");
            get_pdb_coordnum(in.stream(), &numAtoms);
        }
        prepareForReading(numAtoms, fn2ftp(path_) == efPDB);
        char*   title = nullptr;
        matrix  box;
        PbcType pbcType = PbcType::Unset;
        gmx_pdb_read_conf(
                path_, symtabForReading(), &title, atomsForReading(), positionsForReading(), &pbcType, box);
        finishReading(title, false, box, pbcType);
    }
    void        write() override { writePdbModels(false); }
    const char* formatName() const override { return "PDB"; }

    bool readNextFrame(t_trxframe* frame, int /*flags*/) override
    {
        if (numFrames_ == 0)
        {
            fprintf(stderr, "Reading frames from pdb file");
            get_pdb_coordnum(frameStream(), &frame->natoms);
            if (frame->natoms == 0)
            {
                GMX_THROW(InvalidInputError("No coordinates in pdb file"));
            }
            frewind(frameStream());
            snew(frame->x, frame->natoms);
        }
        return readFrame(frame);
    }
    void writeFrame(const t_trxframe& frame, ArrayRef<const int> index, gmx_conect_t* connections) override
    {
        requireToWrite(frame, true);
        const std::string title   = frameTitle(frame);
        const PbcType     pbcType = frame.bPBC ? frame.pbcType : PbcType::Unset;
        if (index.empty())
        {
            write_pdbfile(
                    frameStream(), title.c_str(), frame.atoms, frame.x, pbcType, frame.box, ' ', frame.step, connections);
        }
        else
        {
            write_pdbfile_indexed(frameStream(),
                                  title.c_str(),
                                  frame.atoms,
                                  frame.x,
                                  pbcType,
                                  frame.box,
                                  ' ',
                                  frame.step,
                                  index.ssize(),
                                  index.data(),
                                  connections,
                                  false);
        }
    }

    TopologyCompleteness     getCompleteness() const override { return TopologyCompleteness::None; }
    std::vector<std::string> getMoleculeTypes() const override { return {}; }
    int getNumMolecules(const std::string& /*moleculeTypeName*/) const override { return 0; }

private:
    //! Reads the next model of \p frame->natoms particles into \p frame.
    bool readFrame(t_trxframe* frame)
    {
        t_atoms atoms;
        atoms.nr      = frame->natoms;
        atoms.atom    = nullptr;
        atoms.pdbinfo = nullptr;
        t_symtab symtab;
        open_symtab(&symtab);
        // A model number of -1 rather than NOTSET saves checking whether one was read.
        int       modelNumber = -1;
        char      title[STRLEN];
        PbcType   pbcType;
        matrix    box;
        const int numAtoms = read_pdbfile(
                frameStream(), title, &modelNumber, &atoms, &symtab, frame->x, &pbcType, box, nullptr);
        free_symtab(&symtab);
        setTrxFramePbcType(frame, pbcType);
        if (numFrames_ == 0)
        {
            fprintf(stderr, " '%s', %d atoms\n", title, frame->natoms);
        }
        frame->bPrec = TRUE;
        frame->prec  = 10000;
        frame->bX    = TRUE;
        frame->bBox  = (box[XX][XX] != 0.0);
        if (frame->bBox)
        {
            copy_mat(box, frame->box);
        }

        frame->step          = 0;
        const char* stepText = std::strstr(title, " step= ");
        frame->bStep = ((stepText != nullptr) && sscanf(stepText + 7, "%" SCNd64, &frame->step) == 1);
        double      time     = 0.0;
        const char* timeText = std::strstr(title, " t= ");
        frame->bTime         = ((timeText != nullptr) && sscanf(timeText + 4, "%lf", &time) == 1);
        frame->timeIsDouble  = true;
        frame->time          = time;

        if (numAtoms == 0)
        {
            return false;
        }
        if (numAtoms != frame->natoms)
        {
            GMX_THROW(InvalidInputError(formatString(
                    "Number of atoms in pdb frame %d is %d instead of %d", numFrames_, numAtoms, frame->natoms)));
        }
        numFrames_++;
        return true;
    }
};

/*! \brief
 * An mmCIF file.
 *
 * Several models are written with their model numbers. struct_conn records
 * are not read yet, so the file reports no topology.
 */
class CifFile : public BuiltinTrajectoryStructureFile<BuiltinMultimodelStructureFile>, public TopologyFile
{
public:
    using BuiltinTrajectoryStructureFile::BuiltinTrajectoryStructureFile;

    void read() override
    {
        int numAtoms = 0;
        CifIO::countAtoms(path_, &numAtoms);
        prepareForReading(numAtoms, true);
        char*   title = nullptr;
        matrix  box;
        PbcType pbcType = PbcType::Unset;
        CifIO::read(path_, symtabForReading(), &title, atomsForReading(), positionsForReading(), &pbcType, box);
        finishReading(title, false, box, pbcType);
    }
    void write() override
    {
        matrix box;
        boxForWriting(box);
        const std::vector<int> particles = particlesToWrite();
        if (numModels() == 1)
        {
            CifIO::writeIndexed(path_,
                                titleForWriting(),
                                atomsForWriting(),
                                positionsForWriting(0),
                                getPbcType(),
                                box,
                                gmx::ssize(particles),
                                particles.data());
            return;
        }
        IOHelper out(path_, "w");
        for (int model = 0; model < numModels(); model++)
        {
            CifIO::writeTrajectoryFrame(out.stream(),
                                        titleForWriting(),
                                        atomsForWriting(),
                                        positionsForWriting(model),
                                        box,
                                        gmx::ssize(particles),
                                        particles.data(),
                                        model + 1,
                                        model == 0);
        }
    }
    const char* formatName() const override { return "mmCIF"; }

    void openForReadingFrames() override
    {
        BuiltinTrajectoryStructureFile::openForReadingFrames();
        // The models come from the whole parsed file, one per frame.
        fprintf(stderr, "Reading frames from cif file");
        cifModels_ = CifIO::openTrajectory(path_);
    }
    bool readNextFrame(t_trxframe* frame, int /*flags*/) override
    {
        if (cifModels_ == nullptr)
        {
            frame->natoms = 0;
            return false;
        }
        if (numFrames_++ == 0)
        {
            frame->natoms = cifModels_->numAtoms();
            if (frame->natoms == 0)
            {
                GMX_THROW(InvalidInputError("No coordinates in cif file"));
            }
            snew(frame->x, frame->natoms);
        }
        return cifModels_->readNextFrame(frame);
    }
    void writeFrame(const t_trxframe& frame, ArrayRef<const int> index, gmx_conect_t* connections) override
    {
        requireToWrite(frame, true);
        std::vector<int> all;
        if (index.empty())
        {
            all.resize(frame.natoms);
            std::iota(all.begin(), all.end(), 0);
            index = all;
        }
        CifIO::writeTrajectoryFrame(frameStream(),
                                    nullptr,
                                    frame.atoms,
                                    frame.x,
                                    frame.box,
                                    index.ssize(),
                                    index.data(),
                                    numFrames_ + 1,
                                    numFrames_ == 0,
                                    connections);
        numFrames_++;
    }

    TopologyCompleteness     getCompleteness() const override { return TopologyCompleteness::None; }
    std::vector<std::string> getMoleculeTypes() const override { return {}; }
    int getNumMolecules(const std::string& /*moleculeTypeName*/) const override { return 0; }

private:
    //! The parsed models being read, or nullptr if the file has none.
    std::unique_ptr<CifTrajectory> cifModels_;
};

//! An Espresso file.
class EspFile : public BuiltinStructureFile<StructureFile>
{
public:
    using BuiltinStructureFile::BuiltinStructureFile;

    void read() override
    {
        const int numAtoms = get_espresso_coordnum(path_);
        prepareForReading(numAtoms, false);
        char*      title = nullptr;
        matrix     box;
        const bool haveVelocities = gmx_espresso_read_conf(path_,
                                                           symtabForReading(),
                                                           &title,
                                                           atomsForReading(),
                                                           positionsForReading(),
                                                           velocitiesForReading(),
                                                           box);
        finishReading(title, haveVelocities, box, PbcType::Unset);
    }
    void write() override
    {
        matrix box;
        boxForWriting(box);
        const std::vector<int> particles = particlesToWrite();
        IOHelper               out(path_, "w");
        write_espresso_conf_indexed(out.stream(),
                                    titleForWriting(),
                                    atomsForWriting(),
                                    gmx::ssize(particles),
                                    haveSelection() ? particles.data() : nullptr,
                                    positionsForWriting(0),
                                    velocitiesForWriting(),
                                    box);
    }
    bool supportsMode(FileMode mode) const override
    {
        return mode == FileMode::Read || mode == FileMode::Write;
    }
    const char* formatName() const override { return "Espresso"; }
};

//! A PQR file, which is written with the PDB writer and cannot be read.
class PqrFile : public BuiltinStructureFile<StructureFile>
{
public:
    using BuiltinStructureFile::BuiltinStructureFile;

    [[noreturn]] void read() override
    {
        GMX_THROW(NotImplementedError("PQR files cannot be read"));
    }
    void        write() override { writePdbModels(true); }
    bool        supportsMode(FileMode mode) const override { return mode == FileMode::Write; }
    const char* formatName() const override { return "PQR"; }
};

} // namespace

FilePointer createGroFile(const std::filesystem::path& path, FileMode mode)
{
    return std::make_unique<GroFile>(path, mode);
}

FilePointer createG96File(const std::filesystem::path& path, FileMode mode)
{
    return std::make_unique<G96File>(path, mode);
}

FilePointer createPdbFile(const std::filesystem::path& path, FileMode mode)
{
    return std::make_unique<PdbFile>(path, mode);
}

FilePointer createCifFile(const std::filesystem::path& path, FileMode mode)
{
    return std::make_unique<CifFile>(path, mode);
}

FilePointer createEspFile(const std::filesystem::path& path, FileMode /*mode*/)
{
    return std::make_unique<EspFile>(path);
}

FilePointer createPqrFile(const std::filesystem::path& path, FileMode /*mode*/)
{
    return std::make_unique<PqrFile>(path);
}

} // namespace gmx
