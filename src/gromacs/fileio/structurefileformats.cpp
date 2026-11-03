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
 * and writes.
 *
 * Each format reads and writes with the existing functions of groio.h,
 * g96io.h, pdbio.h, cifio.h and espio.h, so its files are the same as before.
 *
 * \author Amr Alhossary <amr_alhossary@hotmail.com>
 * \ingroup module_fileio
 */
#include "gmxpre.h"

#include "structurefileformats.h"

#include <cstdio>

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
#include "gromacs/fileio/pdbio.h"
#include "gromacs/fileio/structurefile.h"
#include "gromacs/fileio/topologyfile.h"
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
        FILE*                  out       = gmx_fio_fopen(path_, "w");
        for (int model = 0; model < numModels(); model++)
        {
            write_pdbfile_indexed(out,
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
        gmx_fio_fclose(out);
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

//! A GRO file. Several models are written one after the other.
class GroFile : public BuiltinMultimodelStructureFile
{
public:
    using BuiltinMultimodelStructureFile::BuiltinMultimodelStructureFile;

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
        FILE*                  out       = gmx_fio_fopen(path_, "w");
        for (int model = 0; model < numModels(); model++)
        {
            write_hconf_indexed_p(out,
                                  titleForWriting(),
                                  atomsForWriting(),
                                  gmx::ssize(particles),
                                  particles.data(),
                                  positionsForWriting(model),
                                  model == 0 ? velocitiesForWriting() : nullptr,
                                  box);
        }
        gmx_fio_fclose(out);
    }
    bool supportsMode(FileMode mode) const override
    {
        return mode == FileMode::Read || mode == FileMode::Write;
    }
    const char* formatName() const override { return "GRO"; }
};

//! A G96 file.
class G96File : public BuiltinStructureFile<StructureFile>
{
public:
    using BuiltinStructureFile::BuiltinStructureFile;

    void read() override
    {
        t_trxframe frame;
        char       line[STRLEN + 1];
        frame.natoms       = -1;
        frame.atoms        = nullptr;
        frame.x            = nullptr;
        frame.v            = nullptr;
        frame.f            = nullptr;
        FILE*     in       = gmx_fio_fopen(path_, "r");
        const int numAtoms = read_g96_conf(in, path_, nullptr, &frame, nullptr, line);
        gmx_fio_fclose(in);

        prepareForReading(numAtoms, false);
        char* title  = nullptr;
        frame.natoms = numAtoms;
        frame.atoms  = atomsForReading();
        frame.x      = positionsForReading();
        frame.v      = velocitiesForReading();
        frame.f      = nullptr;
        in           = gmx_fio_fopen(path_, "r");
        read_g96_conf(in, path_, &title, &frame, symtabForReading(), line);
        gmx_fio_fclose(in);
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
        FILE*                  out       = gmx_fio_fopen(path_, "w");
        write_g96_conf(out,
                       titleForWriting(),
                       &frame,
                       haveSelection() ? gmx::ssize(particles) : -1,
                       haveSelection() ? particles.data() : nullptr);
        gmx_fio_fclose(out);
    }
    bool supportsMode(FileMode mode) const override
    {
        return mode == FileMode::Read || mode == FileMode::Write;
    }
    const char* formatName() const override { return "G96"; }
};

/*! \brief
 * A PDB file, also used for BRK and ENT files.
 *
 * Several models are written with MODEL records. CONECT records are not
 * read yet, so the file reports no topology.
 */
class PdbFile : public BuiltinMultimodelStructureFile, public TopologyFile
{
public:
    using BuiltinMultimodelStructureFile::BuiltinMultimodelStructureFile;

    void read() override
    {
        FILE* in       = gmx_fio_fopen(path_, "r");
        int   numAtoms = 0;
        get_pdb_coordnum(in, &numAtoms);
        gmx_fio_fclose(in);
        prepareForReading(numAtoms, fn2ftp(path_) == efPDB);
        char*   title = nullptr;
        matrix  box;
        PbcType pbcType = PbcType::Unset;
        gmx_pdb_read_conf(
                path_, symtabForReading(), &title, atomsForReading(), positionsForReading(), &pbcType, box);
        finishReading(title, false, box, pbcType);
    }
    void write() override { writePdbModels(false); }
    bool supportsMode(FileMode mode) const override
    {
        return mode == FileMode::Read || mode == FileMode::Write;
    }
    const char* formatName() const override { return "PDB"; }

    TopologyCompleteness     getCompleteness() const override { return TopologyCompleteness::None; }
    std::vector<std::string> getMoleculeTypes() const override { return {}; }
    int getNumMolecules(const std::string& /*moleculeTypeName*/) const override { return 0; }
};

/*! \brief
 * An mmCIF file.
 *
 * Several models are written with their model numbers. struct_conn records
 * are not read yet, so the file reports no topology.
 */
class CifFile : public BuiltinMultimodelStructureFile, public TopologyFile
{
public:
    using BuiltinMultimodelStructureFile::BuiltinMultimodelStructureFile;

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
        FILE* out = gmx_fio_fopen(path_, "w");
        for (int model = 0; model < numModels(); model++)
        {
            CifIO::writeTrajectoryFrame(out,
                                        titleForWriting(),
                                        atomsForWriting(),
                                        positionsForWriting(model),
                                        box,
                                        gmx::ssize(particles),
                                        particles.data(),
                                        model + 1,
                                        model == 0);
        }
        gmx_fio_fclose(out);
    }
    bool supportsMode(FileMode mode) const override
    {
        return mode == FileMode::Read || mode == FileMode::Write;
    }
    const char* formatName() const override { return "mmCIF"; }

    TopologyCompleteness     getCompleteness() const override { return TopologyCompleteness::None; }
    std::vector<std::string> getMoleculeTypes() const override { return {}; }
    int getNumMolecules(const std::string& /*moleculeTypeName*/) const override { return 0; }
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
        FILE*                  out       = gmx_fio_fopen(path_, "w");
        write_espresso_conf_indexed(out,
                                    titleForWriting(),
                                    atomsForWriting(),
                                    gmx::ssize(particles),
                                    haveSelection() ? particles.data() : nullptr,
                                    positionsForWriting(0),
                                    velocitiesForWriting(),
                                    box);
        gmx_fio_fclose(out);
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

    void        read() override { GMX_THROW(NotImplementedError("PQR files cannot be read")); }
    void        write() override { writePdbModels(true); }
    bool        supportsMode(FileMode mode) const override { return mode == FileMode::Write; }
    const char* formatName() const override { return "PQR"; }
};

} // namespace

FilePointer createGroFile(const std::filesystem::path& path, FileMode /*mode*/)
{
    return std::make_unique<GroFile>(path);
}

FilePointer createG96File(const std::filesystem::path& path, FileMode /*mode*/)
{
    return std::make_unique<G96File>(path);
}

FilePointer createPdbFile(const std::filesystem::path& path, FileMode /*mode*/)
{
    return std::make_unique<PdbFile>(path);
}

FilePointer createCifFile(const std::filesystem::path& path, FileMode /*mode*/)
{
    return std::make_unique<CifFile>(path);
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
