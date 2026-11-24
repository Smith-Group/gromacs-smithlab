/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright 1991- The GROMACS Authors
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
#include "gmxpre.h"

#include "gromacs/fileio/confio.h"

#include <cstdio>
#include <cstring>

#include <algorithm>
#include <filesystem>
#include <memory>
#include <optional>
#include <string>
#include <vector>

#include "gromacs/fileio/boxcontainingfile.h"
#include "gromacs/fileio/fileformatregistry.h"
#include "gromacs/fileio/filetypes.h"
#include "gromacs/fileio/gmxfio.h"
#include "gromacs/fileio/groio.h"
#include "gromacs/fileio/structurefile.h"
#include "gromacs/fileio/tpxio.h"
#include "gromacs/mdtypes/md_enums.h"
#include "gromacs/topology/atoms.h"
#include "gromacs/topology/block.h"
#include "gromacs/topology/chainid.h"
#include "gromacs/topology/mtop_util.h"
#include "gromacs/topology/symtab.h"
#include "gromacs/topology/topology.h"
#include "gromacs/utility/arrayref.h"
#include "gromacs/utility/basedefinitions.h"
#include "gromacs/utility/cstringutil.h"
#include "gromacs/utility/exceptions.h"
#include "gromacs/utility/fatalerror.h"
#include "gromacs/utility/gmxassert.h"
#include "gromacs/utility/matrix.h"
#include "gromacs/utility/smalloc.h"
#include "gromacs/utility/stringutil.h"
#include "gromacs/utility/vec.h"
#include "gromacs/utility/vectypes.h"

namespace
{

//! Writes a structure with the file object of the format of \p outfile.
void writeStructure(const std::filesystem::path&            outfile,
                    const char*                             title,
                    const t_atoms*                          atoms,
                    const rvec                              x[],
                    const rvec*                             v,
                    PbcType                                 pbcType,
                    const matrix                            box,
                    std::optional<gmx::ArrayRef<const int>> selection)
{
    std::unique_ptr<gmx::StructureFile> file =
            gmx::openFileAs<gmx::StructureFile>(outfile, gmx::FileMode::Write);
    file->setStructure((title != nullptr) ? title : "",
                       *atoms,
                       gmx::constArrayRefFromArray(reinterpret_cast<const gmx::RVec*>(x), atoms->nr),
                       (v != nullptr) ? gmx::constArrayRefFromArray(
                                                reinterpret_cast<const gmx::RVec*>(v), atoms->nr)
                                      : gmx::ArrayRef<const gmx::RVec>());
    if (selection.has_value())
    {
        file->setSelection(selection.value());
    }
    if (auto* boxFile = dynamic_cast<gmx::BoxContainingFile*>(file.get()))
    {
        if (box != nullptr)
        {
            boxFile->setBoundaries(gmx::createMatrix3x3FromLegacyMatrix(box));
        }
        boxFile->setPbcType(pbcType);
    }
    file->write();
}

} // namespace

void write_sto_conf_indexed(const std::filesystem::path& outfile,
                            const char*                  title,
                            const t_atoms*               atoms,
                            const rvec                   x[],
                            const rvec*                  v,
                            PbcType                      pbcType,
                            const matrix                 box,
                            int                          nindex,
                            int                          index[])
{
    writeStructure(outfile, title, atoms, x, v, pbcType, box, gmx::constArrayRefFromArray(index, nindex));
}

void write_sto_conf(const std::filesystem::path& outfile,
                    const char*                  title,
                    const t_atoms*               atoms,
                    const rvec                   x[],
                    const rvec*                  v,
                    PbcType                      pbcType,
                    const matrix                 box)
{
    writeStructure(outfile, title, atoms, x, v, pbcType, box, std::nullopt);
}

void write_sto_conf_mtop(const std::filesystem::path& outfile,
                         const char*                  title,
                         const gmx_mtop_t&            mtop,
                         const rvec                   x[],
                         const rvec*                  v,
                         PbcType                      pbcType,
                         const matrix                 box)
{
    if (fn2ftp(outfile) == efGRO)
    {
        /* The GRO writer reads mtop directly, without the global atoms. */
        FILE* out = gmx_fio_fopen(outfile, "w");
        write_hconf_mtop(out, title, mtop, x, v, box);
        gmx_fio_fclose(out);
    }
    else
    {
        /* This is a brute force approach which requires a lot of memory.
         * We should implement mtop versions of all writing routines.
         */
        t_atoms atoms = gmx_mtop_global_atoms(mtop);

        write_sto_conf(outfile, title, &atoms, x, v, pbcType, box);

        done_atom(&atoms);
    }
}

/*! \brief
 * Copies \p source into \p atoms, with every name put in \p symtab.
 *
 * Names are put in the order the readers put them: the name of a residue
 * before the names of its particles.
 */
static void copyAtoms(const t_atoms& source, t_symtab* symtab, t_atoms* atoms)
{
    init_t_atoms(atoms, source.nr, source.pdbinfo != nullptr);
    atoms->nres = source.nres;
    std::vector<bool> residueCopied(source.nres, false);
    auto              copyResidue = [&](int residue)
    {
        atoms->resinfo[residue]      = source.resinfo[residue];
        atoms->resinfo[residue].name = put_symtab(symtab, *source.resinfo[residue].name);
        if (source.resinfo[residue].rtp != nullptr)
        {
            atoms->resinfo[residue].rtp = put_symtab(symtab, *source.resinfo[residue].rtp);
        }
        residueCopied[residue] = true;
    };
    if (source.atomtype != nullptr)
    {
        snew(atoms->atomtype, source.nr);
    }
    if (source.atomtypeB != nullptr)
    {
        snew(atoms->atomtypeB, source.nr);
    }
    for (int i = 0; i < source.nr; i++)
    {
        atoms->atom[i]    = source.atom[i];
        const int residue = source.atom[i].resind;
        if (residue < source.nres && !residueCopied[residue])
        {
            copyResidue(residue);
        }
        atoms->atomname[i] = put_symtab(symtab, *source.atomname[i]);
        if (source.atomtype != nullptr)
        {
            atoms->atomtype[i] = put_symtab(symtab, *source.atomtype[i]);
        }
        if (source.atomtypeB != nullptr)
        {
            atoms->atomtypeB[i] = put_symtab(symtab, *source.atomtypeB[i]);
        }
        if (source.pdbinfo != nullptr)
        {
            atoms->pdbinfo[i] = source.pdbinfo[i];
        }
    }
    for (int residue = 0; residue < source.nres; residue++)
    {
        if (!residueCopied[residue])
        {
            copyResidue(residue);
        }
    }
    atoms->haveMass    = source.haveMass;
    atoms->haveCharge  = source.haveCharge;
    atoms->haveType    = source.haveType;
    atoms->haveBState  = source.haveBState;
    atoms->havePdbInfo = source.havePdbInfo;
}

void readConfAndAtoms(const std::filesystem::path& infile,
                      t_symtab*                    symtab,
                      char**                       name,
                      t_atoms*                     atoms,
                      PbcType*                     pbcType,
                      rvec**                       x,
                      rvec**                       v,
                      matrix                       box)
{
    if (fn2ftp(infile) == efTPR)
    {
        bool       haveTopology;
        gmx_mtop_t mtop;
        readConfAndTopology(infile, &haveTopology, &mtop, pbcType, x, v, box);
        *symtab                          = mtop.symtab;
        *name                            = gmx_strdup(*mtop.name);
        *atoms                           = gmx_mtop_global_atoms(mtop);
        gmx::RangePartitioning molecules = gmx_mtop_molecules(mtop);
        gmx::makeChainIdentifiersAfterTprReading(atoms, molecules);

        /* Inelegant solution to avoid all char pointers in atoms becoming
         * invalid after destruction of mtop.
         * This will be fixed soon by converting t_symtab to C++.
         */
        mtop.symtab.symbuf = nullptr;
        mtop.symtab.nr     = 0;

        return;
    }

    std::unique_ptr<gmx::StructureFile> file =
            gmx::openFileAs<gmx::StructureFile>(infile, gmx::FileMode::Read);
    file->read();
    const int natoms = file->getNumAtoms();
    if (natoms == 0)
    {
        fprintf(stderr, "Warning: Number of atoms in %s is 0\n", infile.string().c_str());
    }
    if (file->getAtoms() == nullptr)
    {
        GMX_THROW(gmx::InvalidInputError(
                gmx::formatString("%s files carry no particle names, so %s cannot be used here",
                                  file->formatName(),
                                  infile.string().c_str())));
    }
    copyAtoms(*file->getAtoms(), symtab, atoms);
    if (name != nullptr)
    {
        *name = gmx_strdup(file->getTitle().c_str());
    }
    if (x != nullptr)
    {
        snew(*x, natoms);
        std::copy(file->getPositions().begin(),
                  file->getPositions().end(),
                  reinterpret_cast<gmx::RVec*>(*x));
    }
    if (v != nullptr)
    {
        snew(*v, natoms);
        std::copy(file->getVelocities().begin(),
                  file->getVelocities().end(),
                  reinterpret_cast<gmx::RVec*>(*v));
    }
    clear_mat(box);
    if (pbcType != nullptr)
    {
        *pbcType = PbcType::Unset;
    }
    if (const auto* boxFile = dynamic_cast<const gmx::BoxContainingFile*>(file.get()))
    {
        if (boxFile->getBoundaries().has_value())
        {
            gmx::fillLegacyMatrix(boxFile->getBoundaries().value(), box);
        }
        if (pbcType != nullptr)
        {
            *pbcType = boxFile->getPbcType();
        }
    }
}

void readConfAndTopology(const std::filesystem::path& infile,
                         bool*                        haveTopology,
                         gmx_mtop_t*                  mtop,
                         PbcType*                     pbcType,
                         rvec**                       x,
                         rvec**                       v,
                         matrix                       box)
{
    GMX_RELEASE_ASSERT(mtop != nullptr, "readConfAndTopology requires mtop!=NULL");

    if (pbcType != nullptr)
    {
        *pbcType = PbcType::Unset;
    }

    *haveTopology = fn2bTPX(infile);
    if (*haveTopology)
    {
        TpxFileHeader header = readTpxHeader(infile, true);
        if (x)
        {
            snew(*x, header.natoms);
        }
        if (v)
        {
            snew(*v, header.natoms);
        }
        int     natoms;
        PbcType pbcType_tmp = read_tpx(
                infile, nullptr, box, &natoms, (x == nullptr) ? nullptr : *x, (v == nullptr) ? nullptr : *v, mtop);
        if (pbcType != nullptr)
        {
            *pbcType = pbcType_tmp;
        }
    }
    else
    {
        t_symtab symtab;
        char*    name;
        t_atoms  atoms;

        open_symtab(&symtab);

        readConfAndAtoms(infile, &symtab, &name, &atoms, pbcType, x, v, box);

        convertAtomsToMtop(&symtab, put_symtab(&symtab, name), &atoms, mtop);
        sfree(name);
    }
}

gmx_bool read_tps_conf(const std::filesystem::path& infile,
                       t_topology*                  top,
                       PbcType*                     pbcType,
                       rvec**                       x,
                       rvec**                       v,
                       matrix                       box,
                       gmx_bool                     requireMasses)
{
    bool       haveTopology;
    gmx_mtop_t mtop;

    readConfAndTopology(infile, &haveTopology, &mtop, pbcType, x, v, box);

    *top = gmx_mtop_t_to_t_topology(&mtop, true);

    // We can only generate chain identifiers when we have a list of molecules (haveTopology).
    // This conditional has the side effect of preserving chain identifiers read from a PDB file.
    if (haveTopology)
    {
        gmx::tpx_make_chain_identifiers(&top->atoms, &top->mols);
    }

    if (requireMasses && !top->atoms.haveMass)
    {
        atomsSetMassesBasedOnNames(&top->atoms, TRUE);

        if (!top->atoms.haveMass)
        {
            gmx_fatal(FARGS,
                      "Masses were requested, but for some atom(s) masses could not be found in "
                      "the database. Use a tpr file as input, if possible, or add these atoms to "
                      "the mass database.");
        }
    }

    return haveTopology;
}
