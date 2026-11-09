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
 * Tests for gmx::FileFormatRegistry and the capability interfaces.
 *
 * \author Amr Alhossary <amr_alhossary@hotmail.com>
 * \ingroup module_fileio
 */
#include "gmxpre.h"

#include "gromacs/fileio/fileformatregistry.h"

#include <algorithm>
#include <filesystem>
#include <fstream>
#include <optional>
#include <set>
#include <string>
#include <utility>
#include <vector>

#include <gmock/gmock.h>
#include <gtest/gtest.h>

#include "gromacs/fileio/boxcontainingfile.h"
#include "gromacs/fileio/confio.h"
#include "gromacs/fileio/filetypes.h"
#include "gromacs/fileio/structurefile.h"
#include "gromacs/fileio/topologyfile.h"
#include "gromacs/mdtypes/md_enums.h"
#include "gromacs/pbcutil/pbc.h"
#include "gromacs/topology/atoms.h"
#include "gromacs/topology/symtab.h"
#include "gromacs/utility/exceptions.h"
#include "gromacs/utility/stringutil.h"

#include "testutils/testasserts.h"
#include "testutils/testfilemanager.h"

namespace gmx
{
namespace test
{
namespace
{

/*! \brief A format defined in this file, as a format outside fileio would be.
 *
 * It stores one position per line and implements two capability interfaces.
 * Registering it edits no file of the fileio module.
 */
class XyzTestFile : public StructureFile, public BoxContainingFile
{
public:
    XyzTestFile(std::filesystem::path path, FileMode /*mode*/) : path_(std::move(path)) {}

    void read() override
    {
        std::ifstream in(path_);
        if (!in)
        {
            GMX_THROW(FileIOError("Cannot open " + path_.string()));
        }
        positions_.clear();
        real x = 0, y = 0, z = 0;
        while (in >> x >> y >> z)
        {
            positions_.emplace_back(x, y, z);
        }
    }
    void write() override
    {
        std::ofstream out(path_);
        for (int i = 0; i < gmx::ssize(positions_); i++)
        {
            if (selection_.has_value()
                && std::find(selection_->begin(), selection_->end(), i) == selection_->end())
            {
                continue;
            }
            out << positions_[i][XX] << ' ' << positions_[i][YY] << ' ' << positions_[i][ZZ] << '\n';
        }
    }
    bool        supportsMode(FileMode mode) const override { return mode != FileMode::Append; }
    const char* formatName() const override { return "XYZ test file"; }

    int                  getNumAtoms() const override { return gmx::ssize(positions_); }
    ArrayRef<const RVec> getPositions() const override { return positions_; }
    const t_atoms*       getAtoms() const override { return nullptr; }
    ArrayRef<const RVec> getVelocities() const override { return {}; }
    const std::string&   getTitle() const override { return title_; }
    void                 setStructure(const std::string& title,
                                      const t_atoms& /*atoms*/,
                                      ArrayRef<const RVec> positions,
                                      ArrayRef<const RVec> /*velocities*/) override
    {
        title_ = title;
        positions_.assign(positions.begin(), positions.end());
        selection_.reset();
    }
    void setSelection(ArrayRef<const int> selection) override { selection_ = selection; }

    std::optional<Matrix3x3> getBoundaries() const override { return box_; }
    void                     setBoundaries(const Matrix3x3& box) override { box_ = box; }
    PbcType getPbcType() const override { return box_.has_value() ? PbcType::Xyz : PbcType::No; }
    void    setPbcType(PbcType /*pbcType*/) override {}

private:
    std::filesystem::path              path_;
    std::string                        title_;
    std::vector<RVec>                  positions_;
    std::optional<ArrayRef<const int>> selection_;
    std::optional<Matrix3x3>           box_;
};

//! Creates an XyzTestFile.
FilePointer createXyzTestFile(const std::filesystem::path& path, FileMode mode)
{
    return std::make_unique<XyzTestFile>(path, mode);
}

//! Registers the test format, which states its capabilities correctly.
const FileFormatRegistrar c_xyzTestRegistrar({ ".gmxtestxyz",
                                               "XYZ test file",
                                               FileCapability::Structure | FileCapability::Box,
                                               true,
                                               true,
                                               createXyzTestFile,
                                               std::nullopt });

//! Registers the same class as a format that declares an interface it does not implement.
const FileFormatRegistrar c_wrongTestRegistrar({ ".gmxtestwrong",
                                                 "Wrongly described test file",
                                                 FileCapability::Structure | FileCapability::Topology,
                                                 true,
                                                 true,
                                                 createXyzTestFile,
                                                 std::nullopt });

//! Returns the GromacsFileType of every registered format providing \p required that can be used in \p mode.
std::set<int> legacyTypesProviding(FileCapability required, FileMode mode)
{
    std::set<int> types;
    for (const FileFormatInfo* format : FileFormatRegistry::instance().formatsProviding(required))
    {
        const bool canUseMode = (mode == FileMode::Read) ? format->canRead : format->canWrite;
        if (canUseMode && format->legacyType.has_value())
        {
            types.insert(format->legacyType.value());
        }
    }
    return types;
}

//! Returns the members of the generic file type \p genericType from filetypes.cpp.
std::set<int> genericTypeMembers(int genericType)
{
    const int* list = ftp2generic_list(genericType);
    return std::set<int>(list, list + ftp2generic_count(genericType));
}

TEST(FileFormatRegistryTest, TrajectoryFormatsMatchGenericTypes)
{
    // BRK and ENT files are PDB files, so trxio.h reads and writes their
    // frames, but efTRX and efTRO do not offer them.
    std::set<int> readable = genericTypeMembers(efTRX);
    readable.insert({ efBRK, efENT });
    EXPECT_EQ(legacyTypesProviding(FileCapability::Trajectory, FileMode::Read), readable);
    std::set<int> writable = genericTypeMembers(efTRO);
    writable.insert({ efBRK, efENT });
    EXPECT_EQ(legacyTypesProviding(FileCapability::Trajectory, FileMode::Write), writable);
}

TEST(FileFormatRegistryTest, StructureFormatsMatchGenericTypes)
{
    EXPECT_EQ(legacyTypesProviding(FileCapability::Structure, FileMode::Read), genericTypeMembers(efSTX));
    // write_sto_conf_indexed() writes PQR, but efSTO does not offer it as a
    // generic structure output.
    std::set<int> writable = genericTypeMembers(efSTO);
    writable.insert(efPQR);
    EXPECT_EQ(legacyTypesProviding(FileCapability::Structure, FileMode::Write), writable);
}

TEST(FileFormatRegistryTest, VelocityFormatsMatchTrjconvTest)
{
    // The list in TrjconvFormatToFormatTest::outputFormatSupportsVelocities()
    // in src/gromacs/tools/tests/trjconv.cpp.
    const std::set<int> trjconvList = { efGRO, efTRR, efH5MD, efTNG, efG96 };
    EXPECT_EQ(legacyTypesProviding(FileCapability::Trajectory | FileCapability::Velocities, FileMode::Write),
              trjconvList);
}

TEST(FileFormatRegistryTest, LegacyFormatsAgreeWithFiletypes)
{
    for (const FileFormatInfo* format : FileFormatRegistry::instance().allFormats())
    {
        if (!format->legacyType.has_value())
        {
            continue;
        }
        SCOPED_TRACE(format->extension);
        EXPECT_EQ(fn2ftp(format->extension), format->legacyType.value());
        EXPECT_EQ(FileFormatRegistry::instance().byLegacyType(format->legacyType.value()), format);
    }
}

TEST(FileFormatRegistryTest, FindsFormatByExtensionIgnoringCase)
{
    const FileFormatInfo* format = FileFormatRegistry::instance().formatOf("conf.GRO");
    ASSERT_NE(format, nullptr);
    EXPECT_EQ(format->legacyType, efGRO);
    EXPECT_EQ(FileFormatRegistry::instance().formatOf("conf.unknownextension"), nullptr);
}

TEST(FileFormatRegistryTest, RefusesDuplicateExtension)
{
    FileFormatInfo duplicate;
    duplicate.extension = ".GRO";
    EXPECT_THROW_GMX(FileFormatRegistry::instance().registerFormat(duplicate), InternalError);
}

TEST(FileFormatRegistryTest, RefusesUnknownExtension)
{
    EXPECT_THROW_GMX(FileFormatRegistry::instance().open("conf.unknownextension", FileMode::Read),
                     InvalidInputError);
}

TEST(FileFormatRegistryTest, RefusesModeTheFormatCannotDo)
{
    EXPECT_THROW_GMX(FileFormatRegistry::instance().open("topol.tpr", FileMode::Write), InvalidInputError);
    EXPECT_THROW_GMX(FileFormatRegistry::instance().open("conf.pqr", FileMode::Read), InvalidInputError);
    EXPECT_THROW_GMX(FileFormatRegistry::instance().open("conf.gmxtestxyz", FileMode::Append),
                     InvalidInputError);
}

TEST(FileFormatRegistryTest, LegacyFormatIsNotYetAFileObject)
{
    EXPECT_THROW_GMX(FileFormatRegistry::instance().open("topol.tpr", FileMode::Read), NotImplementedError);
}

TEST(FileFormatRegistryTest, RefusesMissingCapabilityBeforeOpening)
{
    // Neither file exists: the request is refused from the format alone.
    EXPECT_THROW_GMX(openFileAs<StructureFile>("topol.top", FileMode::Read), InvalidInputError);
    EXPECT_THROW_GMX(openFileAs<TopologyFile>("conf.gmxtestxyz", FileMode::Read), InvalidInputError);
}

TEST(FileFormatRegistryTest, DetectsFormatDeclaringInterfaceItDoesNotImplement)
{
    EXPECT_THROW_GMX(FileFormatRegistry::instance().open("conf.gmxtestwrong", FileMode::Read), InternalError);
}

TEST(FileFormatRegistryTest, OpensRegisteredFormatWithoutDoingIO)
{
    // The file does not exist; creating the object does not touch it.
    std::unique_ptr<StructureFile> file = openFileAs<StructureFile>("missing.gmxtestxyz", FileMode::Read);
    ASSERT_NE(file, nullptr);
    EXPECT_STREQ(file->formatName(), "XYZ test file");
    EXPECT_NE(dynamic_cast<BoxContainingFile*>(file.get()), nullptr);
    EXPECT_EQ(dynamic_cast<TopologyFile*>(file.get()), nullptr);
    EXPECT_THROW_GMX(file->read(), FileIOError);
}

TEST(FileFormatRegistryTest, WritesAndReadsThroughCapabilityInterfaces)
{
    TestFileManager             fileManager;
    const std::filesystem::path path = fileManager.getTemporaryFilePath("coordinates.gmxtestxyz");
    const std::vector<RVec>     positions = { { 0.5, 1.0, 1.5 }, { 2.0, 2.5, 3.0 } };
    {
        t_atoms atoms;
        init_t_atoms(&atoms, gmx::ssize(positions), false);
        std::unique_ptr<StructureFile> file = openFileAs<StructureFile>(path, FileMode::Write);
        file->setStructure("", atoms, positions, {});
        file->write();
        done_atom(&atoms);
    }
    std::unique_ptr<StructureFile> file = openFileAs<StructureFile>(path, FileMode::Read);
    file->read();
    ASSERT_EQ(file->getNumAtoms(), 2);
    EXPECT_REAL_EQ(3.0, file->getPositions()[1][ZZ]);
    EXPECT_FALSE(dynamic_cast<BoxContainingFile&>(*file).getBoundaries().has_value());
}

TEST(FileFormatRegistryTest, ConfioWritesFormatRegisteredOutsideFileio)
{
    TestFileManager             fileManager;
    const std::filesystem::path path = fileManager.getTemporaryFilePath("selected.gmxtestxyz");
    const rvec                  x[3] = { { 0.5, 1.0, 1.5 }, { 2.0, 2.5, 3.0 }, { 3.5, 4.0, 4.5 } };
    matrix                      box  = { { 5, 0, 0 }, { 0, 5, 0 }, { 0, 0, 5 } };
    int                         index[2] = { 2, 0 };
    t_atoms                     atoms;
    init_t_atoms(&atoms, 3, false);
    write_sto_conf_indexed(path, "title", &atoms, x, nullptr, PbcType::Xyz, box, 2, index);
    done_atom(&atoms);

    std::unique_ptr<StructureFile> file = openFileAs<StructureFile>(path, FileMode::Read);
    file->read();
    ASSERT_EQ(file->getNumAtoms(), 2);
    EXPECT_REAL_EQ(0.5, file->getPositions()[0][XX]);
    EXPECT_REAL_EQ(4.5, file->getPositions()[1][ZZ]);
}

TEST(FileFormatRegistryTest, ConfioRefusesFormatWithoutParticleNames)
{
    TestFileManager             fileManager;
    const std::filesystem::path path = fileManager.getTemporaryFilePath("unnamed.gmxtestxyz");
    std::ofstream(path) << "0.5 1.0 1.5\n";
    t_symtab symtab;
    open_symtab(&symtab);
    char*   name = nullptr;
    t_atoms atoms;
    rvec*   x = nullptr;
    matrix  box;
    EXPECT_THROW_GMX(readConfAndAtoms(path, &symtab, &name, &atoms, nullptr, &x, nullptr, box),
                     InvalidInputError);
    done_symtab(&symtab);
}

} // namespace
} // namespace test
} // namespace gmx
