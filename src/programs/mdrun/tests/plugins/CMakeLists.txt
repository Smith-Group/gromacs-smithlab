#
# This file is part of the GROMACS molecular simulation package.
#
# Copyright 2026- The GROMACS Authors
# and the project initiators Erik Lindahl, Berk Hess and David van der Spoel.
# Consult the AUTHORS/COPYING files and https://www.gromacs.org for details.
#
# GROMACS is free software; you can redistribute it and/or
# modify it under the terms of the GNU Lesser General Public License
# as published by the Free Software Foundation; either version 2.1
# of the License, or (at your option) any later version.
#
# GROMACS is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
# Lesser General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public
# License along with GROMACS; if not, see
# https://www.gnu.org/licenses, or write to the Free Software Foundation,
# Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA.
#
# If you want to redistribute modifications to GROMACS, please
# consider that scientific software is very special. Version
# control is crucial - bugs must be traceable. We will be happy to
# consider code for inclusion in the official distribution, but
# derived work must not be called official GROMACS. Details are found
# in the README & COPYING files - if they are missing, get the
# official version at https://www.gromacs.org.
#
# To help us fund GROMACS development, we humbly ask that you cite
# the research papers on the package. Check out https://www.gromacs.org.


# Builds the electric-field test plugin against an installed GROMACS, to
# show that a plugin needs only the installed headers. This is a project
# of its own, not part of the GROMACS build:
#
#   cmake -S src/programs/mdrun/tests/plugins -B <build> -DGROMACS_DIR=<prefix>
#
# GROMACS has to be configured with GMX_MDMODULE_PLUGINS=ON and
# GMX_INSTALL_LEGACY_API=ON. Name the plugin directory in
# GMX_MDMODULE_PLUGIN_PATH to load the plugin.

cmake_minimum_required(VERSION 3.28)
project(electric-field-plugin CXX)

set(GROMACS_DIR "" CACHE PATH "The GROMACS installation to build against")
set(GMX_SUFFIX "" CACHE STRING "The suffix of that installation, e.g. _mpi")
find_package(gromacs${GMX_SUFFIX} REQUIRED CONFIG HINTS ${GROMACS_DIR}/share/cmake/gromacs${GMX_SUFFIX})
gromacs_check_compiler(CXX)

# Of the GROMACS source tree, the plugin uses only its own three files. The
# header of the module is copied, so that no other header of the tree is
# found. gmxpre.h, which every GROMACS source includes first, is not
# installed, and the plugin needs nothing from it.
set(_moduleSourceDirectory ${CMAKE_CURRENT_SOURCE_DIR}/../../../../gromacs/applied_forces)
set(_includeDirectory ${CMAKE_CURRENT_BINARY_DIR}/include)
configure_file(${_moduleSourceDirectory}/electricfield.h
               ${_includeDirectory}/gromacs/applied_forces/electricfield.h COPYONLY)
file(WRITE ${_includeDirectory}/gmxpre.h "")

add_library(electric-field-plugin MODULE
            electricfieldplugin.cpp
            ${_moduleSourceDirectory}/electricfield.cpp)
set_target_properties(electric-field-plugin PROPERTIES PREFIX "")
target_include_directories(electric-field-plugin PRIVATE ${_includeDirectory})
target_link_libraries(electric-field-plugin PRIVATE Gromacs::libgromacs)
if(NOT APPLE)
    # Call the plugin's own copy of the module, not the one in libgromacs
    target_link_options(electric-field-plugin PRIVATE LINKER:-Bsymbolic-functions)
endif()
