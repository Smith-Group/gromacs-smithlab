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


# The version of the interface that MDModules plugins are built against.
# See "Plugin API version" in docs/doxygen/lib/mdmodules.md.
#
# A plugin built for one major runs on every host with the same major.
# Within a major, the interface only grows: fields are appended to the
# structs that the host fills, and new hooks are new interfaces. Any other
# change increments the major.
set(GMX_MDMODULE_API_MAJOR 1)

# The headers of the interface, relative to src/.
set(GMX_MDMODULE_API_HEADERS
    gromacs/mdrunutility/mdmodulesnotifier.h
    gromacs/mdrunutility/mdmodulesnotifiers.h
    gromacs/mdtypes/iforceprovider.h
    gromacs/mdtypes/imdmodule.h
    gromacs/mdtypes/imdoutputprovider.h
    gromacs/mdtypes/imdpoptionprovider.h
    gromacs/mdtypes/mdmoduleplugin.h)

# The hash of these headers that the major was last checked for. A test
# fails when the headers change without this hash being updated.
set(GMX_MDMODULE_API_RECORDED_HASH "d42d23d0ab7bcb09")
