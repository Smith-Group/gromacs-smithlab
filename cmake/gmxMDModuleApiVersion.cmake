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


# Computes the version hash of the MDModules plugin interface.
#
# The hash covers the content of the interface headers, so any change to
# them changes it, including changes to comments. Line endings are
# normalized, so the hash is the same on all platforms.
#
# Usage:
#   gmx_mdmodule_api_hash(<output variable> <base directory> <header> ...)
# The headers are given relative to the base directory; their paths are
# part of the hash. The result is 16 hexadecimal digits.
function(gmx_mdmodule_api_hash _output _baseDirectory)
    set(_content "")
    foreach(_header ${ARGN})
        file(READ "${_baseDirectory}/${_header}" _headerContent)
        string(REPLACE "\r" "" _headerContent "${_headerContent}")
        string(APPEND _content "${_header}\n${_headerContent}")
    endforeach()
    string(SHA256 _hash "${_content}")
    string(SUBSTRING "${_hash}" 0 16 _hash)
    set(${_output} "${_hash}" PARENT_SCOPE)
endfunction()
