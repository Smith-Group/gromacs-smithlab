#
# This file is part of the GROMACS molecular simulation package.
#
# Copyright 2015- The GROMACS Authors
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

# Set up the module library
add_library(mdtypes INTERFACE)

set(MDTYPES_SOURCES
    awh_history.cpp
    checkpointdata.cpp
    commrec.cpp
    df_history.cpp
    energyhistory.cpp
    forcebuffers.cpp
    group.cpp
    iforceprovider.cpp
    inputrec.cpp
    interaction_const.cpp
    md_enums.cpp
    mdatoms.cpp
    mdmoduleplugin.cpp
    multipletimestepping.cpp
    observableshistory.cpp
    observablesreducer.cpp
    output_control.cpp
    state.cpp
    threaded_force_buffer.cpp)

if(GMX_GPU)
    gmx_add_libgromacs_sources(
       state_propagator_data_gpu_impl_gpu.cpp
       )
   if(GMX_GPU_CUDA)
       _gmx_add_files_to_property(CUDA_SOURCES
           state_propagator_data_gpu_impl_gpu.cpp
           )
   endif()
   if(GMX_GPU_SYCL)
       _gmx_add_files_to_property(SYCL_SOURCES
           state_propagator_data_gpu_impl_gpu.cpp
           )
   endif()
   if(GMX_GPU_HIP)
       _gmx_add_files_to_property(HIP_SOURCES
           state_propagator_data_gpu_impl_gpu.cpp
           )
   endif()
else()
    gmx_add_libgromacs_sources(
      state_propagator_data_gpu_impl.cpp
      )
endif()

# Source files have the following private module dependencies.
target_link_libraries(mdtypes PRIVATE
                      #                      gmxlib
                      #                      math
                      #                      mdtypes
                      #                      tng_io
                      )

# Public interface for modules, including dependencies and interfaces
#target_include_directories(mdtypes PUBLIC
#                           $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>)
#target_link_libraries(mdtypes PUBLIC
target_link_libraries(mdtypes INTERFACE
                      legacy_api
                      gpu_utils
                      math
                      )
# Interface conveys transitive dependence on other interface(s)
target_link_libraries(mdtypes INTERFACE utility)

# TODO: when fileio is an OBJECT target
#target_link_libraries(mdtypes PUBLIC legacy_api)
#target_link_libraries(mdtypes PRIVATE common)

# Source files have the following private module dependencies.
#target_link_libraries(mdtypes PRIVATE tng_io)

# Headers that are installed as part of the legacy API. The module owns its
# part of the installed include tree and installs it through a target file set.
target_include_directories(mdtypes INTERFACE
                           $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>)

target_sources(mdtypes
    INTERFACE
        FILE_SET public_headers
        TYPE HEADERS
        BASE_DIRS ${CMAKE_CURRENT_SOURCE_DIR}/include
        FILES
            include/gromacs/mdtypes/inputrec.h
            include/gromacs/mdtypes/md_enums.h
)

if(GMX_INSTALL_LEGACY_API)
    install(TARGETS mdtypes FILE_SET public_headers DESTINATION include)
endif()

# The version of the MDModules plugin interface
include(${CMAKE_CURRENT_SOURCE_DIR}/mdmoduleapi.cmake)
include(gmxMDModuleApiVersion)
gmx_mdmodule_api_hash(GMX_MDMODULE_API_HASH ${PROJECT_SOURCE_DIR}/src ${GMX_MDMODULE_API_HEADERS})
# Reconfigure when the interface headers change, so the hash stays current
set(_mdmoduleApiHeaderPaths ${GMX_MDMODULE_API_HEADERS})
list(TRANSFORM _mdmoduleApiHeaderPaths PREPEND "${PROJECT_SOURCE_DIR}/src/")
set_property(DIRECTORY APPEND PROPERTY CMAKE_CONFIGURE_DEPENDS
             ${_mdmoduleApiHeaderPaths} ${CMAKE_CURRENT_SOURCE_DIR}/mdmoduleapi.cmake)
if(GMX_LIB_MPI)
    set(GMX_MDMODULE_API_MPI 2)
elseif(GMX_THREAD_MPI)
    set(GMX_MDMODULE_API_MPI 1)
else()
    set(GMX_MDMODULE_API_MPI 0)
endif()
configure_file(mdmoduleapiversion.h.cmakein
               ${PROJECT_BINARY_DIR}/src/include/gromacs/mdtypes/mdmoduleapiversion.h)

# Headers that MDModules plugins need beyond the legacy API, installed
# when plugins can be loaded
if(GMX_MDMODULE_PLUGINS)
    target_sources(mdtypes
        INTERFACE
            FILE_SET mdmodule_plugin_headers
            TYPE HEADERS
            BASE_DIRS ${PROJECT_SOURCE_DIR}/src
            FILES
                forceoutput.h
                iforceprovider.h
                imdmodule.h
                imdoutputprovider.h
                imdpoptionprovider.h
                mdmoduleplugin.h
                mdmodulepluginentrypoints.h
    )
    if(GMX_INSTALL_LEGACY_API)
        install(TARGETS mdtypes FILE_SET mdmodule_plugin_headers DESTINATION include)
        install(FILES ${PROJECT_BINARY_DIR}/src/include/gromacs/mdtypes/mdmoduleapiversion.h
                DESTINATION include/gromacs/mdtypes)
    endif()
endif()

if (BUILD_TESTING)
    add_subdirectory(tests)
endif()

list(TRANSFORM MDTYPES_SOURCES PREPEND "${CMAKE_CURRENT_SOURCE_DIR}/")
set(LIBGROMACS_SOURCES ${LIBGROMACS_SOURCES} ${MDTYPES_SOURCES} PARENT_SCOPE)
