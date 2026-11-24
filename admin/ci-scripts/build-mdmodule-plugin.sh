#!/usr/bin/env bash

# Build the electric-field test plugin against an installed GROMACS, which
# shows that an MDModules plugin needs only the installed headers, and run
# a few steps of water in which the installed mdrun loads it.
# Input variables:
#   CMAKE: Executable (path) for the CMake binary. (E.g. `$(which cmake)`)
#   GROMACS_DIR: absolute path to a GROMACS installation built with
#       GMX_MDMODULE_PLUGINS=ON and GMX_INSTALL_LEGACY_API=ON, and without
#       an MPI library
#   GMX_SUFFIX: (optional) suffix of that installation
#   TMP_BUILD: (optional) the name of a directory to create and use for the build

set -e

SOURCE_DIR=$PWD
CLIENT_CMAKE_ARGS="-C ${GROMACS_DIR}/share/cmake/gromacs${GMX_SUFFIX}/gromacs-hints${GMX_SUFFIX}.cmake \
 -Werror=dev \
 -DCMAKE_COMPILE_WARNING_AS_ERROR=ON \
 -DGROMACS_DIR=${GROMACS_DIR}"

if [ "${GMX_SUFFIX}" ] ; then
  CLIENT_CMAKE_ARGS="${CLIENT_CMAKE_ARGS} -DGMX_SUFFIX=${GMX_SUFFIX}"
fi

BUILD=${TMP_BUILD:-"mdmodule-plugin-build"}
mkdir -p ${BUILD}
pushd ${BUILD}
  ${CMAKE} ${SOURCE_DIR}/src/programs/mdrun/tests/plugins ${CLIENT_CMAKE_ARGS}
  ${CMAKE} --build . --target all -- VERBOSE=1

  mkdir -p run
  pushd run
    cp ${SOURCE_DIR}/src/testutils/simulationdatabase/spc216.gro \
       ${SOURCE_DIR}/src/testutils/simulationdatabase/spc216.top .
    cat > plugin.mdp <<MDP
integrator              = md
nsteps                  = 4
cutoff-scheme           = Verlet
coulombtype             = PME
rcoulomb                = 0.7
rvdw                    = 0.7
mdmodule-plugins        = electric-field-plugin
electric-field-plugin-z = 1 0 0 0
MDP
    GMX=${GROMACS_DIR}/bin/gmx${GMX_SUFFIX}
    # grompp and mdrun both load the plugin, which is in the build directory
    export GMX_MDMODULE_PLUGIN_PATH=$(dirname $PWD)
    ${GMX} grompp -f plugin.mdp -c spc216.gro -p spc216.top -o plugin.tpr
    ${GMX} mdrun -s plugin.tpr -ntmpi 1 -ntomp 1 -deffnm plugin
    grep "Loaded MDModules plugin electric-field-plugin" plugin.log
  popd
popd
