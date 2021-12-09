#!/bin/bash
VERSION_FILE=$1
VERSION_HEADER=$2
BUILDDIR=$3

update_version(){
cat > ${BUILDDIR}${VERSION_HEADER} << EOF
#ifndef OVVERSION_H
#define OVVERSION_H

#define VER_MAJOR ${MAJOR}
#define VER_MINOR ${MINOR}
#define VER_REVISION ${REVISION}
#define VER_BUILD "${BUILD}"
#define OV_VERSION "${MAJOR}.${MINOR}.${REVISION}"

#define OV_STR(s) #s
#define OV_VERSION_STR(maj,min,rev,build) OV_STR(maj)OV_STR(.)OV_STR(min)OV_STR(.)OV_STR(rev)OV_STR(-)build

#endif // OVVERSION_H
EOF
}

# Parse file given in $1 argument
VERSION=($(sed 's/\./\ /g' ${VERSION_FILE}))

MAJOR=${VERSION[0]}
MINOR=${VERSION[1]}
REVISION=${VERSION[2]}

# Try to recover short SHA1 of current commit
BUILD=$(git rev-parse --short HEAD 2>/dev/null || echo "")

# Check if uncommited changes exists
GIT_INDEX=$(git update-index --refresh || echo "")

if [ ! -z "${GIT_INDEX}" ]; then
  BUILD="${BUILD}-dirty"
fi

PREVIOUS_VERSION=""
PREVIOUS_BUILD=""

if [ ! -f ${BUILDDIR}${VERSION_HEADER} ]; then
  echo "ovversion.h does not exist."
  echo "Creating ovversion.h..."
  update_version
else
  echo "ovversion.h already exists."
  PREVIOUS_VERSION=$((awk '/#define OV_VERSION "*"/ { print $3 }' ${BUILDDIR}${VERSION_HEADER} || echo "") | sed 's/"//g')
  PREVIOUS_BUILD=$((awk '/#define VER_BUILD "*"/ { print $3 }' ${BUILDDIR}${VERSION_HEADER} || echo "") | sed 's/"//g')
  if [ ${PREVIOUS_VERSION} != ${MAJOR}.${MINOR}.${REVISION} ] || [ ${PREVIOUS_BUILD} != ${BUILD} ]; then
    echo "Updating ovversion.h..."
    update_version
  else
    echo "Version already up to date. ovversion.h remains unchanged."
  fi
fi




# source ${BUILDDIR}config.sh
# prefix=${INSTALL_PREFIX}
# libdir=${prefix}/lib
# includedir=${prefix}/include
#
# cat > ${BUILDDIR}libopenvvc.pc << EOF
# prefix=${INSTALL_PREFIX}
# libdir=${prefix}/lib
# includedir=${prefix}/include
# extra_ld_flags=${EXTRA_LD_FLAGS}
# extra_cflags=${EXTRA_CFLAGS}
#
# Name: OpenVVC Library
# Description: Open Source VVC Decoder Library
# Version: ${MAJOR}.${MINOR}.${REVISION}
# Requires:
# Conflicts:
# Libs: -L\${libdir}/libopenvvc -Wl,-rpath=\${libdir}/libopenvvc -lovvc \${extra_ld_flags}
# Libs.private: -lpthreads
# Cflags: -I\${includedir}/libopenvvc -I\${libdir}/libopenvvc \${extra_cflags}
# EOF
