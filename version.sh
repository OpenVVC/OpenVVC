#!/bin/bash

# Parse file given in $1 argument
VERSION=($(sed 's/\./\ /g' $1))

MAJOR=${VERSION[0]}
MINOR=${VERSION[1]}
REVISION=${VERSION[2]}

# Try to recover short SHA1 of current commit
BUILD=$(git rev-parse --short HEAD 2>/dev/null || echo "")

cat > $2 << EOF
#ifndef OVVERSION_H
#define OVVERSION_H

#define VER_MAJOR ${MAJOR}
#define VER_MINOR ${MINOR}
#define VER_REVISION ${REVISION}
#define VER_BUILD "${BUILD}"

#endif // OVVERSION_H
EOF

source config.sh
prefix=${INSTALL_PREFIX}
libdir=${prefix}/lib
includedir=${prefix}/include

cat > libopenvvc.pc << EOF
prefix=${INSTALL_PREFIX}
libdir=${prefix}/lib
includedir=${prefix}/include

Name: OpenVVC Library
Description: Open Source VVC Decoder Library
Version: ${MAJOR}.${MINOR}.${REVISION}
Requires:
Conflicts:
Libs: -L\${libdir}/libopenvvc -lovvc
Libs.private: -lpthreads
Cflags: -I\${includedir}/libopenvvc -I\${libdir}/libopenvvc
EOF
