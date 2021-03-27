#!/bin/bash

# Parse file given in $1 argument
VERSION=($(sed 's/\./\n/g' $1))

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
