#!/bin/bash

# Parse file given in $1 argument
VERSION=($(sed 's/\./\n/g' $1))

MAJOR=${VERSION[0]}
MINOR=${VERSION[1]}
REVISON=${VERSION[2]}

# Try to recover short SHA1 of current commit
BUILD=$(git rev-parse --short HEAD 2>/dev/null || echo "")

echo "#ifndef OVVERSION_H
#define OVVERSION_H

#define VER_MAJOR ${MAJOR}
#define VER_MINOR ${MINOR}
#define VER_REVISION ${REVISON}
#define VER_BUILD \"${BUILD}\"

#endif // OVVERSION_H" > $2
