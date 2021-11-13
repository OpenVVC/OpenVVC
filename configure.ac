AC_INIT([openVVC], m4_esyscmd([printf $(cat VERSION)]))

# Define directories
AC_CONFIG_AUX_DIR([build-aux])
AC_CONFIG_MACRO_DIR([m4])

AM_INIT_AUTOMAKE([foreign -Wall])
AC_PROG_CC
AM_PROG_AS
AM_PROG_AR

# Enable Libtool
LT_INIT

AC_CONFIG_HEADERS([ovconfig.h])
AC_CONFIG_FILES([Makefile libovvc/Makefile libovvc/x86/Makefile libovvc/arm/Makefile examples/Makefile libovvc/ovvc.pc])

# We need to force AX_PTHREAD to check -pthread -lpthread since otherwise
# it only outputs -pthread for GCC. Without -lpthread GCC does not link the
# shared library against the pthread library (even though it does link the
# executable).
PTHREAD_CFLAGS=-pthread
PTHREAD_LIBS=-lpthread

# This does workarounds for pthreads on various compilers.
AX_PTHREAD([],[AC_MSG_ERROR([POSIX threads not found])])

CFLAGS="$PTHREAD_CFLAGS $CFLAGS"
LIBS="$PTHREAD_LIBS $LIBS"
CC="$PTHREAD_CC"

# --enable-werror
AC_ARG_ENABLE([werror], [AS_HELP_STRING([--enable-werror], [treat warnings as errors [no]])],
              [CFLAGS="-Werror $CFLAGS"], []
)

# --disable-simd
AC_ARG_ENABLE([simd], [AS_HELP_STRING([--disable-simd], [disable all simd optimisations [no]])],
              [flag_disable_simd="true" CFLAGS="-DNO_SIMD $CFLAGS"], [flag_disable_simd="false"]
)


# --enable-arm-asm
AC_ARG_ENABLE([arm-asm], [AS_HELP_STRING([--enable-arm-asm], [enables arm assembly optimisation [no]])],
              [CFLAGS="-DARCH_AARCH64_ASSEMBLY=1 $CFLAGS"], []
)

# host and cpu specific settings
AS_CASE([$host_cpu],
        [i?86], [BITS="32" X86="true"],
        [x86_64], [BITS="64" X86="true"],
        [aarch32|arm32], [BITS="32" ARM="true"],
        [aarch64|arm64], [BITS="64" ARM="true"]
)

AS_CASE([$host_os],
        [darwin*], [
                CFLAGS="$CFLAGS -DHAVE_POSIX_MEMALIGN -DPREFIX=1"
        ],
        [cygwin*|msys*|mingw*], [
                CFLAGS="$CFLAGS -DHAVE_ALIGNED_MALLOC"
        ],
        [*android*], [
                CFLAGS="$CFLAGS -DHAVE_POSIX_MEMALIGN"
        ],
        [linux*], [
                CFLAGS="$CFLAGS -DHAVE_POSIX_MEMALIGN"
        ], [
                CFLAGS="$CFLAGS -DHAVE_POSIX_MEMALIGN"
        ]
)

# Check compiler build capability
AS_IF([test x"$X86" = x"true"], [AX_CHECK_COMPILE_FLAG([-msse4.1], [flag_sse4_1="true" CFLAGS="-msse4.1 $CFLAGS"])],[])

# Define Makefile Conditionals
AM_CONDITIONAL([HAVE_SSE4_1], [test x"$flag_sse4_1" = x"true" -a x"$X86" = x"true" -a x"$flag_disable_simd" = x"false"])

AM_CONDITIONAL([HAVE_NEON], [test x"$ARM" = x"true" -a x"$flag_disable_simd" = x"false"])

AC_OUTPUT
