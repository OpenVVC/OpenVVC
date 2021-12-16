AC_INIT([openVVC], m4_esyscmd([printf $(cat VERSION)]))

# Define directories
AC_CONFIG_AUX_DIR([build-aux])
AC_CONFIG_MACRO_DIR([m4])

AM_INIT_AUTOMAKE([foreign -Wall])
AM_SILENT_RULES([yes])

AC_PROG_CC(gcc, clang)
AM_PROG_AS(gcc, clang)
AM_PROG_AR

# Enable Libtool
LT_INIT

AC_CONFIG_HEADERS([ovconfig.h])
AC_CONFIG_FILES([Makefile libovvc/Makefile libovvc/x86/Makefile libovvc/arm/Makefile examples/Makefile libovvc/libopenvvc.pc])

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

# --disable-slhdr
AC_ARG_ENABLE([slhdr], [AS_HELP_STRING([--disable-slhdr], [disable slhdr [no]])], [], [enable_slhdr=yes])

#Check for SLHDR lib
if test "x$enable_slhdr" = "xyes"
then
PKG_CHECK_MODULES([SLHDR], [slhdr_cwrap],[
        CPPFLAGS="$CPPFLAGS -DENABLE_SLHDR"
        CFLAGS="$SLHDR_CFLAGS $CFLAGS"
        LIBS="$SLHDR_LIBS $LIBS"
],
[
        AC_MSG_WARN([SLHDR library not found. SLHDR will be deactivated.])
])
fi

# --enable-werror
AC_ARG_ENABLE([werror], [AS_HELP_STRING([--enable-werror], [treat warnings as errors [no]])],
              [CFLAGS="-Werror $CFLAGS"], []
)

# --disable-simd
AC_ARG_ENABLE([simd], [AS_HELP_STRING([--disable-simd], [disable all simd optimisations [no]])],
              [AS_IF([test $enable_simd = "yes"], [CPPFLAGS="$CPPFLAGS -DNO_SIMD"])], [enable_simd="yes" CPPFLAGS="$CPPFLAGS -DNO_SIMD"]
)


# --disable-arm-asm
AC_ARG_ENABLE([arm-asm], [AS_HELP_STRING([--disable-arm-asm], [disables arm assembly optimisation [no]])],
              [AS_IF([test $enable_arm_asm = "yes"], [CPPFLAGS="$CPPFLAGS -DARCH_AARCH64_ASSEMBLY=1"])], [enable_arm_asm="yes" CPPFLAGS="$CPPFLAGS -DARCH_AARCH64_ASSEMBLY=1"]
)

# --enable-arm-simde
AC_ARG_ENABLE([arm-simde], [AS_HELP_STRING([--enable-arm-simde], [enables arm simde optimisation [no]])],
              [AS_IF([test $enable_arm_simde = "yes"], [CPPFLAGS="$CPPFLAGS -DARM_SIMDE=1 -DSIMDE_NO_CHECK_IMMEDIATE_CONSTANT"])], [$enable_arm_simde = "yes" CPPFLAGS="$CPPFLAGS -DARM_SIMDE=1 -DSIMDE_NO_CHECK_IMMEDIATE_CONSTANT"]
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
                CPPFLAGS="$CPPFLAGS -DPREFIX=1 -DHAVE_POSIX_MEMALIGN"
        ],
        [cygwin*|msys*|mingw*], [
                CPPFLAGS="$CPPFLAGS -DHAVE_ALIGNED_MALLOC"
        ],
        [*android*], [
                CPPFLAGS="$CPPFLAGS -DHAVE_POSIX_MEMALIGN -DEXCLUDE_FOR_ANDROID=1"
                isAndroid="true"
        ],
        [linux*], [
                CPPFLAGS="$CPPFLAGS -DHAVE_POSIX_MEMALIGN"
        ], [
                CPPFLAGS="$CPPFLAGS -DHAVE_POSIX_MEMALIGN"
        ]
)

# Check compiler build capability
AS_IF([test x"$X86" = x"true"], [AX_CHECK_COMPILE_FLAG([-msse4.1], [flag_sse4_1="true" CFLAGS="-msse4.1 $CFLAGS"])],[])

# Define Makefile Conditionals
AM_CONDITIONAL([HAVE_SSE4_1], [test x"$flag_sse4_1" = x"true" -a x"$X86" = x"true" -a x"$enable_simd" = x"yes"])

AM_CONDITIONAL([HAVE_NEON], [test x"$ARM" = x"true" -a x"$enable_simd" = x"yes"])

AM_CONDITIONAL([HAVE_NEON_ASM], [test x"$ARM" = x"true" -a x"$enable_simd" = x"yes" -a x"$enable_arm_asm" = x"yes"])

AM_CONDITIONAL([HAVE_NEON_SIMDE], [test x"$ARM" = x"true" -a x"$enable_simd" = x"yes" -a x"$enable_arm_simde" = x"yes"])

AM_CONDITIONAL([HAVE_EXCLUDE_ANDROID], [test x"$ARM" = x"true" -a x"$isAndroid" = x"true"])
AC_OUTPUT
