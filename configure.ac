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

AC_SUBST([LIB_VERSION], m4_esyscmd([printf $(cat VERSION | sed  's/\./:/g' )]))

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

#Check for SLHDR lib

PKG_HAVE_DEFINE_WITH_MODULES([SLHDR], [slhdr_cwrap], [ Use SL HDR ], [yes])

# --enable-werror
AC_ARG_ENABLE([werror], [AS_HELP_STRING([--enable-werror], [treat warnings as errors [no]])],
                [CFLAGS="-Werror $CFLAGS"], []
)

# --disable-simd
AC_ARG_ENABLE([simd], [AS_HELP_STRING([--disable-simd], [disable all simd optimisations [no]])],
                [
                AS_IF([test $enable_simd = "no"], [
                                AC_DEFINE([NO_SIMD],[1], [If define, it will remove SIMD compilation])
                        ]
                )
                ], 
                [enable_simd="yes"]
)


# --disable-arm-asm
AC_ARG_ENABLE([arm-asm], [AS_HELP_STRING([--disable-arm-asm], [disables arm assembly optimisation [no]])],
                [
                AS_IF([test $enable_arm_asm = "yes"], [
                                AC_DEFINE([ARCH_AARCH64_ASSEMBLY],[1], [Define the use of AARCH64 ASSEMBLY optimisation])
                        ]
                )
                ],
                [
                        enable_arm_asm="yes"
                        AC_DEFINE([ARCH_AARCH64_ASSEMBLY],[1], [Define the use of AARCH64 ASSEMBLY optimisation])
                ]
)

# --enable-arm-simde
AC_ARG_ENABLE([arm-simde], [AS_HELP_STRING([--enable-arm-simde], [enables arm simde optimisation [no]])],
                [
                AS_IF([test $enable_arm_simde = "yes"], [
                                AC_DEFINE([ARM_SIMDE],[1], [Define the use of SIMDE optimisation conversion])
                                AC_DEFINE([SIMDE_NO_CHECK_IMMEDIATE_CONSTANT],[1], [Define SIMDE_NO_CHECK_IMMEDIATE_CONSTANT])
                        ])
                ], 
                [
                        enable_arm_simde="yes" 
                        AC_DEFINE([ARM_SIMDE],[1], [Define the use of SIMDE optimisation conversion])
                        AC_DEFINE([SIMDE_NO_CHECK_IMMEDIATE_CONSTANT],[1], [Define SIMDE_NO_CHECK_IMMEDIATE_CONSTANT])
                ]
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
                AC_DEFINE([PREFIX],[1], [Define if assembly function should be prefixed])
                AC_DEFINE([HAVE_POSIX_MEMALIGN],[1], [Define the use of posix memalign malloc])
                AC_DEFINE([EXCLUDE_FOR_CLANG],[1], [Exclude files from clang compilation])
                isClang="true"
        ],
        [mingw*], [
                AC_DEFINE([HAVE_ALIGNED_MALLOC],[1], [Define the use of align malloc])
                
        ],
        [*android*], [
                AC_DEFINE([EXCLUDE_FOR_CLANG],[1], [Exclude files from clang compilation])
                isClang="true"
        ],
        [linux*|cygwin*|msys*], [
                AC_DEFINE([HAVE_POSIX_MEMALIGN],[1], [Define the use of posix memalign malloc])
        ], [
                AC_DEFINE([HAVE_POSIX_MEMALIGN],[1], [Define the use of posix memalign malloc])
        ]
)

# Check compiler build capability
AS_IF([test x"$X86" = x"true"], [
                AX_CHECK_COMPILE_FLAG([-msse4.1], [
                        flag_sse4_1="true" CFLAGS="-msse4.1 $CFLAGS"
                        AC_DEFINE([HAVE_SSE4_1],[1], [Define the use of SSE4_1 optimisation])
                ])
        ])
AM_CONDITIONAL([HAVE_SSE4_1], [test x"$enable_simd" = x"yes" -a x"$flag_sse4_1" = x"true"])

AS_IF([test x"$X86" = x"true"], [
                AX_CHECK_COMPILE_FLAG([-mavx2], [
                        flag_avx2="true" CFLAGS="-mavx2 $CFLAGS"
                        AC_DEFINE([HAVE_AVX2],[1], [Define the use of AVX2 optimisation])
                ])
        ])
AM_CONDITIONAL([HAVE_AVX2], [test x"$enable_simd" = x"yes" -a x"$flag_avx2" = x"true"])

AS_IF([test x"$flag_sse4_1" = x"true" -o x"$flag_avx2" = x"true"], [
                AC_DEFINE([HAVE_X86_OPTIM],[1], [Define the use of X86 optimisation])
        ])
AM_CONDITIONAL([HAVE_X86_OPTIM], [test x"$X86" = x"true" -a x"$enable_simd" = x"yes"])



# Define Makefile Conditionals

AM_CONDITIONAL([HAVE_NEON], [test x"$ARM" = x"true" -a x"$enable_simd" = x"yes"])

AM_CONDITIONAL([HAVE_NEON_ASM], [test x"$ARM" = x"true" -a x"$enable_simd" = x"yes" -a x"$enable_arm_asm" = x"yes"])

AM_CONDITIONAL([HAVE_NEON_SIMDE], [test x"$ARM" = x"true" -a x"$enable_simd" = x"yes" -a x"$enable_arm_simde" = x"yes"])

AM_CONDITIONAL([HAVE_EXCLUDE_CLANG], [test x"$ARM" = x"true" -a x"$isClang" = x"true"])

AC_OUTPUT
