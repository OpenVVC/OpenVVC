OpenVVC
=======
An open-source VVC decoder library licensed under LGPLv2.1. OpenVVC is still under development.

|     | Windows   | Linux     | MacOS     | Android     |
|-----|-----------|-----------|-----------|-------------|
| x86 | Supported | Supported | Supported |  -          |
| ARM | -         | Supported | Supported |  Supported  |

A detailed list of supported tools is available on our [website](https://openvvc.github.io/#supported-tools).

## Table of Contents
- [Changelog](#changelog)
- [Using OpenVVC](#using-openvvc-test-program)
  - [Example](#example)
  - [Parameters](#parameters)
- [Compiling OpenVVC](#compiling-openvvc)
  - [Minimal Compilation Steps](#minimal-compilation-steps)
  - [For Windows](#for-windows)
  - [Configure Parameters](#configure-parameters)

## Changelog
### v1.0.0
* Support decoding of 8-bit bitstreams.
* Support Frame + Tile parallelism.
* Support RPR.
* Build system was changed to use Autotool.

## Using OpenVVC Test Program

### Example:

This library includes a simple executable to test the library. The executable can used as follows:

    dectest -i <file>.266 -o out.yuv

### Parameters:

```
usage: dectest [options]
options:
	-h, --help                                Show this message.
	-v, --version                             Show version information.
	-l <level>, --log-level=<level>           Define the level of verbosity. Value between 0 and 6. (Default: 2)
	-i <file>, --infile=<file>                Path to the file to be decoded (Default: test.266).
	-o <file>, --outfile=<file>               Path to the output file (Default: test.yuv).
	-t <nbthreads>, --framethr=<nbthreads>    Number of simultaneous frames decoded (Default: 0).
	-e <nbthreads>, --entrythr=<nbthreads>    Number of simultaneous entries decoded per frame (Default: 0).

```

## Compiling OpenVVC
### Prerequisites
```
sudo apt install build-essential curl autoconf automake libtool
```

### Minimal Compilation Steps:
```
autoreconf -if
./configure
make
```
To test the library, you can perform:
```
./CI/checkMD5.sh ./CI/test_bitstreams/random_access ./examples/dectest http://openvvc.insa-rennes.fr/bitstreams/JVET_CTC/random_access/
```

### For Windows
Build on windows was tested using Cygwin and MSYS2. After installing the dependencies, the compilation steps are identical to previous [section](#minimal-compilation-steps).

### Configure Parameters:
```
`configure' configures openVVC 1.0.0 to adapt to many kinds of systems.

Usage: ../../configure [OPTION]... [VAR=VALUE]...

To assign environment variables (e.g., CC, CFLAGS...), specify them as
VAR=VALUE.  See below for descriptions of some of the useful variables.

Defaults for the options are specified in brackets.

Configuration:
  -h, --help              display this help and exit
      --help=short        display options specific to this package
      --help=recursive    display the short help of all the included packages
  -V, --version           display version information and exit
  -q, --quiet, --silent   do not print `checking ...' messages
      --cache-file=FILE   cache test results in FILE [disabled]
  -C, --config-cache      alias for `--cache-file=config.cache'
  -n, --no-create         do not create output files
      --srcdir=DIR        find the sources in DIR [configure dir or `..']

Installation directories:
  --prefix=PREFIX         install architecture-independent files in PREFIX
                          [/usr/local]
  --exec-prefix=EPREFIX   install architecture-dependent files in EPREFIX
                          [PREFIX]

By default, `make install' will install all the files in
`/usr/local/bin', `/usr/local/lib' etc.  You can specify
an installation prefix other than `/usr/local' using `--prefix',
for instance `--prefix=$HOME'.

For better control, use the options below.

Fine tuning of the installation directories:
  --bindir=DIR            user executables [EPREFIX/bin]
  --sbindir=DIR           system admin executables [EPREFIX/sbin]
  --libexecdir=DIR        program executables [EPREFIX/libexec]
  --sysconfdir=DIR        read-only single-machine data [PREFIX/etc]
  --sharedstatedir=DIR    modifiable architecture-independent data [PREFIX/com]
  --localstatedir=DIR     modifiable single-machine data [PREFIX/var]
  --runstatedir=DIR       modifiable per-process data [LOCALSTATEDIR/run]
  --libdir=DIR            object code libraries [EPREFIX/lib]
  --includedir=DIR        C header files [PREFIX/include]
  --oldincludedir=DIR     C header files for non-gcc [/usr/include]
  --datarootdir=DIR       read-only arch.-independent data root [PREFIX/share]
  --datadir=DIR           read-only architecture-independent data [DATAROOTDIR]
  --infodir=DIR           info documentation [DATAROOTDIR/info]
  --localedir=DIR         locale-dependent data [DATAROOTDIR/locale]
  --mandir=DIR            man documentation [DATAROOTDIR/man]
  --docdir=DIR            documentation root [DATAROOTDIR/doc/openvvc]
  --htmldir=DIR           html documentation [DOCDIR]
  --dvidir=DIR            dvi documentation [DOCDIR]
  --pdfdir=DIR            pdf documentation [DOCDIR]
  --psdir=DIR             ps documentation [DOCDIR]

Program names:
  --program-prefix=PREFIX            prepend PREFIX to installed program names
  --program-suffix=SUFFIX            append SUFFIX to installed program names
  --program-transform-name=PROGRAM   run sed PROGRAM on installed program names

System types:
  --build=BUILD     configure for building on BUILD [guessed]
  --host=HOST       cross-compile to build programs to run on HOST [BUILD]

Optional Features:
  --disable-option-checking  ignore unrecognized --enable/--with options
  --disable-FEATURE       do not include FEATURE (same as --enable-FEATURE=no)
  --enable-FEATURE[=ARG]  include FEATURE [ARG=yes]
  --enable-silent-rules   less verbose build output (undo: "make V=1")
  --disable-silent-rules  verbose build output (undo: "make V=0")
  --enable-dependency-tracking
                          do not reject slow dependency extractors
  --disable-dependency-tracking
                          speeds up one-time build
  --enable-shared[=PKGS]  build shared libraries [default=yes]
  --enable-static[=PKGS]  build static libraries [default=yes]
  --enable-fast-install[=PKGS]
                          optimize for fast installation [default=yes]
  --disable-libtool-lock  avoid locking (might break parallel builds)
  --disable-slhdr         disable slhdr [no]
  --enable-werror         treat warnings as errors [no]
  --disable-simd          disable all simd optimisations [no]
  --disable-arm-asm       disables arm assembly optimisation [no]
  --enable-arm-simde      enables arm simde optimisation [no]

Optional Packages:
  --with-PACKAGE[=ARG]    use PACKAGE [ARG=yes]
  --without-PACKAGE       do not use PACKAGE (same as --with-PACKAGE=no)
  --with-pic[=PKGS]       try to use only PIC/non-PIC objects [default=use
                          both]
  --with-aix-soname=aix|svr4|both
                          shared library versioning (aka "SONAME") variant to
                          provide on AIX, [default=aix].
  --with-gnu-ld           assume the C compiler uses GNU ld [default=no]
  --with-sysroot[=DIR]    Search for dependent libraries within DIR (or the
                          compiler's sysroot if not specified).

Some influential environment variables:
  CC          C compiler command
  CFLAGS      C compiler flags
  LDFLAGS     linker flags, e.g. -L<lib dir> if you have libraries in a
              nonstandard directory <lib dir>
  LIBS        libraries to pass to the linker, e.g. -l<library>
  CPPFLAGS    (Objective) C/C++ preprocessor flags, e.g. -I<include dir> if
              you have headers in a nonstandard directory <include dir>
  CCAS        assembler compiler command (defaults to CC)
  CCASFLAGS   assembler compiler flags (defaults to CFLAGS)
  LT_SYS_LIBRARY_PATH
              User-defined run-time library search path.
  CPP         C preprocessor
  PKG_CONFIG  path to pkg-config utility
  PKG_CONFIG_PATH
              directories to add to pkg-config's search path
  PKG_CONFIG_LIBDIR
              path overriding pkg-config's built-in search path
  SLHDR_CFLAGS
              C compiler flags for SLHDR, overriding pkg-config
  SLHDR_LIBS  linker flags for SLHDR, overriding pkg-config

Use these variables to override the choices made by `configure' or to help
it to find libraries and programs with nonstandard names/locations.

Report bugs to the package provider.

```
