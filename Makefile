include config.mak

# Set defaults
CC?=gcc
VERBOSITY?=0
CFLAGS_COMMON?= --pedantic -fPIC
CFLAGS_DEBUG?=  -O0 -g
CFLAGS_RELEASE?= -O3 -Werror
SHARED_LIBSUFF:=.so
STATIC_LIBSUFF:=.a
LD_FLAGS=-lpthread
BUILDDIR?=.
BUILD_TYPE?=RELEASE
ARCH?=x86


# Compiler Verbosity Control
USER_CC := $(CC)
CC_0 = @echo "$(USER_CC) $@"; $(USER_CC)
CC_1 = $(USER_CC)
CC = $(CC_$(VERBOSITY))

USER_AR := $(AR)
AR_0 = @echo "$(USER_AR) $@"; $(USER_AR)
AR_1 = $(USER_AR)
AR = $(AR_$(VERBOSITY))

AT_0 = @
AT_1 =
AT = $(AT_$(VERBOSITY))

BUILDDIR_TYPE:=$(addprefix $(BUILDDIR)/, $(shell echo $(BUILD_TYPE) | tr A-Z a-z)/)
# Handle flags depending of BUILD_TYPE
CFLAGS=$(CFLAGS_COMMON)
CFLAGS+=$(CFLAGS_$(BUILD_TYPE))

# Find Sources
include libovvc.mak
LIB_SRC:=$(addprefix $(SRC_FOLDER),$(LIB_SRC))
LIB_HEADER:=$(addprefix $(SRC_FOLDER),$(LIB_HEADER))
LIB_OBJ:=$(addprefix $(BUILDDIR_TYPE),$(LIB_SRC:%.c=%.o))
LIB_FILE:=$(LIB_HEADER) $(LIB_SRC)

include $(ARCH)libovvc.mak
$(ARCH)_LIB_SRC:=$(addprefix $($(ARCH)_SRC_FOLDER),$($(ARCH)_LIB_SRC))
$(ARCH)_LIB_OBJ:=$(addprefix $(BUILDDIR_TYPE),$($(ARCH)_LIB_SRC:%.c=%.o))
BUILDDIR_TYPE_ARCH:=$(addprefix $(BUILDDIR_TYPE), $($(ARCH)_SRC_FOLDER))

LIB_NAME:= libovvc

SHARED_LIBSUFF?=.so
STATIC_LIBSUFF?=.a
DEFAULT_LIBSUFF?=$(SHARED_LIBSUFF)


PROG=examples/dectest

ALL_OBJS=$(LIB_OBJ) $(addprefix $(BUILDDIR_TYPE),$(addsuffix .o, $(PROG))) $($(ARCH)_LIB_OBJ)

all: version libs examples

version:
	$(AT)./version.sh RELEASE $(SRC_FOLDER)$(LIB_VERSION_HEADER)

libs: version $(BUILDDIR_TYPE)$(LIB_NAME)$(STATIC_LIBSUFF) $(BUILDDIR_TYPE)$(LIB_NAME)$(SHARED_LIBSUFF)

examples: version $(BUILDDIR_TYPE)$(PROG)

$(BUILDDIR_TYPE)$(PROG):  $(BUILDDIR_TYPE)$(PROG).o $(BUILDDIR_TYPE)$(LIB_NAME)$(DEFAULT_LIBSUFF)
	$(CC) $^ -o $@


$(BUILDDIR_TYPE)$(LIB_NAME)$(STATIC_LIBSUFF): $(LIB_OBJ) $($(ARCH)_LIB_OBJ)
	$(AR) rcD $@ $^
	ranlib $@

$(BUILDDIR_TYPE)$(LIB_NAME)$(SHARED_LIBSUFF): $(LIB_OBJ) $($(ARCH)_LIB_OBJ)
	$(CC) -shared $^ -o $@ $(LD_FLAGS)

$(BUILDDIR_TYPE_ARCH)%_sse.o: $($(ARCH)_SRC_FOLDER)%_sse.c
	echo $(BUILDDIR_TYPE_ARCH)
	$(AT)mkdir -p $(@D)
	$(CC) -c $< -o $@ -MMD -MF $(@:.o=.d) -MT $@ $(CFLAGS) $(SSE_CFLAGS) -I$(SRC_FOLDER)

$(BUILDDIR_TYPE)%.o: %.c
	$(AT)mkdir -p $(@D)
	$(CC) -c $< -o $@ -MMD -MF $(@:.o=.d) -MT $@ $(CFLAGS) -I$(SRC_FOLDER)


.PHONY: style check-style tidy version
FILE_TO_STYLE:=$(shell find . -type f -name "*.[ch]")
style:
	$(AT)for src in $(FILE_TO_STYLE) ; do \
		echo "Formatting $$src..." ; \
		clang-format -i "$$src" ; \
	done
	$(AT)echo "Done"


check-style:
	$(AT)for src in $(FILE_TO_STYLE) ; do \
		var=`clang-format "$$src" | diff "$$src" - | wc -l` ; \
		# clang-tidy -checks='-*,readability-identifier-naming'\
		#     -config="{CheckOptions: [ \
		#     { key: readability-identifier-naming.NamespaceCase, value: lower_case },\
		#     { key: readability-identifier-naming.ClassCase, value: CamelCase  },\
		#     { key: readability-identifier-naming.StructCase, value: CamelCase  },\
		#     { key: readability-identifier-naming.FunctionCase, value: lower_case },\
		#     { key: readability-identifier-naming.VariableCase, value: lower_case },\
		#     { key: readability-identifier-naming.GlobalConstantCase, value: UPPER_CASE }\
		#     ]}" "$$src" ; \
		if [ $$var -ne 0 ] ; then \
			echo "$$src does not respect the coding style (diff: $$var lines)" ; \
			exit 1 ; \
		fi ; \
	done
	$(AT)echo "Style check passed"

.PHONY: clean mrproper

# Force .o files to depend on the content of their associated .d file
# if it already exists which will ensure the .o is rebuild when one of
# its previous dependencies are modified
$(ALL_OBJS):
include $(wildcard $(ALL_OBJS:.o=.d))

clean:
	$(AT)rm -f $(SRC_FOLDER)$(LIB_VERSION_HEADER)
	$(AT)rm -f $(ALL_OBJS) $(ALL_OBJS:.o=.d) $(addprefix $(BUILDDIR_TYPE),$(PROG)) $(BUILDDIR_TYPE)$(LIB_NAME)$(STATIC_LIBSUFF)

mrproper:
	$(AT)rm -f $(SRC_FOLDER)$(LIB_VERSION_HEADER)
	$(AT)rm -rf $(BUILDDIR)
