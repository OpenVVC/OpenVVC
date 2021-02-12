include config.mak

# Set defaults
CC?=gcc
VERBOSITY?=0
CFLAGS_COMMON?= --pedantic -fPIC
CFLAGS_DEBUG?=  -O0 -g
CFLAGS_RELEASE?= -O3 -Werror
SHARED_LIBSUFF:=.so
STATIC_LIBSUFF:=.a
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
LIB_NAME:= libovvc

SHARED_LIBSUFF:=.so
STATIC_LIBSUFF:=.a

PROG=examples/dectest

ALL_OBJS=$(LIB_OBJ) $(addprefix $(BUILDDIR_TYPE),$(addsuffix .o, $(PROG)))


all: version libs examples

version: $(SRC_FOLDER)$(LIB_VERSION_HEADER) RELEASE

libs: version $(BUILDDIR_TYPE)$(LIB_NAME)$(STATIC_LIBSUFF)

examples: version $(BUILDDIR_TYPE)$(PROG)

$(SRC_FOLDER)$(LIB_VERSION_HEADER): RELEASE
	$(AT)./version.sh $^ $@

$(BUILDDIR_TYPE)$(PROG):  $(BUILDDIR_TYPE)$(PROG).o $(BUILDDIR_TYPE)$(LIB_NAME)$(STATIC_LIBSUFF)
	$(CC) $^ -o $@


$(BUILDDIR_TYPE)$(LIB_NAME)$(STATIC_LIBSUFF): $(LIB_OBJ)
	$(AR) rcD $@ $^
	ranlib $@


$(BUILDDIR_TYPE)%.o: %.c
	$(AT)mkdir -p $(@D)
	$(CC) -c $< -o $@ -MMD -MF $(@:.o=.d) -MT $@ $(CFLAGS) -I$(SRC_FOLDER)

.PHONY: style check-style tidy version
style:
	$(AT)for src in $(LIB_FILE) ; do \
		echo "Formatting $$src..." ; \
		clang-format -i "$$src" ; \
	done
	$(AT)echo "Done"


check-style:
	$(AT)for src in $(LIB_FILE) ; do \
		var=`clang-format "$$src" | diff "$$src" - | wc -l` ; \
		clang-tidy -checks='-*,readability-identifier-naming' \
		    -config="{CheckOptions: [ \
		    { key: readability-identifier-naming.NamespaceCase, value: lower_case },\
		    { key: readability-identifier-naming.ClassCase, value: CamelCase  },\
		    { key: readability-identifier-naming.StructCase, value: CamelCase  },\
		    { key: readability-identifier-naming.FunctionCase, value: lower_case },\
		    { key: readability-identifier-naming.VariableCase, value: lower_case },\
		    { key: readability-identifier-naming.GlobalConstantCase, value: UPPER_CASE }\
		    ]}" "$$src" ; \
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
