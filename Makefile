include config.mak

# Compiler Verbosity Control
USER_CC := $(CC)
CC_0 = @echo "Compiling $<..."; $(USER_CC)
CC_1 = $(USER_CC)
CC = $(CC_$(VERBOSITY))

# Set default BUILD_TYPE
ifndef BUILD_TYPE
	BUILD_TYPE=release
endif

BUILDDIR_TYPE:=$(addprefix $(BUILDDIR)/, $(BUILD_TYPE)/)

# Set default ARCH
ifndef ARCH
ARCH=x86
endif

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

all: libs examples
	# all: $(BUILDDIR_TYPE)$(LIB_NAME)$(STATIC_LIBSUFF) $(BUILDDIR_TYPE)${PROG}
libs: $(BUILDDIR_TYPE)$(LIB_NAME)$(STATIC_LIBSUFF)

examples: $(BUILDDIR_TYPE)$(PROG)


$(BUILDDIR_TYPE)$(PROG):  $(BUILDDIR_TYPE)$(PROG).o $(BUILDDIR_TYPE)$(LIB_NAME)$(STATIC_LIBSUFF)
	$(CC) $^ -o $@


$(BUILDDIR_TYPE)$(LIB_NAME)$(STATIC_LIBSUFF): $(LIB_OBJ)
	$(AR) rcD $@ $^
	ranlib $@


$(BUILDDIR_TYPE)%.o: %.c
	@mkdir -p $(@D)
	$(CC) -c $< -o $@ -MMD -MF $(@:.o=.d) -MT $@ $(CFLAGS) -I$(SRC_FOLDER)

.PHONY: style check-style tidy
style:
	@for src in $(LIB_FILE) ; do \
		echo "Formatting $$src..." ; \
		clang-format -i "$$src" ; \
	done
	@echo "Done"


check-style:
	@for src in $(LIB_FILE) ; do \
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
	@echo "Style check passed"

.PHONY: clean

clean:
	rm -r $(BUILDDIR)
