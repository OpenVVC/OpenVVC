include config.mak

# Virtual path points to src dir so the files search is done
# from there
vpath $(SRC_DIR):../
BUILD_DIR := build/

# Disable suffixes rules to ensure we are using our rules
.SUFFIXES:

# Defintion of program dependencies
PROG=dectest
PROG_LIBS=ovvcdec ovvcdmx ovvcutils
ovvcdec_LIBS = ovvcdmx ovvcutils
ovvcdmx_LIBS = ovvcutils

# Various glabal defintions to be moved in a config file
CFLAGS=--pedantic -Wall -fPIC -O0 -g
SHARED_LIBSUFF:=.so
STATIC_LIBSUFF:=.a

# list to keep track of all objects
ALL_OBJS :=
ALL_LIBOBJS :=

# List libs available for build
ALL_LIBS := ovvcutils
ALL_LIBS += ovvcdec
ALL_LIBS += ovvcdmx
# rpath-link so linker knows where to find libs
LD_FLAGS = -Wl,-rpath-link=$(BUILD_DIR):$(BUILD_DIR)libovvcdmx:$(BUILD_DIR)libovvcdec:$(BUILD_DIR)libovvcutils
# rpath so loader knows where to find  libs
PLD_FLAGS = -Wl,-rpath=$(BUILD_DIR):$(BUILD_DIR)libovvcdmx:$(BUILD_DIR)libovvcdec:$(BUILD_DIR)libovvcutils

# Create a list of openvvc libs flags used by test program
#TODO put these under a macro so we can create more than one
# program and permit to select either static or shared
PROG_LDFLAGS    :=$(foreach lib, $(PROG_LIBS), -l$(lib))
PROG_DEPLIBOBJS :=$(foreach lib, $(PROG_LIBS),$(BUILD_DIR)lib$(lib)/lib$(lib).so)
PROG_DEPLIBDIRS :=$(foreach lib, $(PROG_LIBS), -L$(BUILD_DIR)lib$(lib))

PROG_OBJS:=$(BUILD_DIR)$(PROG:%=%.o)

ALL_OBJS += $(PROG_OBJS)

all: $(BUILD_DIR)$(PROG) ${ALL_LIBOBJS}

# Macro called for each lib available for build
# It creates per libs specific list of objects
#
define lib_rules
$(ltarget_shared): $(OBJS) $(DEPLIBS)
	$(CC)  $(LD_FLAGS) --shared $^ -o $@ $(IFLAGS)

$(warning $(ltarget_shared): $(OBJS)  $(DEPLIBS))
$(warning	$(CC)  $(LD_FLAGS) --shared $^ -o $@ $(IFLAGS) )

# FIXME do we need dependencies besides OBJS in static lib
$(ltarget_static): $(OBJS)
	$(AR) rcD $@ $^
	ranlib $@
endef

define lib_objs_rules_template

# lsubdir evaluates to the sub directory containing the lib src
# $(1) is the root name of the lib
lname  := lib$(1)
lsubdir:= lib$(1)/

ltarget_shared :=$$(lsubdir)$$(lname)$$(SHARED_LIBSUFF)
ltarget_static :=$$(lsubdir)$$(lname)$$(STATIC_LIBSUFF)

include $(SRC_DIR)/$$(lsubdir)Makefile

# OBJS is a list of .o objects defined in each lib sub Makefiles
# we sort it and store it prefixed with corresponding subdir name
# in a OBJS variable prefixed with libname root pattern
OBJS      := $$(sort $$(OBJS:%=$$(BUILD_DIR)$$(lsubdir)%))
$(1)_OBJS := $$(OBJS)

# Append per subdir objects to general OBJS list
ALL_OBJS  += $$(OBJS)


# Rule to build dynamic library
#DEPLIBS := $$(foreach l,$$($(1)_LIBS),lib$$(l)/lib$$(l).so)
#IFLAGS  := $$(foreach l,$$($(1)_LIBS),-I./lib$$(l))

include $(SRC_DIR)/shared.mak


# Append lib targets to general lib objects list
ALL_LIBOBJS += $$(ltarget_shared)
ALL_LIBOBJS += $$(ltarget_static)
endef

# Create libs related variables and rules by calling lib macro
$(foreach lib,$(ALL_LIBS),$(eval $(call lib_objs_rules_template,$(lib))))

# Suffix rule on .c objects to build corresponding .o and .d objects
$(BUILD_DIR)%.o: %.c
	@mkdir -p $(@D)
	$(CC) -c $< -o $@ -MMD -MF $(@:.o=.d) -MT $@ $(CFLAGS) -I.

$(BUILD_DIR)${PROG}: ${PROG_DEPLIBOBJS} $(PROG_OBJS)
	$(CC) $(PLD_FLAGS) $(LD_FLAGS) -o $@  $(filter %.o, $^) $(PROG_LDFLAGS) $(PROG_DEPLIBDIRS)

.PHONY: all clean

# Force .o files to depend on the content of their associated .d file
# if it already exists which will ensure the .o is rebuild when one of
# its previous dependencies are modified
$(ALL_OBJS):
include $(wildcard $(ALL_OBJS:.o=.d))

#TODO Create a clean rule based on every objects and suffixes expansion
# to be sure not to forget anything
clean:
	$(RM) -r $(BUILD_DIR)
# $(RM) $(ALL_OBJS) $(ALL_OBJS:.o=.d) $(ALL_LIBOBJS) $(PROG)
# $(RM) -rf $(BUILD_DIR)**.o $(BUILD_DIR)**.d $(BUILD_DIR)**.a $(BUILD_DIR)**.so
