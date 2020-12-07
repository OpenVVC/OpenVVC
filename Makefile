include config.mak

vpath $(SRC_DIR)

.SUFFIXES:
# Defintion of program dependencies
PROG=dectest
PROG_LIBS=ovvcdec ovvcdmx ovvcutils

# Various glabal defintions to be moved in a config file
CFLAGS=-ansi --pedantic -Wall -fPIC -O0 -g
SHARED_LIBSUFF:=.so
STATIC_LIBSUFF:=.a

# list to keep track of all objects
ALL_OBJS :=
ALL_LIBOBJS :=

# List libs available for build
ALL_LIBS := ovvcutils
ALL_LIBS += ovvcdec
ALL_LIBS += ovvcdmx
LD_FLAGS = -Wl,-rpath-link=:libovvcdmx:libovvcdec:libovvcutils
PLD_FLAGS = -Wl,-rpath=:libovvcdmx:libovvcdec:libovvcutils

# Create a list of openvvc libs flags used by test program
#TODO put these under a macro so we can create more than one 
# program and permit to select either static or shared
PROG_LDFLAGS    :=$(foreach lib, $(PROG_LIBS), -l$(lib))
PROG_DEPLIBOBJS :=$(foreach lib, $(PROG_LIBS),lib$(lib)/lib$(lib).so)
PROG_DEPLIBDIRS :=$(foreach lib, $(PROG_LIBS), -Llib$(lib))

PROG_OBJS:=$(PROG:=.o)

ALL_OBJS += $(PROG_OBJS)

all: $(PROG)

# Macro called for each lib available for build
# It creates per libs specific list of objects 
# 
define lib_rules
$(warning lib_rules: $(1): $$($(2)))
$(1): $$($(2))
	$(warning lib_rules:$$(CC) --shared $$^ -o $$@ )
	@echo "(1) $(1)  (2) $(2)"
	$$(CC) --shared $$^ -o $$@

endef

define lib_objs_rules_template

$(warning objs_rules: $(1))
# lsubdir evaluates to the sub directory containing the lib src
# $(1) is the root name of the lib
lname  := lib$(1)
lsubdir:= lib$(1)/

ltarget_shared :=$$(lsubdir)$$(lname)$$(SHARED_LIBSUFF)
ltarget_static :=$$(lsubdir)$$(lname)$$(STATIC_LIBSUFF)

$(warning warn
ltarget_shared :=$$(lsubdir)$$(lname)$$(SHARED_LIBSUFF)
ltarget_static :=$$(lsubdir)$$(lname)$$(STATIC_LIBSUFF)
)
$(warning objs_rules:$(lsubdir), $(lname), $(SHARED_LIBSUFF))
$(warning objs_rules:$(1))

include $(SRC_DIR)/$$(lsubdir)Makefile

# OBJS is a list of .o objects defined in each lib sub Makefiles
# we sort it and store it prefixed with corresponding subdir name
# in a OBJS variable prefixed with libname root pattern
OBJS      := $$(sort $$(OBJS:%=$$(lsubdir)%))
$(1)_OBJS := $$(OBJS)
tmpne:= $(1)_OBJS
$(warning objs_rules:$(tmpne))

# Append per subdir objects to general OBJS list
ALL_OBJS  += $$(OBJS)

#$$(call lib_rules,$$(ltarget_shared),$$($(1)_OBJS))

# Rule to build dynamic library
#$(warning
$(1)_test_dep := $$(foreach l,$$($(1)_LIBS),$$(call $$(lib$$(l)/lib$$(l).so)$(warning $(l))))
lib$(1)/lib$(1).so: $$($(1)_OBJS) $$(foreach l,$$($(1)_LIBS),$$(lib$$(l)/lib$$(l).so))
	$$(CC) --shared $$^ -o $$@
#	)
$(1)_DEPLIBS := $$(foreach l,$$($(1)_LIBS),lib$$(l)/lib$$(l).so)
$(1)_IFLAGS := $$(foreach l,$$($(1)_LIBS),-I./lib$$(l))

lib$(1)/lib$(1).so: $$($(1)_OBJS)  $$($(1)_DEPLIBS)
	$$(CC)  $$(LD_FLAGS) --shared $$^ -o $$@ $$($(1)_IFLAGS) 

#$(warning
$$(ltarget_shared): $$(OBJS)
	$$(CC) --shared $$^ -o $$@
#	)
#$$(ltarget_shared): $$(OBJS)
#	$$(CC) --shared $$^ -o $$@

# Rule to build static library
$$(ltarget_static)$(warning static lib ruel): $$($(1)_OBJS) $(warning static lib preq)
	$(warning AR)
	$$(AR) rcD $$@ $$^
	ranlib $$@

# Append lib targets to general lib objects list
ALL_LIBOBJS += $$(ltarget_shared)
ALL_LIBOBJS += $$(ltarget_static)

endef

# Create libs related variables and rules by calling lib macro
$(foreach lib,$(ALL_LIBS),$(eval $(call lib_objs_rules_template,$(lib))))

$(warning $(ovvcdec_OBJS))
$(warning $(ovvcdec_test_dep))
# Suffix rule on .c objects to build corresponding .o and .d objects
%.o: %.c
	$(CC) -o $@ -c $(patsubst $(SRC_DIR)/%,%,$<) -MMD -MF $(@:.o=.d) -MT $@ $(CFLAGS) -I.

${PROG}: ${PROG_DEPLIBOBJS} $(PROG_OBJS)
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
	$(RM) $(ALL_OBJS) $(ALL_OBJS:.o=.d) $(ALL_LIBOBJS) $(PROG)
