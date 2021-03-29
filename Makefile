include config.mak

# Set defaults
VERBOSITY?=0
LD_FLAGS?=-lpthread
SRC_FOLDER:=libovvc/


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

# Quick hack to avoid missing / in builddir 
BUILDDIR:=$(BUILDDIR)/

# Find Sources
include $(SRC_FOLDER)/libobj.mak
LIB_SRC:=$(addprefix $(SRC_FOLDER),$(LIB_SRC))
LIB_HEADER:=$(addprefix $(SRC_FOLDER),$(LIB_HEADER))
LIB_OBJ:=$(addprefix $(BUILDDIR),$(LIB_SRC:%.c=%.o))
LIB_FILE:=$(LIB_HEADER) $(LIB_SRC)

include $(SRC_FOLDER)/$(ARCH)/$(ARCH)obj.mak
$(ARCH)_LIB_SRC:=$(addprefix $($(ARCH)_SRC_FOLDER),$($(ARCH)_LIB_SRC))
$(ARCH)_LIB_OBJ:=$(addprefix $(BUILDDIR),$($(ARCH)_LIB_SRC:%.c=%.o))
BUILDDIR_TYPE_ARCH:=$(addprefix $(BUILDDIR), $($(ARCH)_SRC_FOLDER))

LIB_NAME:= libovvc

PROG=examples/dectest

ALL_OBJS=$(LIB_OBJ) $(addprefix $(BUILDDIR),$(addsuffix .o, $(PROG))) $($(ARCH)_LIB_OBJ)


.PHONY: all test version libs examples

all: libs examples

test:
	./CI/checkMD5.sh $(TESTSTREAMSDIR) $(BUILDDIR)$(PROG)

version:
	$(AT)./version.sh VERSION $(SRC_FOLDER)$(LIB_VERSION_HEADER)

libs: version $(BUILDDIR)$(LIB_NAME)$(STATIC_LIBSUFF) $(BUILDDIR)$(LIB_NAME)$(SHARED_LIBSUFF)

examples: $(BUILDDIR)$(PROG) $(BUILDDIR)$(PROG)_stat

$(SRC_FOLDER)$(LIB_VERSION_HEADER): VERSION
	$(AT)./version.sh $< $@

$(BUILDDIR)$(PROG):  $(BUILDDIR)$(PROG).o $(BUILDDIR)$(LIB_NAME)$(SHARED_LIBSUFF)
	$(CC) $^ -o $@ $(LD_FLAGS)

$(BUILDDIR)$(PROG)_stat:  $(BUILDDIR)$(PROG).o $(BUILDDIR)$(LIB_NAME)$(STATIC_LIBSUFF)
	$(CC) $^ -o $@ $(LD_FLAGS)


$(BUILDDIR)$(LIB_NAME)$(STATIC_LIBSUFF): $(LIB_OBJ) $($(ARCH)_LIB_OBJ)
	$(AT)./version.sh VERSION $(SRC_FOLDER)$(LIB_VERSION_HEADER)
	$(AR) rc $@ $^
	$(RANLIB) $@

$(BUILDDIR)$(LIB_NAME)$(SHARED_LIBSUFF): $(LIB_OBJ) $($(ARCH)_LIB_OBJ)
	$(AT)./version.sh VERSION $(SRC_FOLDER)$(LIB_VERSION_HEADER)
	$(CC) -shared $^ -o $@ $(LD_FLAGS)

$(BUILDDIR_TYPE_ARCH)%_sse.o: $($(ARCH)_SRC_FOLDER)%_sse.c
	$(AT)mkdir -p $(@D)
	$(CC) -c $< -o $@ -MMD -MF $(@:.o=.d) -MT $@ $(CFLAGS) $(SSE_CFLAGS) -I$(SRC_FOLDER)

$(BUILDDIR)%.o: %.c
	$(AT)mkdir -p $(@D)
	$(CC) -c $< -o $@ -MMD -MF $(@:.o=.d) -MT $@ $(CFLAGS) -I$(SRC_FOLDER)


.PHONY: install install-shared install-headers install-pkgconfig

install-shared: $(BUILDDIR)$(LIB_NAME)$(SHARED_LIBSUFF)
	$(AT)mkdir -p $(INSTALL_LIB)
	cp $< $(INSTALL_LIB)/$(<F)

install-headers: $(LIB_HEADER)
	$(AT)mkdir -p $(INSTALL_INCLUDE)
	cp $^ $(INSTALL_INCLUDE)/

install: install-shared install-headers install-pkgconfig

install-pkgconfig: version
	$(AT)mkdir -p $(INSTALL_PKGCONFIG)
	cp libopenvvc.pc $(INSTALL_PKGCONFIG)/libopenvvc.pc

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
	$(AT)rm -f $(ALL_OBJS) $(ALL_OBJS:.o=.d) $(addprefix $(BUILDDIR),$(PROG)*) $(BUILDDIR)$(LIB_NAME)$(STATIC_LIBSUFF) $(BUILDDIR)$(LIB_NAME)$(SHARED_LIBSUFF)
