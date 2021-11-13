SUBDIRS = libovvc examples
ACLOCAL_AMFLAGS = -Im4

all: update-build-number all-am

.PHONY: update-build-number
update-build-number:
	$(AM_V_GEN)$(srcdir)/version.sh $(srcdir)/VERSION libovvc/ovversion.h ./
