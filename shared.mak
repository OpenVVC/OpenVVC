
	
DEPLIBS := $(foreach l,$($(NAME)_LIBS),$(BUILD_DIR)lib$(l)/lib$(l).so)
IFLAGS  := $(foreach l,$($(NAME)_LIBS),-I./lib$(l))

# Rule to build dynamic library
$(BUILD_DIR)$(ltarget_shared): $(OBJS) $(DEPLIBS)
	$(CC)  $(LD_FLAGS) $(IFLAGS) --shared $^ -o $@ 
	
# FIXME do we need dependencies besides OBJS in static lib

# Rule to build static library
$(BUILD_DIR)$(ltarget_static): $(OBJS)
	$(AR) rcD $@ $^
	ranlib $@
