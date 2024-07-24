# If the user has no environment variable
# called METHOD, he gets optimized mode.
ifeq (x$(METHOD),x)
  METHOD := opt
endif

# Paths to libMesh and libCEED
LIBMESH_DIR = "libmesh"
LIBCEED_DIR = "libceed"


# installed libmesh location of libmesh-config script
libmesh_config := $(LIBMESH_DIR)/bin/libmesh-config

# Use the libmesh-config script to determine the usual make variables.
# Note: passing $METHOD along to the libmesh-config script handles the
# case where METHOD is not set and the user does
# make METHOD=dbg foo
# since in this case, METHOD is never set in the environment and thus
# never passed to libmesh-config correctly.
libmesh_CXX      := $(shell METHOD=$(METHOD) $(libmesh_config) --cxx)
libmesh_INCLUDE  := $(shell METHOD=$(METHOD) $(libmesh_config) --include)
libmesh_CPPFLAGS := $(shell METHOD=$(METHOD) $(libmesh_config) --cppflags)
libmesh_CXXFLAGS := $(shell METHOD=$(METHOD) $(libmesh_config) --cxxflags)
libmesh_LIBS     := $(shell METHOD=$(METHOD) $(libmesh_config) --libs)

CEED_FLAGS ?= -I$(LIBCEED_DIR)/include -O -g
CEED_LIBS ?= -Wl,-rpath,$(abspath $(LIBCEED_DIR)/lib) -L$(LIBCEED_DIR)/lib -lceed -lm


# File management variables.
srcfiles 	:= $(wildcard *.C)
opt_executables := $(patsubst %.C, %-opt, $(srcfiles))

.PHONY: clean

# How to build executables. If you have a source file called foo.cc,
# type 'make foo' to build an executable from it.
%: %.C
	@echo "Compiling" $<
	@$(libmesh_CXX) $(libmesh_INCLUDE) $(libmesh_CPPFLAGS) $(libmesh_CXXFLAGS) $(CEED_FLAGS) $(CEED_LDFLAGS) $< -o $@-$(METHOD) $(libmesh_LIBS) $(libmesh_LDFLAGS) $(CEED_LIBS)

# File management rules.
clean:
	@rm -f *~ $(opt_executables)
