# If the user has no environment variable
# called METHOD, he gets optimized mode.
ifeq (x$(METHOD),x)
  METHOD := opt
endif

# Paths to libMesh and libCEED
LIBMESH_DIR ?= libmesh
LIBCEED_DIR ?= libCEED
INCLUDE_DIR = include
QFUNCION_DIR = qfunctions


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

CEED_FLAGS ?= -I$(LIBCEED_DIR)/include -g
CEED_LIBS ?= -Wl,-rpath,$(abspath $(LIBCEED_DIR)/lib) -L$(LIBCEED_DIR)/lib -lceed -lm

# Add the project's include directory to the CXXFLAGS
libmesh_CXXFLAGS += -I$(INCLUDE_DIR)
# Add the project's Qfunctions directory to the CXXFLAGS
libmesh_CXXFLAGS += -I$(QFUNCION_DIR)

# File management variables.
SRC_DIR := src
srcfiles 	:= $(wildcard $(SRC_DIR)/*.C)
objfiles := $(patsubst $(SRC_DIR)/%.C, $(SRC_DIR)/%.o, $(srcfiles))
#executables := $(patsubst $(SRC_DIR)/%.C, %, $(srcfiles))
executables := main


.PHONY: all clean

# Default target: build all executables
all: $(executables)


# Link object files to create the executable
$(executables): $(objfiles)
	@echo "Linking" $@
	$(libmesh_CXX) -o $@-$(METHOD) $^  $(libmesh_LIBS) $(CEED_LIBS)

# Compile source files to object files
$(SRC_DIR)/%.o: $(SRC_DIR)/%.C
	@echo "Compiling" $<
	$(libmesh_CXX) $(libmesh_INCLUDE) $(libmesh_CPPFLAGS) $(libmesh_CXXFLAGS) $(CEED_FLAGS) -c $< -o $@

# Clean up the generated files
clean:
	rm -f $(SRC_DIR)/*.o $(executables)-$(METHOD)

# # How to build executables. If you have a source file called foo.cc,
# # type 'make foo' to build an executable from it.
# %: $(SRC_DIR)/%.C
# 	@echo "Compiling" $<
# 	$(libmesh_CXX) $(libmesh_INCLUDE) $(libmesh_CPPFLAGS) $(libmesh_CXXFLAGS) $(CEED_FLAGS) $(CEED_LDFLAGS) $< -o $@-$(METHOD) $(libmesh_LIBS) $(libmesh_LDFLAGS) $(CEED_LIBS)

# # File management rules.
# clean:
# 	@rm -f *~ $(opt_executables)
