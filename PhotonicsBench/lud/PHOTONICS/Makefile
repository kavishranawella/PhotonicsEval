# Makefile: C++ version of lud
# Copyright (c) 2024 University of Virginia
# This file is licensed under the MIT License.
# See the LICENSE file in the root of this repository for more details.

PROJ_ROOT = ../../..
include ${PROJ_ROOT}/Makefile.photonics

# make USE_OPENMP=1
USE_OPENMP ?= 0
ifeq ($(USE_OPENMP),1)
	CXXFLAGS += -fopenmp
endif

LDFLAGS += -pthread  # <-- ADD THIS LINE

EXEC := lud.out
SRC := lud.cpp

CONDA_PREFIX ?= $(shell echo $$CONDA_PREFIX)
VIENNACL_DIR := $(HOME)/local_builds/viennacl-dev

USE_CUDA ?= 1        # <-- Enable GPU support (set to 0 to disable)
ifeq ($(USE_CUDA),1)
	CUDA_PATH ?= /usr/local/cuda-12.4
	CXXFLAGS  += -DUSE_GPU_SOLVER=1 -I$(CUDA_PATH)/include -I$(VIENNACL_DIR) -I$(CONDA_PREFIX)/include -O2 -std=c++17 -DVIENNACL_WITH_OPENCL -DCL_TARGET_OPENCL_VERSION=120
	LDFLAGS   += -L$(CUDA_PATH)/lib64 -L$(CONDA_PREFIX)/lib -lOpenCL -ldl -lcublas -lcudart -lnpps -lnppc
else
	CXXFLAGS  += -DUSE_GPU_SOLVER=0
endif

COMPILE_WITH_JPEG ?= 0
ifeq ($(COMPILE_WITH_JPEG), 1)
    CXXFLAGS += -DCOMPILE_WITH_JPEG
    # Add linker flags for libjpeg
	LDFLAGS += -ljpeg
endif

perf dramsim3_integ: $(EXEC)

debug: CXXFLAGS += -g -O0 -fno-omit-frame-pointer
debug: $(EXEC)

$(EXEC): $(SRC) $(DEPS)
	$(CXX) $< $(CXXFLAGS) $(LDFLAGS) -o $@

clean:
	rm -rf $(EXEC) *.dSYM
