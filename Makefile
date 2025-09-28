# Makefile for PIMeval / PIMbench Framework
# Copyright (c) 2024 University of Virginia
# This file is licensed under the MIT License.
# See the LICENSE file in the root of this repository for more details.

PIMLIBDIR := libpimeval
PHOTONICSLIBDIR := libphotonicseval
BITSERIALDIR := bit-serial
PIMAPPDIR := PIMbench
PHOTONICSAPPDIR := PhotonicsBench
TESTDIR := misc-bench tests
ALLDIRS := $(PIMLIBDIR) $(PHOTONICSLIBDIR) $(BITSERIALDIR) $(PIMAPPDIR) $(PHOTONICSAPPDIR) $(TESTDIR)

# Handle dependency between lib and apps to support make -j
DEP_LIBPIMEVAL := $(PIMLIBDIR)/lib/libpimeval.a
DEP_LIBPHOTONICSEVAL := $(PHOTONICSLIBDIR)/lib/libphotonicseval.a

.PHONY: debug perf dramsim3_integ clean $(ALLDIRS)
.DEFAULT_GOAL := perf

debug: $(ALLDIRS)
	@echo "\nINFO: Built PHOTONICSeval Simulator with target = debug\n"

perf: $(ALLDIRS)
	@echo "\nINFO: Built PHOTONICSeval Simulator with target = perf\n"

dramsim3_integ: $(ALLDIRS)
	@echo "\nINFO: Built PHOTONICSeval Simulator with target = dramsim3_integ\n"

clean: $(ALLDIRS)

# Run make with PHOTONICS_SIM_TARGET=<PhotonicsDeviceEnum> to override default simulation target
PIM_SIM_TARGET ?= PIM_DEVICE_NONE
PHOTONICS_SIM_TARGET ?= PHOTONICS_DEVICE_NONE

# Run make with USE_OPENMP=1 to enable OpenMP in some apps
USE_OPENMP ?= 0

# Run make with COMPILE_WITH_JPEG=0 to disable compilation with JPEG. JPEG compilation is needed for VGG apps.
COMPILE_WITH_JPEG ?= 0

$(DEP_LIBPIMEVAL) $(PIMLIBDIR):
	$(MAKE) -C $(PIMLIBDIR) $(MAKECMDGOALS) PIM_SIM_TARGET=$(PIM_SIM_TARGET) USE_OPENMP=$(USE_OPENMP) COMPILE_WITH_JPEG=$(COMPILE_WITH_JPEG)

$(DEP_LIBPHOTONICSEVAL) $(PHOTONICSLIBDIR):
	$(MAKE) -C $(PHOTONICSLIBDIR) $(MAKECMDGOALS) PHOTONICS_SIM_TARGET=$(PHOTONICS_SIM_TARGET) USE_OPENMP=$(USE_OPENMP) COMPILE_WITH_JPEG=$(COMPILE_WITH_JPEG)

$(BITSERIALDIR) $(PIMAPPDIR) $(TESTDIR): $(DEP_LIBPIMEVAL)
	$(MAKE) -C $@ $(MAKECMDGOALS) PIM_SIM_TARGET=$(PIM_SIM_TARGET) USE_OPENMP=$(USE_OPENMP) COMPILE_WITH_JPEG=$(COMPILE_WITH_JPEG)

$(PHOTONICSAPPDIR): $(DEP_LIBPHOTONICSEVAL)
	$(MAKE) -C $@ $(MAKECMDGOALS) PHOTONICS_SIM_TARGET=$(PHOTONICS_SIM_TARGET) USE_OPENMP=$(USE_OPENMP) COMPILE_WITH_JPEG=$(COMPILE_WITH_JPEG)

