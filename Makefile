#
# CDDL HEADER START
#
# The contents of this file are subject to the terms of the Common Development
# and Distribution License Version 1.0 (the "License").
#
# You can obtain a copy of the license at
# http://www.opensource.org/licenses/CDDL-1.0.  See the License for the
# specific language governing permissions and limitations under the License.
#
# When distributing Covered Code, include this CDDL HEADER in each file and
# include the License file in a prominent location with the name LICENSE.CDDL.
# If applicable, add the following below this CDDL HEADER, with the fields
# enclosed by brackets "[]" replaced with your own identifying information:
#
# Portions Copyright (c) [2016] [Nikhil Chandra Admal and Giacomo Po]. All rights reserved.
#
# CDDL HEADER END
#

#
# Copyright (c) 2013--2014, Regents of the University of Minnesota.
# All rights reserved.
#
# Contributors:
#    Nikhil Chandra Admal
#


ifneq (clean,$(MAKECMDGOALS))
  ifeq ($(wildcard ../Makefile.KIM_Config_Helper),)
    KIM_CONFIG_HELPER = kim-api-build-config
  else
    include ../Makefile.KIM_Config_Helper
  endif
  ifeq ($(shell $(KIM_CONFIG_HELPER) --version 2> /dev/null),)
    $(error $(KIM_CONFIG_HELPER) utility is not available.  Something is wrong with your KIM API package setup)
  endif
endif

CC           = $(shell $(KIM_CONFIG_HELPER) --cc)
CXX          = $(shell $(KIM_CONFIG_HELPER) --cxx)
FC           = $(shell $(KIM_CONFIG_HELPER) --fc)
LD           = $(shell $(KIM_CONFIG_HELPER) --ld)
INCLUDES     = $(shell $(KIM_CONFIG_HELPER) --includes)
CFLAGS       = $(shell $(KIM_CONFIG_HELPER) --cflags)
CXXFLAGS     = $(shell $(KIM_CONFIG_HELPER) --cxxflags)
FFLAGS       = $(shell $(KIM_CONFIG_HELPER) --fflags)
LDFLAGS      = $(shell $(KIM_CONFIG_HELPER) --ldflags)
LDLIBS       = $(shell $(KIM_CONFIG_HELPER) --ldlibs)
OBJONLYFLAG  = $(shell $(KIM_CONFIG_HELPER) --objonlyflag)
OUTPUTINFLAG = $(shell $(KIM_CONFIG_HELPER) --outputinflag)

SRC = mod_global.F03  \
	  mod_lattice.F03 \
	  mod_crystal.F03 \
	  mod_atomistic.F03 \
	  mod_kim.F03	  \
	  mod_equilibrium.F03 \
	  mod_fsge.F03 \
	  mod_matrix.F03 \
	  main.F03

FOBJ = mod_global.o mod_fsge.o mod_matrix.o mod_atomistic.o mod_lattice.o mod_crystal.o mod_kim.o mod_equilibrium.o main.o

TEST_NAME := runner

all: $(TEST_NAME)

$(TEST_NAME): $(FOBJ)
	$(LD) $(LDFLAGS) $(FOBJ) $(LDLIBS) $(OUTPUTINFLAG) $@

main.o : main.F03 mod_lattice.o mod_crystal.o mod_equilibrium.o mod_global.o
mod_global.o : mod_global.F03
mod_matrix.o : mod_matrix.F03 mod_global.o
mod_lattice.o : mod_lattice.F03 mod_global.o
mod_crystal.o : mod_crystal.F03 mod_lattice.o mod_global.o
mod_atomistic.o : mod_atomistic.F03 mod_global.o mod_fsge.o mod_matrix.o mod_fsge.o
mod_kim.o : mod_kim.F03 mod_atomistic.o mod_global.o
mod_equilibrium.o : mod_equilibrium.F03 mod_kim.o mod_crystal.o mod_atomistic.o mod_global.o
mod_fsge : mod_fsge.F03 mod_global.o

.PHONY: all clean

%.o:%.F03
	$(FC) $(INCLUDES) $(FFLAGS) $(OBJONLYFLAG) $<

clean:
	rm -f $(TEST_NAME) *.o kim.log *.mod
