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
# Portions Copyright (c) [yyyy] [name of copyright owner]. All rights reserved.
#
# CDDL HEADER END
#

#
# Copyright (c) 2013--2014, Regents of the University of Minnesota.
# All rights reserved.
#
# Contributors:
#    Ryan S. Elliott
#

#
# Release: This file is part of the kim-api-v1.7.1 package.
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


.PHONY: all clean

TEST_NAME := test_elastic_constants

all: $(TEST_NAME)

$(TEST_NAME).o: $(TEST_NAME).F03 Makefile 
	$(FC) $(INCLUDES) $(FFLAGS) $(OBJONLYFLAG) $<

$(TEST_NAME): $(TEST_NAME).o
	$(LD) $(LDFLAGS) $< $(LDLIBS) $(OUTPUTINFLAG) $@

clean:
	rm -f $(TEST_NAME) $(TEST_NAME).o kim.log *.mod
