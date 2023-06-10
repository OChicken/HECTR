#!/usr/bin/make -f
#
# Copyright (C) shouran.ma@rwth-aachen.de
#
# This file is part of HECTR.
#
# GPQHE is free software; you can redistribute it and/or modify
# it under the terms of the GNU Lesser General Public License as
# published by the Free Software Foundation; either version 2.1 of
# the License, or (at your option) any later version.
#
# GPQHE is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU Lesser General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public
# License along with this program; if not, see <http://www.gnu.org/licenses/>.

CC     ?= /usr/bin/gcc
CFLAGS += -Wall -Wextra -Wpedantic -Wredundant-decls -Wshadow -Wpointer-arith -Og -fomit-frame-pointer
LIBDIR ?= $(ROOT)/lib

all:
	cd libpmu && $(MAKE) && mv libpmu.so ../lib && cd ..
	cd GPQHE  && mkdir -p lib && cd src && $(MAKE) && cd .. && mv lib/libgpqhe.so ../lib && rm -r lib && cd ..
	cd src    && $(MAKE)