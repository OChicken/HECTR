#!/usr/bin/make -f
#
# Copyright (C) shouran.ma@rwth-aachen.de
#
# This file is part of HECTR.
#
# HECTR is free software; you can redistribute it and/or modify
# it under the terms of the GNU Lesser General Public License as
# published by the Free Software Foundation; either version 2.1 of
# the License, or (at your option) any later version.
#
# HECTR is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU Lesser General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public
# License along with this program; if not, see <http://www.gnu.org/licenses/>.

ROOT ?= ..

HEADERS = \
	config.h \
	hectr.h
SOURCES = \
	ode.c \
	cstr.c \
	matrices.c \
	ctr.c \
	dare.c \
	dlqr.c \
	dlqe.c \
	mpc.c \
	hempc.c \
	quadprog.c

libhectr.so: $(HEADERS) $(SOURCES)
	$(CC) $(CFLAGS) $(SOURCES) -fPIC -shared -llapacke -lm -o $(LIBDIR)/$@

include $(ROOT)/Makefile
