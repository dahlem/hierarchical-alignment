# Copyright (C) 2015 Dominik Dahlem <Dominik.Dahlem@gmail.com>
#  
# This file is free software; as a special exception the author gives
# unlimited permission to copy and/or distribute it, with or without 
# modifications, as long as this notice is preserved.
# 
# This program is distributed in the hope that it will be useful, but
# WITHOUT ANY WARRANTY, to the extent permitted by law; without even the
# implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.

bin_PROGRAMS = ha

ha_SOURCES =                                                                 \
	main.cc                                                              \
	CL.cc

ha_CPPFLAGS =                                                                \
	$(OPENMP_CXXFLAGS)                                                   \
	$(BOOST_CPPFLAGS)                                                    \
	-I$(top_srcdir)/src/common/includes                                  \
	-I$(top_srcdir)/src/alignment/includes                               \
	-I./includes

ha_LDADD =                                                                   \
	$(BOOST_FILESYSTEM_LIB)                                              \
	$(BOOST_PROGRAM_OPTIONS_LIB)                                         \
	$(BOOST_SYSTEM_LIB)

ha_LDFLAGS =                                                                 \
	$(BOOST_LDFLAGS)

include_HEADERS = includes
MAINTAINERCLEANFILES = Makefile.in

check-syntax:
	gcc -o nul -S ${ha_SOURCES}
