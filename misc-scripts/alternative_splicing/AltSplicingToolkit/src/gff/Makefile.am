# noinst_
noinst_LIBRARIES = libgff.a

libgff_a_SOURCES = BioMartGffHandler.h BioMartGffHandler.cpp GffEventModel.h GffEventModel.cpp GffParser.h GffParser.cpp GffSimpleHandler.h GffSimpleHandler.cpp  SplicingEventGffGenerator.h SplicingEventGffGenerator.cpp

libgff_a_DEPENDENCIES = @top_srcdir@/src/as/libas.a

AM_CPPFLAGS = -I@top_srcdir@/src
#libgff_a__LDFLAGS = -static
