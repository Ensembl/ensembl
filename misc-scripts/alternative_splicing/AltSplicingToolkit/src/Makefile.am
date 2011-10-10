
bin_PROGRAMS = altSpliceFinder

bin_CPPFLAGS = -L@top_srcdir@/src
LDADD = -lAltSplicing -lm

altSpliceFinder_SOURCES = Main.cpp Constants.h
altSpliceFinder_CPPFLAGS = @CPPFLAGS@ -L@top_srcdir@/src

# general use to compile individuals objects
AM_CPPFLAGS = -I@top_srcdir@/src -L@top_srcdir@/src

# -------------------------------
# AltSplicingToolkit library configuration
# -------------------------------

lib_LIBRARIES = libAltSplicing.a

libAltSplicing_a_SOURCES = \
	as/GeneFeature.cpp\
	as/SplicingEventMatrix.cpp\
	as/Transcript.cpp\
	as/Feature.cpp\
	as/Coordinates.cpp\
	as/SplicingEventContainer.cpp\
	as/TranscriptFeature.cpp\
	as/Gene.cpp\
	as/SplicingEvent.cpp\
	as/RegionChunk.cpp\
	as/ExonChunk.cpp\
	gff/BioMartGffHandler.cpp\
	gff/GffParser.cpp\
	gff/SplicingEventGffGenerator.cpp\
	gff/GffEventModel.cpp\
	gff/GffSimpleHandler.cpp\
	util/StringUtil.cpp\
	util/Logger.cpp

libAltSplicing_a_CPPFLAGS = $(AM_CPPFLAGS)

# -------------------------------
# AltSplicingToolkit package distribution
# -------------------------------

pkginclude_HEADERS = \
	as/Coordinates.h\
	as/GeneFeature.h\
	as/SplicingEventMatrix.h\
	as/Transcript.h\
	as/Feature.h\
	as/SplicingEventContainer.h\
	as/TranscriptFeature.h\
	as/Gene.h\
	as/SplicingEvent.h\
	as/RegionChunk.h\
	as/ExonChunk.h\
	gff/BioMartGffHandler.h\
	gff/GffParser.h\
	gff/SplicingEventGffGenerator.h\
	gff/GffEventModel.h\
	gff/GffSimpleHandler.h\
	util/StringUtil.h

EXTRA_DIST = $(pkginclude_HEADERS)
