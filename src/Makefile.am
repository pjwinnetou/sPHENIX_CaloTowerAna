AUTOMAKE_OPTIONS = foreign

AM_CPPFLAGS = \
  -I$(includedir) \
  -I$(OFFLINE_MAIN)/include \
  -I$(ROOTSYS)/include

AM_LDFLAGS = \
  -L$(libdir) \
  -L$(OFFLINE_MAIN)/lib

pkginclude_HEADERS = \
  ClusterIso.h \
  CaloAna.h 

ROOTDICTS = 

pcmdir = $(libdir)
nobase_dist_pcm_DATA = 

lib_LTLIBRARIES = \
  libcaloana_io.la \
  libcaloana.la

libcaloana_io_la_SOURCES = \
  $(ROOTDICTS) 

libcaloana_la_SOURCES = \
	ClusterIso.cc \
  CaloAna.cc 

libcaloana_io_la_LIBADD = \
  -lphool

libcaloana_la_LIBADD = \
  libcaloana_io.la \
	-lcalo_io \
	-ljetbase_io \
	-lfun4all \
	-lg4dst \
	-lg4eval \
	-lg4jets \
	-lmbd_io \
	-lcentrality \
	-lcalotrigger \
  -lSubsysReco -lphg4hit -lffaobjects
	
%_Dict.cc: %.h %LinkDef.h
	rootcint -f $@ @CINTDEFS@ $(DEFAULT_INCLUDES) $(AM_CPPFLAGS) $^

#just to get the dependency
%_Dict_rdict.pcm: %_Dict.cc ;

BUILT_SOURCES = testexternals.cc

noinst_PROGRAMS = \
  testexternals_io \
  testexternals

testexternals_io_SOURCES = testexternals.cc
testexternals_io_LDADD   = libcaloana_io.la

testexternals_SOURCES = testexternals.cc
testexternals_LDADD   = libcaloana.la

testexternals.cc:
	echo "//*** this is a generated file. Do not commit, do not edit" > $@
	echo "int main()" >> $@
	echo "{" >> $@
	echo "  return 0;" >> $@
	echo "}" >> $@

clean-local:
	rm -f $(BUILT_SOURCES)
