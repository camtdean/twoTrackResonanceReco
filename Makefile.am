AUTOMAKE_OPTIONS = foreign

AM_CPPFLAGS = \
  -I$(includedir) \
  -isystem$(OFFLINE_MAIN)/include \
  -isystem$(ROOTSYS)/include

AM_LDFLAGS = \
  -L$(libdir) \
  -L$(OFFLINE_MAIN)/lib \
  -L$(OFFLINE_MAIN)/lib64 \
  `root-config --libs`

pkginclude_HEADERS = \
  twoTrackResonanceReco.h

lib_LTLIBRARIES = \
  libtwoTrackResonanceReco.la

libtwoTrackResonanceReco_la_SOURCES = \
  twoTrackResonanceReco.cc

libtwoTrackResonanceReco_la_LIBADD = \
  -lphool \
  -lSubsysReco \
  -lfun4all \
  -ltrackbase_historic \
  -ltrackbase_historic_io \
  -ltrackbase \
  -lglobalvertex \
  -lg4detectors \
  -lffamodules

BUILT_SOURCES = testexternals.cc

noinst_PROGRAMS = \
  testexternals

testexternals_SOURCES = testexternals.cc
testexternals_LDADD   = libtwoTrackResonanceReco.la

testexternals.cc:
	echo "//*** this is a generated file. Do not commit, do not edit" > $@
	echo "int main()" >> $@
	echo "{" >> $@
	echo "  return 0;" >> $@
	echo "}" >> $@

clean-local:
	rm -f $(BUILT_SOURCES)
