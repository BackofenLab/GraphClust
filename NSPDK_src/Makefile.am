bin_PROGRAMS = NSPDK

NSPDK_SOURCES=NSPDK.cc vectors.cpp vectors.h gzstream.cpp gzstream.h	\
              GraphClass.h GraphClass.cc Histogram.h Histogram.cc	\
              NSPDK_FeatureGenerator.h NSPDK_FeatureGenerator.cc	\
              FlagsService.h FlagsService.cc Utility.h Utility.cc	\
              BaseGraphClass.h BaseGraphClass.cc FeatureGenerator.h \
              wrapper.h

LDADD = -lz -lm

## -mtune=native -march=native optimzes for the machine on which make is run
## -std=c++0x is needed if we use not the tr1 unordered_map include but just unordered_map
# -mtune=native -march=native
AM_CPPFLAGS = -DLOSS=1 -fopenmp
