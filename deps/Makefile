# Figure out OS and architecture
OS = $(shell uname)
ifeq ($(OS), MINGW32_NT-6.1)
  OS=WINNT
endif

# file extensions
ifeq ($(OS), WINNT)
  SHLIB_EXT = dll 
else ifeq ($(OS), Darwin)
  SHLIB_EXT = dylib
else
  SHLIB_EXT = so
endif




all: TSSM/OPENMP/libtssm.$(SHLIB_EXT)


TSSM/Makefile.system:
	if [ ! -d "TSSM" ];then     \
		git clone git://github.com/HaraldHofstaetter/TSSM.git; \
	fi


export JULIA_BUILD = 1
TSSM/OPENMP/libtssm.$(SHLIB_EXT): TSSM/Makefile.system
	make -C TSSM/OPENMP all
	
