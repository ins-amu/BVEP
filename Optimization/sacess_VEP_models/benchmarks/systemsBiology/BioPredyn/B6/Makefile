# CODE VERSION ######################################################

VERSION = 10.0
USRFLAG  = -DUSR=$(USER)
HOSTFLAG = -DMACHINE=$(HOST)
COMPFLAG = -DCOMPILER=$(CC)
FLAGFLAG = -DFLAGS=optimized
VERSFLAG = -DVERS=$(VERSION)
VFLAGS = $(USRFLAG) $(HOSTFLAG) $(COMPFLAG) $(FLAGFLAG) $(VERSFLAG) 

MAKEFILEPATH:=$(abspath $(lastword $(MAKEFILE_LIST)))
TRUNKDIR:=$(patsubst %/,%,$(dir $(MAKEFILEPATH)))


# set default architecture
BITS = 64

MPIFLAGS+= $(CCFLAGS) -DMPI 
DEBUGFLAGS+= $(DEBUGFLAGS) -DMPI
PROFILEFLAGS+= $(PROFILEFLAGS) -DMPI

ifeq ($(BITS),64)
	#64 bit 
	SUNDIALS = $(TRUNKDIR)/lib64/sundials-2.5.0
	GSL = $(TRUNKDIR)/lib64/gsl-1.14
	LGSL = -L$(GSL)/lib
        export LGSL
	M = m64
endif

LSUNDIALS = -L$(SUNDIALS)/lib 

CCFLAGS+=  # -DHAVE_SSE2 -$(M) -fPIC -DPIC -static
PROFILEFLAGS+= -g #-pg #2 -DHAVE_SSE2
LIBS+= -O3  -lm -lgsl -lgslcblas -lsundials_cvode -lsundials_nvecserial $(LSUNDIALS) $(LGSL) 
FLIBS+= $(LIBS)
KCC+= $(CC)
KFLAGS+= $(CCFLAGS)

ifdef DEBUG
	CCFLAGS+= $(DEBUGFLAGS)
	FLAGFLAG+= -DFLAGS=\"debugging\"
else
	DEBUG+= "Off"
endif

ifdef PROFILE
	CCFLAGS+= $(PROFILEFLAGS)
	FLAGFLAG+= -DFLAGS=\"profiling\"
	KCC+= $(CC)
	KFLAGS+=
else
	PROFILE+= "Off"
endif

ifdef USEPROFILE
        CCFLAGS+= $(USEPROFFLAGS)
endif

INCLUDES+= -I$(GSL)/include -I. -I$(TRUNKDIR)/util  -I$(ICC_PATH)/compiler/include -I$(SUNDIALS)/include -I$(TRUNKDIR)/fly 
export INCLUDES
CFLAGS+= $(CCFLAGS) $(INCLUDES) $(CLFLAGS)
export CFLAGS
export VFLAGS
export CC
export KCC
export MPICC
export KFLAGS
export LIBS
export FLIBS
export MPIFLAGS
export FLYEXECS
#define targets
.PHONY: ggn fly util 

#ggn:	ggn.o 
#	$(CC) -o ggn $(CFLAGS) ./fly/*.o ./util/*.o ggn.o $(FLIBS) 
	
ggn.o:	objects
	$(CC) -c  $(CFLAGS) $(VFLAGS) ggn.c

objects: util
	 cd fly && $(MAKE) objects

fly:	util
	cd fly && $(MAKE) all

util:
	cd util && $(MAKE)

deps: 
	cd fly && $(MAKE) -f basic.mk Makefile && chmod +w Makefile

clean:
	rm -f core* *.o *.il
	rm -f */core* util/*.o fly/*.o */*.il
	rm -f ggn
	

veryclean:	clean
	rm -f */*.slog */*.pout */*.uout
	rm -f fly/Makefile
	rm -f fly/zygotic.cmp.c

help:
	@echo "make: this is the Makefile for fly code"
	@echo "      always 'make deps' first after a 'make veryclean'"
	@echo ""
	@echo "      the following targets are available:"
	@echo "      util:      make object files in the util directory only"
	@echo "      fly:       compile the fly code (which is in 'fly')"
	@echo "      clean:     gets rid of cores and object files"
	@echo "      veryclean: gets rid of executables and dependencies too"
	@echo ""
	@echo "      your current settings are:"   
	@echo "      compiler:  $(CC)"
	@echo "      flags:     $(CFLAGS)"
	@echo "      debugging: $(DEBUG)"
	@echo "      profiling: $(PROFILE)"
	@echo "      os type:   $(OSTYPE)"
	@echo ""

