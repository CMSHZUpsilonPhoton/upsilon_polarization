### definition of the compiler options ###
#       -I location of directory containing include files
#       -L location of directory containing libraries
#       -lname include the library from -L location called libname.a
#          -lg2c is the library containing info on converting fortran to C
#          -lf   is the library containing the intrinsic for HPUX only.
#       -shared make a shared library as output
#       -fPIC produce position independent code
#        necessary on some platforms (including HPUX) for -shared
#       -fpic ^^ same(?)
#       -O optimizes
#       -g produces output for the debugger
#       -pg produces output for gprof profiler
#       note: if you want to see all warnings and ensure ansi standard
#             compatibility, use:
#             -pipe -ansi -pedantic -fnonnull-objects \
#             -W -Wall -Wwrite-strings -Wpointer-arith -Wnested-externs \
#             -Woverloaded-virtual -Wbad-function-cast -fnonnull-objects
#       The proper order for cernlib libraries is:
#       -lpawlib -lgraflib -lgrafX11 -lmathlib -lkernlib -lpacklib -ljetset74
#
# makefile syntax:
#        for target thedir/target.suf from source anotherdir/source.suf2
#        ${*D}  = thedir
#        ${*F}  = target
#        $^     = name of all prerequisites
#        $+     = like $^ with duplicated prerequisites if mentioned more than once
#        $*     = thedir/target
#        $@     = thedir/target.suf
#        $<     = anotherdir/source.suf2
#

# define source files
SOURCES := $(wildcard upsilon_polarization/src/*.cxx)

# define header, object and dependency files based on source files
HEADERS := $(SOURCES:.cxx=.h)
OBJECTS := $(SOURCES:.cxx=.o)
DEPENDS := $(SOURCES:.cxx=.d)


########################################
# compiler and flags

# ROOT libraries:
ROOT_CFLAGS  := $(shell root-config --cflags)
ROOT_LDFLAGS := $(shell root-config --ldflags)
ROOT_GLIBS   := $(shell root-config --libs)

CC := g++
CFLAGS	:= -O3 -Wall -fPIC -Wextra

LD := g++
LDFLAGS := -O3fast -lz -shared -fPIC

# Debug flags?
ifdef DEBUG
CFLAGS = -Og -Wall -fPIC -g -pg -fprofile-generate
LDFLAGS = -Og -g -pg -lz -fprofile-generate
endif

CFLAGS  += -I. $(ROOT_CFLAGS) 
LDFLAGS += $(ROOT_LDFLAGS) $(ROOT_GLIBS) 


########################################
# define colors for better overview

RED=$(shell tput setaf 1)
YELLOW=$(shell tput setaf 2)
GREEN=$(shell tput setaf 3)
BOLD=$(shell tput bold)
NORMAL=$(shell tput sgr0)

########################################
# targets

# target: dependency
#	@shellcommand
#	command

all: upsilon_polarization

clean:
	rm -f $(OBJECTS) $(DEPENDS) upsilon_polarization/lib/*
	rm -rf lib
	rm -rf validation_outputs

test:
	pytest --verbosity=99

validation: validation_extreme validation_full validation_with_MC_sample

validation_extreme: 
	mkdir -p validation_outputs
	rm -rf validation_outputs/extreme*
	python3 -m upsilon_polarization.validation.val_extreme_scenarios

validation_full: 
	python3 -m upsilon_polarization.validation.val_full_method

validation_with_MC_sample: 
	mkdir -p validation_outputs
	root -l -b -q "upsilon_polarization/validation/validation_Z.C+(1)"

upsilon_polarization: $(OBJECTS)
	echo "Building shared library $@.so ..."
	$(LD) $^ $(LDFLAGS) -o $@.so
	echo "$@.so done"
	mkdir -p upsilon_polarization/lib
	mv $@.so upsilon_polarization/lib/.

########################################
# rules
%.o : %.cxx
	$(CC) -MD -MP $(CFLAGS) -g -c -o $@ $<

-include $(DEPENDS)


