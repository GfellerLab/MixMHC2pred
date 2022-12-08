CPPFLAGS = -std=gnu++17
CPPFLAGS += -O3
# Add some compilation instructions here.
# -O2 for optimisation of level 2; -O3 for even higher optimisation
CPPFLAGS += -Wall # To generate more compilation warnings of possible issues.

LDLIBS = -lm # Write the libraries to add, eg. -lm
LDFLAGS = # For flags used when linking the .o files

VPATH = ..:../src:../external

target = MixMHC2pred

# Add/name different options based on OS (such that 'make' creates and 
# searches for the correct target name and possible use static libraries).
ifeq ($(OS),Windows_NT)  # is Windows_NT on XP, 2000, 7, Vista, 10...
    detected_OS = Windows
	target := $(target).exe
else
    detected_OS = $(shell sh -c 'uname 2>/dev/null || echo Unknown')
	ifeq ($(detected_OS),Linux)
		target := $(target)_unix
	endif
	ifeq ($(detected_OS),Unknown)
		target = $(error Unknown OS)
	endif
endif
ifneq ($(detected_OS),Darwin)
	# LDFLAGS += -static-libstdc++
	LDFLAGS += -static
	# To include the libraries as static to avoid some issues when sharing the
	# program (especially on windows; but it makes bigger files).
	# This didn't work on my mac compiler (it didn't find one of the needed
	# libraries), that's why I'm only doing it for the other OS.
endif


all: $(target)


# When calling "make debug" it'll compile MixMHC2pred based on "debug
# configuration".
debug: CPPFLAGS += -g
debug: $(target)


$(target): MixMHC2pred.o helper_MixMHC2pred.o helper_nnets.o alleles_def.o \
		whereami.o helper_general.o
	$(LINK.cc) -o ../$@ $^


helper_MixMHC2pred.o: helper_MixMHC2pred.hpp ## helper_general.hpp
# No need to indicate the .cpp when using an implicit compile rule.
helper_nnets.o: helper_nnets.hpp
alleles_def.o: alleles_def.hpp
helper_general.o: helper_general.hpp

MixMHC2pred.o: helper_MixMHC2pred.hpp helper_nnets.hpp alleles_def.hpp argparse/argparse.hpp
# To tell that MixMHC2pred.o is also depending on these headers to
# recompile the prog in case the files change.


whereami.o: whereami/src/whereami.c
	gcc -c -o $@ $^

clean:
	rm -f ../$(target) *.o
