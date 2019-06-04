CC=cc
NVCC=nvcc

SILENCE= -Wno-unused-variable
#-Wmissing-field-initializers # I think this was for nvcc or something
CFLAGS= -Wall -Wextra -Wsign-conversion -pedantic $(SILENCE) -g -O3
CLIBS=
CFLAGS+= $(CLIBS)
#CFLAGS += $(INCLUDES)
#-Wsign-conversion (some nvidia libs can make this a noisy warning, might be fixed now)
LDFLAGS= -shared
NVCCFLAGS= --compiler-options '$(CFLAGS)' -g -O3 --use_fast_math
NVCCLDFLAGS= --compiler-options '-shared -fPIC '
#NVOPT= --maxrregcount=32
#NVCCFLAGS += $(NVOPT)


.PHONY: clean all

all: test.x

# libcuorbit.so:
# 	$(NVCC) $(NVCCLDFLAGS) $(NVCCFLAGS) $^ -o $@

%.o: %.c %.h
	$(CC) $(CFLAGS) -c $< -o $@

%.o: %.cu %h
	$(NVCC) -x cu -dc $(NVCCLDFLAGS) $(NVCCFLAGS) -c $< -o $@

liborbit.so: orbit_config.o  orbit_particles.o orbit_equilibrium.o \
             orbit_perturbation.o orbit_util.o inih/ini.o
	$(CC) $(CFLAGS) $(LDFLAGS) $^ -o $@

test.x: test.o liborbit.so
	$(CC) $(CFLAGS) $^ -o $@ -L. -lorbit

clean:
	-rm -f inih/*.o
	-rm -f *.o
	-rm -f *.oo
	-rm -f *.so
	-rm -f *.x
	-rm -rf test.x.dSYM
