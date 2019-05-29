CC=cc
NVCC=nvcc

CFLAGS= -Wall -Wextra -Wno-missing-field-initializers -Wsign-conversion -g -O3
CLIBS=
CFLAGS+= $(CLIBS)
#CFLAGS += $(INCLUDES)
#-Wsign-conversion (some nvidia libs can make this a noisy warning, might be fixed now)
#LDFLAGS= --shared
NVCCFLAGS= --compiler-options '$(CFLAGS)' -g -O3 --use_fast_math
NVCCLDFLAGS= --compiler-options '-shared -fPIC '
#NVOPT= --maxrregcount=32
#NVCCFLAGS += $(NVOPT)


.PHONY: clean all

all: test.x

# libcuorbit.so:
# 	$(NVCC) $(NVCCLDFLAGS) $(NVCCFLAGS) $^ -o $@

%.o: %.c
	$(CC) $(CFLAGS) -c $< -o $@

%.o: %.cu
	$(NVCC) -x cu -dc $(NVCCLDFLAGS) $(NVCCFLAGS) -c $< -o $@

test.x: test.o orbit_config.o  orbit_particles.o orbit_equilibrium.o orbit_perturbation.o \
orbit_util.o
	$(CC) $(CFLAGS) $^ -o $@

clean:
	-rm -f *.o
	-rm -f *.oo
	-rm -f *.so
	-rm -f *.x
