CC=cc

SILENCE= -Wno-unused-variable
INCLUDES=
LIBRARIES= -L. -lm
#-Wmissing-field-initializers # I think this was for nvcc or something -pedantic
CFLAGS= -fPIC -std=gnu99 -Wall -Wextra -Wsign-conversion $(SILENCE) -g -O3
CLIBS=
CFLAGS+= $(CLIBS)
#CFLAGS += $(INCLUDES)
#-Wsign-conversion (some nvidia libs can make this a noisy warning, might be fixed now)
LDFLAGS= -shared

ifeq ($(CC),nvcc)
	CFLAGS := --compiler-options '$(CFLAGS) -Wno-c++11-extensions -Wno-c++11-long-long' -g -O3 --use_fast_math
	LDFLAGS := --compiler-options '-fPIC ' $(LDFLAGS)
	INCLUDES += -I/usr/local/cuda/include
	#NVOPT= --maxrregcount=32
	#NVCCFLAGS += $(NVOPT)
endif


.PHONY: clean all

all: test.x

%.o: %.c *.h
	$(CC) $(CFLAGS) $(INCLUDES) -c $< -o $@

%.o: %.cu *.h
	$(NVCC) -x cu -dc $(NVCCLDFLAGS) $(NVCCFLAGS) $(INCLUDES) -c $< -o $@

liborbit.so: orbit_main.o orbit_config.o  orbit_particles.o orbit_equilibrium.o \
             orbit_perturbation.o orbit_deposition.o \
             orbit_util.o inih/ini.o \
             cuda_helpers.o
	$(CC) $(CFLAGS) $(LDFLAGS) $^ -o $@

test.x: test.o liborbit.so
	$(CC) $(CFLAGS) $^ -o $@ $(LIBRARIES) -lorbit

clean:
	-rm -f inih/*.o
	-rm -f *.o
	-rm -f *.oo
	-rm -f *.so
	-rm -f *.x
	-rm -rf test.x.dSYM
