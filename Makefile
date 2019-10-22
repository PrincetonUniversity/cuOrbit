CC=cc

SILENCE= -Wno-unused-variable
INCLUDES=
LIBRARIES= -L. -lm
#-Wmissing-field-initializers # I think this was for nvcc or something -pedantic
CFLAGS= -fPIC -std=gnu99 $(SILENCE) -Wall -Wextra -Wsign-conversion -g -O3 
#CFLAGS += $(INCLUDES)
#-Wsign-conversion (some nvidia libs can make this a noisy warning, might be fixed now)
LDFLAGS=

ifeq ($(CC),nvcc)
	# treat files as CU, treat files as device-code-relocatable (dc)
	NVCCFLAGS += -x=cu -dc
	CFLAGS := $(NVCCFLAGS) --compiler-options '$(CFLAGS) ' -g -O3 --use_fast_math
	LDFLAGS := --compiler-options '-fPIC ' $(LDFLAGS)
	INCLUDES += -I/usr/local/cuda/include
	#NVOPT= --maxrregcount=32

else
	CFLAGS += -pedantic -Wno-c++11-extensions -Wno-c++11-long-long
	LDFLAGS +=$(CFLAGS)
endif


.PHONY: clean all

all: test.x

%.o: %.c *.h
	$(CC) $(CFLAGS) $(INCLUDES) -c $< -o $@

# %.o: %.cu *.h
# 	$(NVCC) -x cu -dc $(NVCCLDFLAGS) $(NVCCFLAGS) $(INCLUDES) -c $< -o $@

liborbit.so: orbit_main.o orbit_config.o  orbit_particles.o orbit_equilibrium.o \
             orbit_perturbation.o orbit_deposition.o \
             orbit_util.o inih/ini.o \
             cuda_helpers.o
	$(CC) -shared $(LDFLAGS) $^ -o $@

test.x: test.o liborbit.so
	$(CC) $(LDFLAGS) $< -o $@ $(LIBRARIES) -lorbit

clean:
	-rm -f inih/*.o
	-rm -f *.o
	-rm -f *.oo
	-rm -f *.so
	-rm -f *.x
	-rm -rf test.x.dSYM
