CC=cc

SILENCE = -Wno-unused-variable
INCLUDES =
LIBRARIES = -L. -lm
CFLAGS= -fPIC $(SILENCE) -Wall -Wextra -Wsign-conversion -g -O3
CFLAGS += $(INCLUDES)
LDFLAGS= -Wl,-rpath,.

# if we are compiling towards CUDA or Host, we alter our flags accordingly
ifeq ($(CC),nvcc)
	# treat files as CU, treat files as device-code-relocatable (dc)
	NVCCFLAGS += -x=cu -dc
	CFLAGS := $(NVCCFLAGS) --compiler-options '$(CFLAGS) ' -g -O3 --use_fast_math
	LDFLAGS := --compiler-options '-fPIC ' -Xcompiler \"$(LDFLAGS)\"
	INCLUDES += -I/usr/local/cuda/include
else
	CFLAGS += -std=gnu99 -pedantic
	LDFLAGS += $(CFLAGS)
endif


all: unpack cuOrbit.x

# alternatively you can yourself just invoke: make CC=nvcc
gpu:
	$(MAKE) all CC=nvcc

%.o: %.c *.h
	$(CC) $(CFLAGS) $(INCLUDES) -c $< -o $@

## We don't currently need, but maybe one day if you add CUDA only optimizations
%.o: %.cu *.h
	$(NVCC) -x cu -dc $(NVCCLDFLAGS) $(NVCCFLAGS) $(INCLUDES) -c $< -o $@

libcuorbit.so: inih/ini.o \
		orbit_config.o \
		orbit_deposition.o \
		orbit_equilibrium.o \
		orbit_main.o \
		orbit_particles.o \
		orbit_perturbation.o \
		orbit_util.o \
		cuda_helpers.o
	$(CC) -shared $(LDFLAGS) $^ -o $@

cuOrbit.x: cuOrbit.o libcuorbit.so
	$(CC) $(LDFLAGS) $< -o $@ -lcuorbit $(LIBRARIES)

GZREFERENCES = reference_test_case/displ_4_5.dat \
		reference_test_case/pDEDP.AEP \
		reference_test_case/spdata \
		INPUT/spdata \
		INPUT/pDEDP.AEP

unpack: $(GZREFERENCES)

%: %.gz
	gunzip -c $< > $@

clean:
	-rm -f inih/*.o
	-rm -f *.o
	-rm -f libcuorbit.so
	-rm -f cuOrbit.x
	-rm -rf cuOrbit.x.dSYM
	-rm -f $(GZREFERENCES)

.PHONY: all
.PHONY: clean
.PHONY: unpack
