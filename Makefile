CC=cc
NVCC=nvcc

CFLAGS=' -Wall -Wextra -Wno-missing-field-initializers -Wsign-conversion -g -O3 '
CLIBS=
CFLAGS+= $(CLIBS)
#CFLAGS += $(INCLUDES)
#-Wsign-conversion (some nvidia libs can make this a noisy warning, might be fixed now)
#LDFLAGS= --shared
NVCCFLAGS= --compiler-options $(CFLAGS) -g -O3 --use_fast_math
NVCCLDFLAGS= --compiler-options '-shared -fPIC '
#NVOPT= --maxrregcount=32
#NVCCFLAGS += $(NVOPT)


.PHONY: clean all

all: orbit_structures.o

# libcuorbit.so:
# 	$(NVCC) $(NVCCLDFLAGS) $(NVCCFLAGS) $^ -o $@

%.o: %.c
	$(NVCC) -x cu -dc $(NVCCLDFLAGS) $(NVCCFLAGS) -DSKIPMAIN -c $< -o $@

%.o: %.cu
	$(NVCC) -x cu -dc $(NVCCLDFLAGS) $(NVCCFLAGS) -DSKIPMAIN -c $< -o $@


clean:
	-rm -f *.o
	-rm -f *.oo
	-rm -f *.so
	-rm -f *.x
