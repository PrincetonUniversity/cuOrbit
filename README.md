# cuOrbit

CUDA Adaptation of an Orbit based KICK Model.

This code is currently in active debugging towards validation.

*It is demonstrably over an order of magnitude faster.*

## About

This is a derived work of R White's ORBIT code (FORTRAN IV/66/77).
Specifically this was based on a reductionist OpenMP implementation
by Mario Podesta, where the code was focused towards study of the KICK operator.
See:

https://iopscience.iop.org/article/10.1088/1361-6587/aa7977

https://doi.org/10.1063/1.864527

I think the bulk of computations are probably present, while bulk of diagnostics
and extended abilities of the ORBIT code are intentionally omitted at this time.

## Getting Started

### Requirements

You will need `gcc`. The code has been built with a few versions,
and is _mostly_ ANSI compliant, though I admittedly use GNU extensions
(ie, -std=gnu99).

If you want to use the CUDA GPU acceleration, you will need `nvcc`
and a capable CUDA GPU.

_NVCC and a GPU is not at all required to build the host version._

If your systems use modules, you can probably load what you need that way...
For example, at PPPL, you could probably do something like the following:

```
module load gcc/8.1.0
module load cuda/10.1.243
```

### Building

For the host version:

```
make
```

For the gpu version:

```
make gpu
```

The makefile is pretty simple, and a lot of commonly edited attributes are
already abstracted to variables near top of file.

### Running

Either of the two build methods will create a library `libcuorbit.so`,
and an executable `cuOrbit.x`.
The executable runs using only that library with an associated `orbit_config_api.h` header.

```
./cuOrbit.x
```

On some machines the local directory is not searched for libraries,
in which case you can simply `export LD_LIBRARY_PATH=.:$LD_LIBRARY_PATH`
or equivalent.

### Configuration

You can use the configuration file `INPUT/config.ini` to change many parameters
without recompiling any code. It is a simple INI format, `Key = Value`.

Input data is suggestively located in `INPUT/` along with the config.ini.
During inititialization the input data is read from there into the Config_t.
You can alter input file names in the configration file to point at different
input decks, and also configure outputs to a different path as you wish.

I have a note to add a feature that will optionally take a config by filename
on the commandline for more convenience.

### Using in other programs

Peeking inside the `cuOrbit.c` program, we instantiate a `orbin_Config_t`,
`orbit_initialize_Config` it, and invoke the `orbit_main_loop`.

The idea here is that a calling program would do something similar,
potentially altering or packing the `Config_t`,
instead of loading exclusively from the ini file.
That is, in a program, you may choose to override data by changing the values
in the configuration types (C structures).
Some helper code coud be devised for that application,
particularly if you are in a different language (F, python).
As that use case is explored, supporting methods can be implemented.

## License

GPLv3 is in discussion with R White.

## Financial Support

This project was supported by TRANSP development at the
Princeton Plasma Physics Laboratory via U.S. Department of Energy
(DE-AC02-09CH11466).

Additionally, the code was developed and debugged utilizing a
small donation of GPU compute resources by Garrett Wright.
