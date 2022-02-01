# minepoch

Minepoch is a mini-app version of the relativistic electromagnetic particle in cell (PIC) code
[EPOCH](https://github.com/Warwick-Plasma/epoch). Minepoch is designed to act both as a proxy-app
for performance testing as well as a test-bed for PIC algorithm development. Minepoch contains the
same core physics module as EPOCH, in which Maxwell's equations are solved on an Eulerian mesh,
coupled to Lagrangian particles, as well as extensions to allow drift-kinetic particles to be
treated.

## Compiling the code

### Make

The code can be compiled using make, and an appropriately set compiler variable, for example:

```
$ make COMPILER=gfortran
```

The requirements are minimal, only a Fortran compiler and MPI are required to compile and run
the code. The code is routinely run on both Linux and OS X, and tested with gfortran, but
other compilers (e.g. `make COMPILER=intel`) ought to work.

### CMake

Alternatively minepoch can be compiled using cmake. For example:

```
$ mkdir build && cd build
$ cmake ..
$ make
$ make install
```

By default minepoch will be installed to the `bin` subdirectory of the minepoch directory.

After installation the temporary `build` directory may be deleted.

### Compiling implicit PIC

The implicit solver requires an installation of Trilinos. The implicit development of
minEPOCH was carried out using [v12.6.3](https://github.com/trilinos/Trilinos/releases/tag/trilinos-release-12-6-3)
but other versions should be compatible. To compile using the Makefile system, the
`Makefile.export.Trilinos` file should be copied (or symlinked) to the mineEPOCH directory:

```
ln -s /trilinos/12.6.3/include/Makefile.export.Trilinos Makefile.export.Trilinos
```

The Trilinos support is then enabled at compile time for minEPOCH:

```
make COMPILER=gfortran DEF=-DTRILINOS
```

A Docker image containing all the required dependencies is available:

```
docker pull tomgoffrey/minepoch_deps:latest
```

## Running the code

Minepoch is controlled using a Fortran namelist, which is read from `Data/input.deck`. The file
`Data/example_input.deck` contains all the run-time configurable options and their default values.
The code contains a number of predefined set-ups, which can be selected by changing the `problem`
variable of the `CONTROL` namelist block.

A number of example input decks are included with the code:

 - `Data/two_stream.deck`: A two-stream instability, run with the standard PIC algorithm
 - `Data/two_stream_substep.deck`: The same two-stream instability, but with particle
    substepping enabled.
 - `Data/one_stream.deck`: A 1D force-free problem illustrating the use of low-noise control-variates.
 - `Data/drift_kinetic.deck`: An demonstration of the use of drift-kinetics in a quasi-1D setting.
 - `Data/em_wave.deck`: A simple electromagnetic wave set-up in a vacuum. Also demonstrates use
    of the field probe diagnostic.

To run one of these problems the input deck must be copied to the correct location before running
minepoch, for example:

```
$ cp Data/two_stream.deck Data/input.deck
$ ./bin/epoch3d
```

### Running Implicit minEPOCH

In order to enable the implicit algorithm in minEPOCH one should add/set the following line in the
`Data/input.deck`:

```
explicit_pic = F,
```

To control the use of the pseudocurrent correction to Gauss's law the following lines should be added/modified:

```
use_pseudo_current = T,
pseudo_current_fac = 1e-8,
```

The tolerance of the solver can be set/modified as follows:

```
linear_tolerance = 1e-1,
nonlinear_tolerance = 1e-2,
```

The CFL number can be adjusted by changing/adding the following line:

```
dt_multiplier = 1.0,
```

Finally, it is possible to provide more detailed solver output by setting:

```
verbose_solver = T,
```

## Analysing the results

As minepoch is primarily designed to test performance, only minimal diagnostics are included. These
are written to `Data/output.dat`. Some example scripts for analysing the output can be found in
the directory `minepoch_py`. For example, from the root minEPOCH directory:

```
$ python3 minepoch_py/energy.py --help

usage: energy.py [-h] [--energy_conservation [ENERGY_CONSERVATION]]
                 [--disable_plots]

Energy conservation analysis for minEPOCH

optional arguments:
  -h, --help            show this help message and exit
  --disable_plots       Disable plot generation (default: False)

Tolerances:
  --energy_conservation [ENERGY_CONSERVATION]
                        Allowed fractional error in energy conservation
                        (default: 0.0001)
```

one can obtain information in the use of the energy conservation analysis tool. Similar information is
available for `minepoch_py/test.py` and `minepoch_py/two_stream.py`.

The file `Data/output.dat` contains human-readable energy diagnostics analysis. Additionally the user may
request time-history of the electromagnetic fields at various points. This is demonstrated in the example
set-up, `Data/em_wave.deck`. Examining it's contents:

```
&CONTROL
  problem = 'em_wave',
  stdout_frequency = 25,
  n_field_probes = 2,
/

! List of positions to output field time-histories
! Two files will be produced:
! 1. ix=10, iy=2, iz=2
! 2. ix=180, iy=2, iz=2
&FIELD_PROBE_POSITIONS
  x_probes = 10, 180,
  y_probes = 2, 2,
  z_probes = 2, 2,
/
```

it can be seen that setting up field probes in a given cell is a two-step process. First the number of probes
must be set in the `CONTROL` namelist, here via the line `  n_field_probes = 2,`. Then the cell indices of the
probes are specified in the `FIELD_PROBE_POSITIONS`, where `x_probes` refers to the x-index of the probe, and
similar for `y_probes` and `z_probes`.

The field probes provide human readable data for the electric and magnetic field components as a function
of time, as well as some basic grid information to allow for convenient analysis.

Additionally the time history of particle momentum and Poynting vector summed over the computational domain
maybe output to a file `Data/momentum.dat` by setting

```
write_momentum = T,
```

in the `CONTROL` block of the `input.deck` file.

## Further Details

 - Information regarding the low noise PIC algorithm can be found in ExCALIBUR-NEPTUNE report 2047355-TN-03
   (in Docs directory of this repository).
 - Documentation on the advanced particle tracing routines (substepping and drift-kinetics) is included
   with the code. The documents (docs/report_implemented.pdf, docs/implicit_pic.pdf) can be generated by
   running `make docs` from the minEPOCH directory.
 - Although minepoch only includes a small subset of the functionality of the original EPOCH code, the
   user manual for EPOCH is included within each [release](https://github.com/ExCALIBUR-NEPTUNE/minepoch/releases).
 - The EPOCH code [paper](http://dx.doi.org/10.1088/0741-3335/57/11/113001) contains a review of the
   core PIC algorithm used in minepoch and EPOCH.
 - The original EPOCH code may be downloaded from the [repository](https://github.com/Warwick-Plasma/epoch).

## Any problems?

Please open an issue on the [issue tracker](https://github.com/ExCALIBUR-NEPTUNE/minepoch/issues).
