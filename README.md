# minepoch

Minepoch is a mini-app version of the relativistic electromagnetic particle in cell (PIC) code
[EPOCH](https://github.com/Warwick-Plasma/epoch). Minepoch is deisgned to act both as a proxy-app
for performance testing as well as a test-bed for PIC algorithm development. Minepoch contains the
same core physics module as EPOCH, in which Maxwell's equations are solved on an Eulerian mesh,
coupled to Lagrangian particles.

## Compiling the code

The code can be compiled using make, and an appropriately set compiler variable, for example:

```
$ make COMPILER=gfortran
```

The requirements are minimal, only a Fortran compiler and MPI are required to compile and run
the code. The code is routinely run on both Linux and OS X, and tested with gfortran, but
other compilers (e.g. `make COMPILER=intel`) ought to work.

## Running the code

Minepoch is controlled using a Fortran namelist, which is read from `Data/input.deck`. The file
`Data/example_input.deck` contains all the run-time configurable options and their default values.
The code contains a number of predefined set-ups, which can be selected by changing the `problem`
variable of the `CONTROL` namelist block.

A number of example input decks are included with the code:

 - `Data/two_stream.deck`: A two-stream instability, run with the standard PIC algorithm
 - `Data/two_stream_substep.deck`: The same two-stream instability, but with particle
    substepping enabled.
 - `Data/one_stream.deck`: ???
 - `example_decks/one_stream.deck`: ???
 
To run one of these problems the input deck must be copied to the correct location before running
minepoch, for example:

```
$ cp Data/two_stream.deck Data/input.deck
$ ./bin/epoch3d
```

## Analysing the results

As minepoch is primarily designed to test performance, only minimal diagnostics are included. These
are written to `Data/output.dat`. Some example scripts for analysing the output can be found in
the directory `minepoch_py`.

## Any problems?

Please open an issue on the [issue tracker](https://github.com/ExCALIBUR-NEPTUNE/minepoch/issues).
