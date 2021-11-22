# polyomino-evolution
Simulation of evolution and complexity in polymino assembly model

Requirements: C, R with `ggplot2` , Bash for the wrapper script

`run-blend.sh` is wrapper code, which compiles and runs a set of experiments in Bash

`sim-poly-blend.c` wraps an individual experiment. It takes command-line arguments (with defaults):
* `--directed 1` (directed evolution? 0 no, 1 size-based fitness, -1 random fitness)
* `--npar 10` (number of individuals in evolutionary population)
* `--mut 0.1000` (mutation rate per genome bit)
* `--targetsize 16` (target size for size-based fitness function)
* `--ntile 16` (number of tiles in an assembly genome)
* `--ncol 64` (number of colours in an assembly genome)
* `--numr 5000` (number of evolutionary runs)
* `--nsamp 100000000` (number of samples of genotype space)
* `--outputall 0` (output all structures found, or just top-fitness ones? 0 just top fitness, 1 all)
* `--confound 0` (apply a confounding step to genomes? 0 no, 1 cyclic permutation based on 1s count, >= 2 random changes with prob (argument-1)/(genome length))
* `--rseed 1` (random seed)

The outputs, depending on the experiment, are files labelled as follows:
* `out-` -- general output, including structures discovered, bitstring and reduced genomes, and evolutionary timescales
* `stats-` -- (CSV) statistics of recorded structures, including size, symmetry, modularity, various complexity measures, frequencies in evolution and sampling
* `lib-` -- shapes of recorded structures, as ASCII characters on a square grid (`.` denotes an empty site)
* `pop-` -- population structure over time and simulation instances, with references referring to the structures in `stats-` and `lib-` (`-1` denotes an unbound or non-deterministic individual)

The other `.c` files take care of things under the hood:
* `assembly.c` -- polyomino assembly
* `ga.c` -- genetic algorithm
* `library.c` -- structure comparison and storage
* `stats-variable.c` -- various summary and complexity statistics for polyomino structures

Visualisation is done in R:
* `plot-blend.R` plots various summaries of the outputs, corresponding roughly to manuscript figures
* `plot-structs.R` , passed a command-line reference to a particular experiment, draws the top polyomino structures from that experiment
* `plot-complexity-hists.R` plots complexity histograms over (simulation) time
* (files with prepended `nz-` produce these plots using the `nnonzero` (number of nonzero genome sites) complexity measure instead of the `ninterface` (number of interface types)

