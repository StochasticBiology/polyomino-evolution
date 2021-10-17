# polyomino-evolution
Simulation of evolution and complexity in polymino assembly model

Requirements: C, R with `ggplot2` 

`sim-poly-blend.c` is the workhorse code. It takes command-line arguments (with defaults):
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

`plot-blend.R` plots various summaries of the outputs

`plot-structs.R` , passed a command-line reference to a particular experiment, draws the top polyomino structures from that experiment

`run-blend.sh` compiles and runs a set of experiments in Bash
