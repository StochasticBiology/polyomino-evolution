#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#include "assembly.c"
#include "ga.c"
#include "stats-variable.c"
#include "library.c"

#define RND drand48()

// simple fitness function (if needed): f = 1/(|s-s*|+1)
// takes the grid representation of a polyomino as input and computes size
float Fitness(int *g, int ARR, int TARGETSIZE)
{
  int i;
  float f = 0;

  for(i = 0; i < ARR*ARR; i++)
    f += (g[i] != -1);

  return 1.0/(fabs(f - TARGETSIZE) + 1);
}

int main(int argc, char *argv[])
{
  int *P;
  int *confoundP;
  float *f;
  int *g, *tmp;
  int i, j, t;
  int size;
  int *Q;
  double *meandtime, *mindtime, *sampledcount;
  int *sizelib;
  int *lib;
  int ref;
  int numlib;
  int newpheno;
  int adaptcount;
  int discovered, adapted;
  int *discoveredpheno;
  int run;
  char string[200];
  char outstr[200];
  int *popnsnap;
  long int *adaptoccur, *occur;
  int *symms;
  float *tmpfitness;
  int modularity, *bestmodularity;
  FILE *output;
  FILE *fp;
  FILE *popfp;
  char popstr[200];
  int symmfreq[6];
  int RSEED;
  int numrecorded;
  ComplexityMeasures *CM;
  ComplexityMeasures tmpCM;

  // global variables (somewhat awkwardly) store parameters for simulation
  int DIRECTED, NPAR, NTILE, NCOL, NBITCOL, TARGETSIZE, LEN;
  long NUMR;
  long NSAMP;
  double MUT;
  int OUTPUTALL;
  int CONFOUND;

  // relatively fixed settings
  int ARR = 16;              // size of assembly grid
  int MAXLIB = 10000;        // max number of structures that can be stored in library
  int MAXT = 20000;          // max simulation timesteps
  int MAXBLOCKS = 100;       // max number of blocks used by a structure (used in modularity calculation)

 
  // default parameterisation
  DIRECTED = 1; NPAR = 10; MUT = 0.1; TARGETSIZE = 16; NTILE = 16; NCOL = 64; NBITCOL = 6; NUMR = 5000; NSAMP = 1e8;
  OUTPUTALL = 0; CONFOUND = 0;
  RSEED = 1;

  // process command-line arguments to alter these parameters
  for(i = 1; i < argc; i++)
    {
      if(strcmp(argv[i], "--help") == 0)
	{
	  printf("Arguments (with defaults):\n  --directed %i\n  --npar %i\n  --mut %.4f\n  --targetsize %i\n  --ntile %i\n  --ncol %i\n  --numr %li\n  --nsamp %li\n  --outputall %i\n  --confound %i\n  --rseed %i\n\n", DIRECTED, NPAR, MUT, TARGETSIZE, NTILE, NCOL, NUMR, NSAMP, OUTPUTALL, CONFOUND, RSEED);
	  exit(0);
	}
      if(strcmp(argv[i], "--directed") == 0 && i+1 < argc) { DIRECTED = atoi(argv[i+1]); i++; }
      if(strcmp(argv[i], "--npar") == 0 && i+1 < argc) { NPAR = atoi(argv[i+1]); i++; }
      if(strcmp(argv[i], "--mut") == 0 && i+1 < argc) { MUT = atof(argv[i+1]); i++; }
      if(strcmp(argv[i], "--targetsize") == 0 && i+1 < argc) { TARGETSIZE = atoi(argv[i+1]); i++; }
      if(strcmp(argv[i], "--ntile") == 0 && i+1 < argc) { NTILE = atoi(argv[i+1]); i++; }
      if(strcmp(argv[i], "--ncol") == 0 && i+1 < argc) { NCOL = atoi(argv[i+1]); i++; }
      if(strcmp(argv[i], "--numr") == 0 && i+1 < argc) { NUMR = atol(argv[i+1]); i++; }
      if(strcmp(argv[i], "--nsamp") == 0 && i+1 < argc) { NSAMP = atol(argv[i+1]); i++; }
      if(strcmp(argv[i], "--outputall") == 0 && i+1 < argc) { OUTPUTALL = atoi(argv[i+1]); i++; }
      if(strcmp(argv[i], "--confound") == 0 && i+1 < argc) { CONFOUND = atoi(argv[i+1]); i++; }
      if(strcmp(argv[i], "--rseed") == 0 && i+1 < argc) { RSEED = atoi(argv[i+1]); i++; }
    }

  switch(NCOL)
    {
    case 2: NBITCOL = 1; break;
    case 4: NBITCOL = 2; break;
    case 8: NBITCOL = 3; break;
    case 16: NBITCOL = 4; break;
    case 32: NBITCOL = 5; break;
    case 64: NBITCOL = 6; break;
    case 128: NBITCOL = 7; break;
    default: printf("Unsupported number of colours!\n"); exit(1); break;
    }
  
  printf("  --directed %i\n  --npar %i\n  --mut %.4f\n  --targetsize %i\n  --ntile %i\n  --ncol %i\n  --nbitcol %i\n  --numr %li\n  --nsamp %li\n  --outputall %i\n  --confound %i\n  --rseed %i\n\n", DIRECTED, NPAR, MUT, TARGETSIZE, NTILE, NCOL, NBITCOL, NUMR, NSAMP, OUTPUTALL, CONFOUND, RSEED);
  
  sprintf(outstr, "out-blend-%i-%.3f-%i-%i-%i-%i-%i-%i-%i-%li-%.0e.csv", DIRECTED, MUT, NPAR, TARGETSIZE, NTILE, NCOL, OUTPUTALL, CONFOUND, RSEED, NUMR, (double)NSAMP);
  sprintf(popstr, "pop-%s", outstr);
  
  srand48(RSEED);
  
  /*      switch(cmdexpt)
	  {
	  case 0: DIRECTED = 0; NUMR = 500; break;
	  case 1: DIRECTED = -1; NUMR = 500; break;
	  case 2: NPAR = 100; break;
	  case 3: MUT = 1; break;
	  case 4: TARGETSIZE = 8; break;
	  case 5: TARGETSIZE = 15; break;
	  case 6: MUT = 0.01; break;
	  case 7: MUT = 0.001; break;
	  case 8: NUMR = 0; NSAMP = NSAMP / 10; break;
	  case 9: NUMR = 0; NSAMP = NSAMP / 100; break;
	  } */
  
  LEN = (4*NTILE*NBITCOL);

  //// allocate memory for storing the various statistics of polyomino structures
  // different complexity measures
  CM = (ComplexityMeasures*)malloc(sizeof(ComplexityMeasures)*MAXLIB);
  // discovery time stats
  meandtime = (double*)malloc(sizeof(double)*MAXLIB);
  mindtime = (double*)malloc(sizeof(double)*MAXLIB);
  // occurrence and discovery counts
  sampledcount = (double*)malloc(sizeof(double)*MAXLIB);
  discoveredpheno = (int*)malloc(sizeof(int)*MAXLIB);
  adaptoccur = (long int*)malloc(sizeof(long int)*MAXLIB);
  occur = (long int*)malloc(sizeof(long int)*MAXLIB);
  // structural statistics
  symms = (int*)malloc(sizeof(int)*MAXLIB);
  bestmodularity = (int*)malloc(sizeof(int)*MAXLIB);

  tmpfitness = (float*)malloc(sizeof(float)*MAXLIB);

  //// allocate memory for the evolutionary simulation
  // population of genomes
  P = (int*)malloc(sizeof(int)*NPAR*LEN);
  confoundP = (int*)malloc(sizeof(int)*LEN);
  // fitnesses
  f = (float*)malloc(sizeof(float)*NPAR);
  // population snapshot
  popnsnap = (int*)malloc(sizeof(int)*NPAR);
  // phenotypes, stored as an ARR*ARR grid
  g = (int*)malloc(sizeof(int)*ARR*ARR);
  tmp = (int*)malloc(sizeof(int)*ARR*ARR);
  // list of coloured edges
  Q = (int*)malloc(sizeof(int)*NTILE*4);

  //// now the big one, the library that will store MAXLIB ARR*ARR grids of phenotypes
  lib = (int*)malloc(sizeof(int)*(MAXLIB+5)*ARR*ARR);
  // and the size of each (used for stats and for quick(er) searching)
  sizelib = (int*)malloc(sizeof(int)*MAXLIB);

  // initialise output file
  output = fopen(outstr, "w");
  fclose(output);

  // initialise statistics
  numlib = 0;
  for(i = 0; i < MAXLIB; i++)
    {
      CM[i].nblock = CM[i].ninterface = CM[i].nnonzero = CM[i].necklace = CM[i].lzc = 10000;
      bestmodularity[i] = -1;
      occur[i] = adaptoccur[i] = 0;
      mindtime[i] = MAXT;
      discoveredpheno[i] = meandtime[i] = 0;
    }

  //// GENETIC ALGORITHM PHASE
  // main loop, while we still have simulations to do and we haven't maxed out our phenotype library space
  if(OUTPUTALL)
    popfp = fopen(popstr, "w");
  for(run = 0; run < NUMR && numlib < MAXLIB; run++)
    {
      // initialise fitness values and simulation properties
      for(i = 0; i < MAXLIB; i++)
	tmpfitness[i] = -1;
      discovered = adapted = 0;
      // initialise population
      for(i = 0; i < NPAR*LEN; i++)
	P[i] = 0;

      output = fopen(outstr, "a");
 
      // a single evolutionary simulation
      // loop while we still have evolutionary time, haven't adapted a top-fitness organism (if appropriate), and still have library space
      for(t = 0; t < MAXT+1 && (DIRECTED != 1 || adapted == 0) && numlib < MAXLIB; t++)
	{
	  // periodic output tracker
	  if(t % 10000 == 0)
	    {
	      fclose(output);
	      output = fopen(outstr, "a");
	      fprintf(output, "%i %i\n", run, t);
	      printf("%i %i\n", run, t);
	    }
	  if(OUTPUTALL && t % 100 == 0)
	    fprintf(popfp, "%i,%i,", run, t);
	  
	  // initialise count of fixed top-fitness phenotypes
	  adaptcount = 0;

	  // loop through population
	  for(i = 0; i < NPAR; i++)
	    {
	      // produce phenotype for this individual
	      ref = -1;
	      newpheno = 0;
	      // if we're confounding genotypes, do so, and convert the result to a edge list for assembly tiles
	      if(CONFOUND)
		{
		  Confound(&(P[i*LEN]), confoundP, CONFOUND, LEN);
		  Convert(confoundP, Q, NTILE, NBITCOL);
		}
	      else
		{
		  Convert(&(P[i*LEN]), Q, NTILE, NBITCOL);
		}
	      // the Grow function grows a polyomino in g from edge list Q and returns nonzero if a UND structure is produced
	      if(Grow(Q, NTILE, g, &size, 0, 0, 20, 0, ARR) != 0)
		{
		  ref = -1;
		  f[i] = 0;
		}
	      else 
		{
		  // we've got a nominally stable phenotype, so assign a fitness (size-based, equal (undirected), or random)
		  switch(DIRECTED)
		    {
		    case 1: f[i] = Fitness(g, ARR, TARGETSIZE);  break;
		    case 0: f[i] = 1; break;
		    case -1: 
		      if(tmpfitness[ref] == -1)
			f[i] = tmpfitness[ref] = RND;
		      else 
			f[i] = tmpfitness[ref];
		      break;
		    }
		
		  // decide whether to add this structure to the library
		  // if we want to store all the structures we find, or if evolution is undirected, do so
		  // otherwise, add it if it's the first top-fitness target structure to be found in this simulation
		  if(OUTPUTALL == 1 || (DIRECTED == 1 && f[i] == 1))
		    {
		      // SymmLibrary returns a reference to an existing structure or adds it and returns the new reference
		      ref = SymmLibrary(g, size, lib, sizelib, &numlib, ARR, &newpheno);
		      if(numlib > MAXLIB)
			{
			  printf("Hit maximum structure capacity!\n");
			  break;
			}
		      // compute number of blocks in this phenotype -- if this gives a better modularity than we had previously, update it
		      modularity = Modularity(g, ARR, MAXBLOCKS);
		      if(bestmodularity[ref] == -1 || modularity < bestmodularity[ref])
			bestmodularity[ref] = modularity;

		      // remember that this structure came up
		      popnsnap[i] = ref;
		      occur[ref]++;
		      if(newpheno != 0)
			{
			  //// we haven't seen this structure before, so:
			  // output it to file
			  fprintf(output, "Library %i\n", ref);

			  for(j = 0; j < LEN; j++)
			    fprintf(output, "%i", P[i*LEN+j]);
			  fprintf(output, "\n");
			  for(j = 0; j < 4*NTILE; j++)
			    fprintf(output, "%i", Q[j]);
			  fprintf(output, "\n");
			  FileOutputGrid(output, g, ARR);

			  // get its complexity and store this
			  GenomeComplexityMultiStats(Q, NTILE, 0, &tmpCM, LEN, NBITCOL, output);
			  //		      fprintf(output, "%.5f %.5f %.5f\n", ic, sc, ic2);
			  if(tmpCM.nblock < CM[ref].nblock) CM[ref].nblock = tmpCM.nblock;
			  if(tmpCM.ninterface < CM[ref].ninterface) CM[ref].ninterface = tmpCM.ninterface;
			  if(tmpCM.nnonzero < CM[ref].nnonzero) CM[ref].nnonzero = tmpCM.nnonzero;
			  if(tmpCM.necklace < CM[ref].necklace) CM[ref].necklace = tmpCM.necklace;
			  if(tmpCM.lzc < CM[ref].lzc) CM[ref].lzc = tmpCM.lzc;
		  
			  // store phenotype in a canonical form (translated in grid)
			  Canonical(&(lib[ref*ARR*ARR]), tmp, ARR);

			  // identify symmetry group
			  switch(SymmCount(tmp, ARR))
			    {
			    case 0: fprintf(output, "D4\n"); break;
			    case 1: fprintf(output, "C4\n"); break;
			    case 2: fprintf(output, "D2\n"); break;
			    case 3: fprintf(output, "C2\n"); break;
			    case 4: fprintf(output, "D1\n"); break;
			    case 5: fprintf(output, "C1\n"); break;
			    }
			  fprintf(output, " %.5f\n", f[i]);
			}

		      // record that, and when, we discovered this structure
		      discoveredpheno[ref]++;
		      meandtime[ref] += t;
		      if(t < mindtime[ref])
			mindtime[ref] = t;

		      if(DIRECTED == 1 && f[i] == 1)
			{
			  if(discovered == 0)
			    {
			      // this is the first top-fitness structure discovered in this simulation
			      discovered = 1;
			      fprintf(output, "Discovered at %i\n", t);
			      printf("Discovered at %i\n", t);
			      FileOutputGrid(output, g, ARR);
			    }
			  // increment the population count of top-fitness structures
			  adaptcount++;

			  // if we're now over 50% adapted but weren't before 
			  if(adaptcount > NPAR/2 && adapted == 0)
			    {
			      // record that this structure was fixed
			      adaptoccur[ref]++;
			      adapted = 1;
			      fprintf(output, "Adapted at %i\n", t);
			      printf("Adapted at %i\n", t);
			      FileOutputGrid(output, g, ARR);
			    }
			}
		    }
		}
	      if(OUTPUTALL && t % 100 == 0)
		fprintf(popfp, "%i%c", ref, (i == NPAR-1 ? '\n' : ','));
	    }

	  // this is the genetic algorithm, breeding and mutating genomes
	  BreedGen(P, f, NPAR, LEN, ((float)MUT)/LEN);
	}

      fclose(output);
    }

  if(OUTPUTALL)
    fclose(popfp);
  
  //// SAMPLING PHASE
  // initialise sampling counts
  for(i = 0; i < MAXLIB; i++)
    sampledcount[i] = 0;

  // if we're using this to refine stats on simulated structures, prevent adding new structures to library
  // alternative: no simulation, just sampling
  if(NUMR != 0)
    numrecorded = numlib;

  // big sampling loop
  for(t = 0; t < NSAMP; t++)
    {
      // periodic output
      if(t % (int) 1e5 == 0)
	printf("Sampling %li of %li\n", (long int) t, NSAMP);
      // generate random edge list
      for(i = 0; i < LEN; i++)
	P[i] = RND < 0.5;
      if(CONFOUND)
	{
	  Confound(P, confoundP, CONFOUND, LEN);
	  Convert(confoundP, Q, NTILE, NBITCOL);
	}
      else
	{
	  Convert(P, Q, NTILE, NBITCOL);
	}

      i = 0; newpheno = 0;
      // build polyomino
      if(Grow(Q, NTILE, g, &size, 0, 0, 20, 0, ARR) != 0)
	f[i] = 0;
      else 
	{
	  // compute statistics of this polyomino grown using this genome
	  // if any statistics are "better" (e.g. lower complexity) than our previously recorded values, replace them
	  f[i] = Fitness(g, ARR, TARGETSIZE);
	  ref = SymmLibrary(g, size, lib, sizelib, &numlib, ARR, &newpheno); 

	  modularity = Modularity(g, ARR, MAXBLOCKS);
	  if(bestmodularity[ref] == -1 || modularity < bestmodularity[ref])
	    bestmodularity[ref] = modularity;

	  GenomeComplexityMultiStats(Q, NTILE, 0, &tmpCM, LEN, NBITCOL, output);
	  //		      fprintf(output, "%.5f %.5f %.5f\n", ic, sc, ic2);
	  if(tmpCM.nblock < CM[ref].nblock) CM[ref].nblock = tmpCM.nblock;
	  if(tmpCM.ninterface < CM[ref].ninterface) CM[ref].ninterface = tmpCM.ninterface;
	  if(tmpCM.nnonzero < CM[ref].nnonzero) CM[ref].nnonzero = tmpCM.nnonzero;
	  if(tmpCM.necklace < CM[ref].necklace) CM[ref].necklace = tmpCM.necklace;
	  if(tmpCM.lzc < CM[ref].lzc) CM[ref].lzc = tmpCM.lzc;

	  sampledcount[ref]++;
	  // recap library count if we're refining simulated stats, or if we've exceeded our library size
	  if(NUMR != 0)
	    numlib = numrecorded;
	  if(numlib > MAXLIB-1)
	    numlib = MAXLIB-1;
	}
    }

  //// OUTPUT PHASE

  // count symmetries discovered in this ensemble
  // initialise symmetry count
  for(i = 0; i < 6; i++)
    symmfreq[i] = 0;
  // loop through library of structures, storing symmetries
  for(i = 0; i < numlib; i++)
    {
      if(occur[i] || NUMR == 0)
	{
	  Canonical(&(lib[i*ARR*ARR]), tmp, ARR);
	  symms[i] = SymmCount(tmp, ARR);
	  symmfreq[symms[i]] += discoveredpheno[i];
	}
    }
  // output symmetry stats
  sprintf(string, "symm-%s", outstr);
  fp = fopen(string, "w");
  fprintf(fp, "SymmGroup,Count\n");
  for(i = 0; i < 6; i++)
    fprintf(fp, "%i,%i\n", i, symmfreq[i]);
  fclose(fp);

  // record stats for this ensemble
  sprintf(string, "stats-%s", outstr);
  fp = fopen(string, "w");
  fprintf(fp, "Label,DiscoveryCount,Size,CMNBlock,CMNInterface,CMNNonzero,CMLZW,CMNecklace,MeanDiscoveryTime,MinDiscoveryTime,SampleCount,Occurrence,AdaptCount,Symmetry,Modularity\n");
  for(i = 0; i < numlib; i++)
    {
      if(NUMR == 0)
	{
	  fprintf(fp, "%i,%i,%i,%i,%i,%i,%i,%.5f,%.5f,%.5f,%.8f,%li,%li,%i,%i\n", i, 0, sizelib[i], CM[i].nblock, CM[i].ninterface, CM[i].nnonzero, CM[i].lzc, CM[i].necklace, 0., 0., sampledcount[i]/NSAMP, (long) 0, (long) 0, symms[i], bestmodularity[i]);
	}
      else if(occur[i] != 0)
	{
	  fprintf(fp, "%i,%i,%i,%i,%i,%i,%i,%.5f,%.5f,%.5f,%.8f,%li,%li,%i,%i\n", i, discoveredpheno[i], sizelib[i], CM[i].nblock, CM[i].ninterface, CM[i].nnonzero, CM[i].lzc, CM[i].necklace, meandtime[i]/discoveredpheno[i], mindtime[i], sampledcount[i]/NSAMP, occur[i], adaptoccur[i], symms[i], bestmodularity[i]);
	} 
    }
  fclose(fp);

  // output phenotypes for this ensemble
  sprintf(string, "lib-%s.txt", outstr);
  fp = fopen(string, "w");
  for(i = 0; i < numlib; i++)
    {
      fprintf(output, "Library %i\n", i);
      FileOutputGrid(output, &lib[ARR*ARR*i], ARR);
    }
  fclose(output);

  return 0;
}
