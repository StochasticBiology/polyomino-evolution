/********* routines for computing standard forms and complexities of rulesets *******/

typedef struct tagComplexityMeasures
{
  int nblock;
  int ninterface;
  int nnonzero;
  int lzc;
  double necklace;
} ComplexityMeasures;

int Modularity(int *g, int ARR, int MAXBLOCKS);
int LZComplexity(char *p_binarySequence, int p_maxTreeNodes);
void Canonical(int *g, int *t, int ARR);
int SymmCount(int *grid, int ARR);
int CompareGene(int *g1, int *g2);
double NecklaceComplexity(int *t, int l, int GENES);
void GenomeComplexityMultiStats(int *gtemp, int n, int nuc, ComplexityMeasures *CM, int LEN, int nbitcol, FILE *fp);
double GenomeReducer(int *r, int l, int *s, int tracker, int GENES, int ARR);
double nfunction(int n, int c);

// get the number of distinct blocks that appear in a polyomino structure
// used to compute modularity in post-processing
// takes the grid representation of a polymoni as input and counts block types
int Modularity(int *g, int ARR, int MAXBLOCKS)
{
  int appear[MAXBLOCKS];
  int i;
  int nappear, size;

  nappear = size = 0;
  for(i = 0; i < MAXBLOCKS; i++) appear[i] = 0;
  for(i = 0; i < ARR*ARR; i++)
    {
      if(g[i] >= 0)
	{
	  if(appear[g[i]] == 0) nappear++;
	  appear[g[i]] = 1;
	  size++;
	}
    }

  return nappear;
}

int LZComplexity(char *p_binarySequence, int p_maxTreeNodes)
{
  void **patternTree;
  void **currentNode;
  void **nextFreeNode;
  int nodeCount;
  int sequenceIndex;
  int currentDigit;

  nodeCount = 0;
  patternTree = malloc(sizeof(void*) * (p_maxTreeNodes << 1));
  currentNode = patternTree;
  nextFreeNode = patternTree + (sizeof(void*) << 1);
  currentNode[0] = NULL;
  currentNode[1] = NULL;
  sequenceIndex = 0;

  while (p_binarySequence[sequenceIndex])
    {
      currentDigit = p_binarySequence[sequenceIndex] - 48;
      if (NULL == currentNode[currentDigit])
	{
	  currentNode[currentDigit] = nextFreeNode;
	  nextFreeNode[0] = NULL;
	  nextFreeNode[1] = NULL;
	  nextFreeNode += (sizeof(void*) << 1);
	  currentNode = patternTree;
	  nodeCount++;
	}
      else
	{
	  currentNode = currentNode[currentDigit];
	}
      sequenceIndex++;
    }

  free(patternTree);
  return nodeCount;
}



// simply puts a grid representation of a polyomino into a standard form using only -1 and 1 pixels
void Canonical(int *g, int *t, int ARR)
{
  int mini, minj;
  int i, j;

  mini = minj = ARR;
  for(i = 0; i < ARR; i++)
    {
      for(j = 0; j < ARR; j++)
	{
	  if(g[j*ARR+i] >= 0)
	    {
	      if(i < mini) mini = i;
	      if(j < minj) minj = j;
	    }
	}
    }

  for(i = 0; i < ARR*ARR; i++)
    t[i] = -1;
  for(i = mini; i < ARR; i++)
    {
      for(j = minj; j < ARR; j++)
	{
	  t[(j-minj)*ARR+(i-mini)] = (g[j*ARR+i] == -1 ? -1 : 1);
	}
    }
}

// assign a symmetry group to a polyomino structure
int SymmCount(int *grid, int ARR)
{
  int T[ARR*ARR];
  int ttry, x, y;
  int success[8];
  int i, j, mini, minj;
  int fail;
  int symm;

  // easier just to comment the overall scheme
  // we apply the 8 possible symmetry transformations to the structure
  // for each, we compare the transformed to the original structure
  // recording where we find identities
  // the pattern of identities then gives the symmetry group
  // T is the transformed structure
  for(ttry = 0; ttry < 8; ttry++)
    {
      for(x = 0; x < ARR; x++)
	{
	  for(y = 0; y < ARR; y++)
	    {
	      switch(ttry)
		{
		case 0: T[x+ARR*y] = grid[x+ARR*y]; break;
		case 1: T[x+ARR*y] = grid[(ARR-1-x)+ARR*y]; break;
		case 2: T[x+ARR*y] = grid[ARR-1-x+ARR*(ARR-1-y)]; break;
		case 3: T[x+ARR*y] = grid[x+ARR*(ARR-1-y)]; break;
		case 4: T[x+ARR*y] = grid[y+ARR*x]; break;
		case 5: T[x+ARR*y] = grid[ARR-1-y+ARR*x]; break;
		case 6: T[x+ARR*y] = grid[ARR-1-y+ARR*(ARR-1-x)]; break;
		case 7: T[x+ARR*y] = grid[y+ARR*(ARR-1-x)]; break;
		}
	    }
	}
      mini = minj = ARR;
      for(i = 0; i < ARR; i++)
	{
	  for(j = 0; j < ARR; j++)
	    {
	      if(T[j*ARR+i] >= 0)
		{
		  if(i < mini) mini = i;
		  if(j < minj) minj = j;
		}
	    }
	}

      
      fail = 0;
      for(i = mini; i < ARR; i++)
	{
	  for(j = minj; j < ARR; j++)
	    {
	      if(T[j*ARR+i] == -1 && grid[(j-minj)*ARR+(i-mini)] != -1 || T[j*ARR+i] != -1 && grid[(j-minj)*ARR+(i-mini)] == -1)
		{
		  fail = 1;
		  break;
		}
	    }
	  if(fail == 1) break;
	}
      
      if(fail == 0)
	success[ttry] = 1;
      else success[ttry] = 0;
    }

  // assign a symmetry group based on the located identities
  if(success[1] && success[2] && success[3] && success[4] && success[5] && success[6] && success[7]) symm = 0;
  else if(success[2] && success[5] && success[7]) symm = 1;
  else if(success[1] && success[2] && success[3]) symm = 2;
  else if(success[2] && success[4] && success[6]) symm = 2;
  else if(success[2]) symm = 3;
  else if(success[4] || success[6]) symm = 4;
  else if(success[1] || success[3]) symm = 4;
  else symm = 5;

  return symm;
}

/******** the CompareGene/Complexity/GenomeReducer system is probably deprecated *******/
// most stuff uses the GenomeComplexity function 
int CompareGene(int *g1, int *g2)
{
  int off, i;
  int fail;

  for(off = 0; off < 4; off++)
    {
      fail = 0;
      for(i = 0; i < 4; i++)
	{
	  if(g1[i] != g2[(i+off)%4]) { fail = 1; break; }
	}
      if(fail == 0) return 0;
    }
  return 1;
}

double NecklaceComplexity(int *t, int l, int GENES)
{
  int used[GENES/4];
  int i, j;
  int numused;
  double complexity = 0;
  int cols[100];
  int highest;
  double ci, nprime;
  int found;

  numused = 0;

  for(i = 0; i < GENES/4; i++)
    used[i] = 0;

  for(i = 0; i < GENES/4; i++)
    {
      found = 0;
      for(j = 0; j < numused; j++)
	{
	  if(CompareGene(&(t[4*i]), &(t[4*used[j]])) == 0)
	    {
	      found = 1;
	      break;
	    }
	}
      if(found == 0)
	{
	  for(j = 0; j < 4; j++)
	    used[numused] = i;
	  numused++;
	}
    }

  highest = 0;

  for(i = 0; i < numused; i++)
    {
      for(j = 0; j < 4; j++)
	{
	  if(t[4*used[i]+j] > highest) highest = t[4*used[i]+j];
	}
    }

  highest++;
  for(i = 0; i < numused; i++)
    {
      if(used[i] == 4 && (l == 3 || l == 4))
	printf(".");
      ci = 0;
      for(j = 0; j < highest; j++)
	cols[j] = 0;
      for(j = 0; j < 4; j++)      
        cols[t[4*used[i]+j]]++;
      for(j = 0; j < highest; j++)
	ci += (cols[j] != 0);

      switch((int)ci)
	{
	case 1: nprime = 1; break;
	case 2: nprime = 4; break;
	case 3: nprime = 9; break;
	case 4: nprime = 6; break;
	}
      if(!(cols[0] > 0 && ci == 1))
        complexity += ci*log2(highest)+log2(nprime);
    }

  return complexity;
}

/***************/

// GenomeComplexity
// takes: g - list of sides (genome); n - number of blocks in genome; nuc - block that acts as nucleus
// returns effective complexity of genome
void GenomeComplexityMultiStats(int *gtemp, int n, int nuc, ComplexityMeasures *CM, int LEN, int nbitcol, FILE *fp)
{
  int used[100];
  int usedblocks[100];
  int onthis[100];
  int promote[100];
  int numonthis;
  int i, j, k, l;
  int change;
  double nc;
  double nprime;
  int partner;
  int g[LEN];
  int nblock, ninterface, nnonzero, lzc;
  double necklace;
  int totalnum;
  char genomestr[LEN*nbitcol];
  int gslen = 0;
  int thisval, currval;
  
  for(i = 0; i < LEN; i++)
    g[i] = gtemp[i];

  for(i = 0; i < 100; i++)
    used[i] = usedblocks[i] = 0;

  // label initial donor bonds from nucleus
  for(i = nuc*4; i < nuc*4+4; i++)
    used[g[i]] = 1;
  usedblocks[nuc] = 1;

  // iterate: each step, "donor" bonds (labelled 1 in used) are investigated. if their partners are present
  // anywhere, the "donor" and "acceptor" are promoted at the end of the iteration to a label 2 in used.
  do{
    change = 0;
    for(i = 1; i < 100; i++) 
      promote[i] = 0;
    for(i = 1; i < 100; i++)
      {
	if(used[i] == 1)
	  {
	    if(i % 2 == 0) partner = i-1;
	    else partner = i+1;
	    // search for "acceptor" partner
            for(j = 0; j < n; j++)
	      {
		for(k = 0; k < 4; k++)
		  {
		    if(g[4*j+k] == partner) 
		      {
			change = 1;
			promote[i] = 1; promote[partner] = 1;
			usedblocks[j] = 1;
			// if an acceptor is found, make all other bonds on that block new donors
			for(l = 0; l < 4; l++)
			  {
			    if(used[g[4*j+l]] == 0 && l != k) 
			      used[g[4*j+l]] = 1;
			  }
		      }
		  }
	      }
	  }
      }
    // promote stored bonds
    for(i = 1; i < 100; i++)
      {
	if(promote[i]) used[i] = 2;
      }
  }while(change == 1);

  // remove unused bonds from genome
  for(i = 0; i < 4*n; i++)
    {
      if(used[g[i]] != 2) g[i] = 0;
    }
 
  // count number of colour pairs
  nc = 0;
  for(i = 0; i < 100; i++)
    nc += (used[i]==2);
  nc /= 2;
  ninterface = nc;
  
  nblock = nnonzero = 0;
  necklace = 0;
  
  for(i = 0; i < n; i++)
    {
      // for each used block, add complexity value to total
      if(usedblocks[i])
	{
	  nblock++;
	  for(j = 0; j < 100; j++)
	    onthis[j] = 0;
	  numonthis = 0;
	  for(j = 0; j < 4; j++)
	    {
	      thisval = g[4*i+j];
	      currval = 1; for(k = 0; k < nbitcol-1; k++) currval *= 2;
	      for(k = 0; k < nbitcol; k++)
		{
		  if(thisval >= currval)
		    {
		      genomestr[gslen++] = '1';
		      thisval -= currval;
		    }
		  else genomestr[gslen++] = '0';
		  currval /= 2;
		}
	      // comment this to stop output
	      fprintf(fp, "%i ", g[4*i+j]);
	      if(g[4*i+j] != 0) { totalnum++; nnonzero++; }
	      if(onthis[g[4*i+j]] == 0)
		{
		  onthis[g[4*i+j]] = 1;
		  numonthis++;
		}
	    }
	  nprime = nfunction(4, numonthis);
	  necklace += numonthis*log2(nc)+log2(nprime);
	  //	  printf("(%i %.5f) ", numonthis, nprime);
	}
    }
  fprintf(fp, "\n");
  //  printf("%.5f ", nc);

  genomestr[gslen] = '\0';
  fprintf(fp, "%s\n", genomestr);
  
  lzc = LZComplexity(genomestr, 1000);
  
  // comment this to stop printout
  //  printf("\n%.5f\n", complexity);

  CM->nblock = nblock;
  CM->ninterface = ninterface;
  CM->nnonzero = nnonzero;
  CM->necklace = necklace;
  CM->lzc = lzc;
  
  //totalnum;//*log2(nc);

  /******** the above bit is fluid -- but these days we're most interested in just the number of colours in the genome. this is nc*2 -- at the moment we're post-processing this by e.g. plotting 2*(iaincomp2) *********/
}


/**** this function belongs to the CompareGene/Complexity/GenomeReducer system ***/
// reduces a genome to minimal effective interactions
// we now favour optimisation bare genomes across a sampled set of structures
double GenomeReducer(int *r, int l, int *s, int tracker, int GENES, int ARR)
{
  int t[GENES];
  int involved[GENES/4];
  int i, j;
  int label, pfound, bond, partner;
  int done[100];
  FILE *fp;
  int size;
  int grid[ARR*ARR];

  if(tracker == 0)
    {
      fp = fopen("interest.gen", "a");
      fprintf(fp, "%i\n", l);
    }

  for(i = 0; i < 100; i++)
    done[i] = 0;

  for(i = 0; i < GENES/4; i++)
    involved[i] = 0;

  for(i = 0; i < GENES; i++)
    {
      s[i] = 0;
      for(j = 0; j < 4; j++)
	{
	  s[i] += r[4*i+j]*pow(2, 3-j);
	}
    }

  if(Grow(s, GENES/4, grid, &size, 0, 0, 10, 0, ARR) == 0)
    {
      for(i = 0; i < ARR*ARR; i++)
	{
	  if(grid[i] >= 0)
	    involved[grid[i]] = 1;
	}
    }
  else
    {
      for(i = 0; i < GENES/4; i++)
	involved[i] = 1;
    }


  for(i = 0; i < GENES; i++)
    {
      //     fprintf(fp, "%i, ", s[i]);
      t[i] = s[i];
      //      if((i+1) % 4 == 0) fprintf(fp, "| ");
    }
  //fprintf(fp, "\n");

  label = 1;
  for(i = 0; i < GENES; i++)
    {
      if(!involved[i/4]) t[i] = 0;
      else if(s[i] && done[s[i]] == 0)
	{
	  pfound = 0;
	  bond = s[i]; partner = (bond % 2 ? bond+1 : bond-1);
	  for(j = 0; j < GENES; j++)
	    {
	      if(s[j] == bond) t[j] = label;
	      if(s[j] == partner) { t[j] = label+1; pfound = 1; }
	    }
          done[bond] = 1; done[partner] = 1;
	  if(pfound == 0)
	    {
	      for(j = 0; j < GENES; j++)
		{
		  if(s[j] == bond) t[j] = 0;
		}
	    }
	  else {  label += 2; }
	}
    }
  if(tracker == 0)
    {
      for(i = 0; i < GENES; i++)
	{
	  fprintf(fp, "%i, ", t[i]);
	  //  if((i+1) % 4 == 0) fprintf(fp, "| ");
	}
      fclose(fp);
    }
  if(tracker == 1)
    {
      for(i = 0; i < GENES; i++)
	s[i] = t[i];
    }
  /*  else
      {
      fp = fopen("tracker.dat", "a");
      for(i = 0; i < GENES; i++)
      {
      fprintf(fp, "%i, ", t[i]);
      //  if((i+1) % 4 == 0) fprintf(fp, "| ");
      }
      fprintf(fp, "\n%.5f\n", Complexity(t));
      fclose(fp);
      }*/

  return NecklaceComplexity(t, l, GENES);
}

double nfunction(int n, int c)
{
  if(n ==  1  && c ==  1 ) return  1. ;
  if(n ==  1  && c ==  2 ) return  2. ;
  if(n ==  1  && c ==  3 ) return  3. ;
  if(n ==  1  && c ==  4 ) return  4. ;
  if(n ==  1  && c ==  5 ) return  5. ;
  if(n ==  1  && c ==  6 ) return  6. ;
  if(n ==  1  && c ==  7 ) return  7. ;
  if(n ==  1  && c ==  8 ) return  8. ;
  if(n ==  1  && c ==  9 ) return  9. ;
  if(n ==  1  && c ==  10 ) return  10. ;
  if(n ==  1  && c ==  11 ) return  11. ;
  if(n ==  1  && c ==  12 ) return  12. ;
  if(n ==  1  && c ==  13 ) return  13. ;
  if(n ==  1  && c ==  14 ) return  14. ;
  if(n ==  1  && c ==  15 ) return  15. ;
  if(n ==  1  && c ==  16 ) return  16. ;
  if(n ==  1  && c ==  17 ) return  17. ;
  if(n ==  1  && c ==  18 ) return  18. ;
  if(n ==  1  && c ==  19 ) return  19. ;
  if(n ==  1  && c ==  20 ) return  20. ;
  if(n ==  2  && c ==  1 ) return  1. ;
  if(n ==  2  && c ==  2 ) return  3. ;
  if(n ==  2  && c ==  3 ) return  6. ;
  if(n ==  2  && c ==  4 ) return  10. ;
  if(n ==  2  && c ==  5 ) return  15. ;
  if(n ==  2  && c ==  6 ) return  21. ;
  if(n ==  2  && c ==  7 ) return  28. ;
  if(n ==  2  && c ==  8 ) return  36. ;
  if(n ==  2  && c ==  9 ) return  45. ;
  if(n ==  2  && c ==  10 ) return  55. ;
  if(n ==  2  && c ==  11 ) return  66. ;
  if(n ==  2  && c ==  12 ) return  78. ;
  if(n ==  2  && c ==  13 ) return  91. ;
  if(n ==  2  && c ==  14 ) return  105. ;
  if(n ==  2  && c ==  15 ) return  120. ;
  if(n ==  2  && c ==  16 ) return  136. ;
  if(n ==  2  && c ==  17 ) return  153. ;
  if(n ==  2  && c ==  18 ) return  171. ;
  if(n ==  2  && c ==  19 ) return  190. ;
  if(n ==  2  && c ==  20 ) return  210. ;
  if(n ==  3  && c ==  1 ) return  1. ;
  if(n ==  3  && c ==  2 ) return  4. ;
  if(n ==  3  && c ==  3 ) return  11. ;
  if(n ==  3  && c ==  4 ) return  24. ;
  if(n ==  3  && c ==  5 ) return  45. ;
  if(n ==  3  && c ==  6 ) return  76. ;
  if(n ==  3  && c ==  7 ) return  119. ;
  if(n ==  3  && c ==  8 ) return  176. ;
  if(n ==  3  && c ==  9 ) return  249. ;
  if(n ==  3  && c ==  10 ) return  340. ;
  if(n ==  3  && c ==  11 ) return  451. ;
  if(n ==  3  && c ==  12 ) return  584. ;
  if(n ==  3  && c ==  13 ) return  741. ;
  if(n ==  3  && c ==  14 ) return  924. ;
  if(n ==  3  && c ==  15 ) return  1135. ;
  if(n ==  3  && c ==  16 ) return  1376. ;
  if(n ==  3  && c ==  17 ) return  1649. ;
  if(n ==  3  && c ==  18 ) return  1956. ;
  if(n ==  3  && c ==  19 ) return  2299. ;
  if(n ==  3  && c ==  20 ) return  2680. ;
  if(n ==  4  && c ==  1 ) return  1. ;
  if(n ==  4  && c ==  2 ) return  6. ;
  if(n ==  4  && c ==  3 ) return  24. ;
  if(n ==  4  && c ==  4 ) return  70. ;
  if(n ==  4  && c ==  5 ) return  165. ;
  if(n ==  4  && c ==  6 ) return  336. ;
  if(n ==  4  && c ==  7 ) return  616. ;
  if(n ==  4  && c ==  8 ) return  1044. ;
  if(n ==  4  && c ==  9 ) return  1665. ;
  if(n ==  4  && c ==  10 ) return  2530. ;
  if(n ==  4  && c ==  11 ) return  3696. ;
  if(n ==  4  && c ==  12 ) return  5226. ;
  if(n ==  4  && c ==  13 ) return  7189. ;
  if(n ==  4  && c ==  14 ) return  9660. ;
  if(n ==  4  && c ==  15 ) return  12720. ;
  if(n ==  4  && c ==  16 ) return  16456. ;
  if(n ==  4  && c ==  17 ) return  20961. ;
  if(n ==  4  && c ==  18 ) return  26334. ;
  if(n ==  4  && c ==  19 ) return  32680. ;
  if(n ==  4  && c ==  20 ) return  40110. ;
  if(n ==  5  && c ==  1 ) return  1. ;
  if(n ==  5  && c ==  2 ) return  8. ;
  if(n ==  5  && c ==  3 ) return  51. ;
  if(n ==  5  && c ==  4 ) return  208. ;
  if(n ==  5  && c ==  5 ) return  629. ;
  if(n ==  5  && c ==  6 ) return  1560. ;
  if(n ==  5  && c ==  7 ) return  3367. ;
  if(n ==  5  && c ==  8 ) return  6560. ;
  if(n ==  5  && c ==  9 ) return  11817. ;
  if(n ==  5  && c ==  10 ) return  20008. ;
  if(n ==  5  && c ==  11 ) return  32219. ;
  if(n ==  5  && c ==  12 ) return  49776. ;
  if(n ==  5  && c ==  13 ) return  74269. ;
  if(n ==  5  && c ==  14 ) return  107576. ;
  if(n ==  5  && c ==  15 ) return  151887. ;
  if(n ==  5  && c ==  16 ) return  209728. ;
  if(n ==  5  && c ==  17 ) return  283985. ;
  if(n ==  5  && c ==  18 ) return  377928. ;
  if(n ==  5  && c ==  19 ) return  495235. ;
  if(n ==  5  && c ==  20 ) return  640016. ;
  if(n ==  6  && c ==  1 ) return  1. ;
  if(n ==  6  && c ==  2 ) return  14. ;
  if(n ==  6  && c ==  3 ) return  130. ;
  if(n ==  6  && c ==  4 ) return  700. ;
  if(n ==  6  && c ==  5 ) return  2635. ;
  if(n ==  6  && c ==  6 ) return  7826. ;
  if(n ==  6  && c ==  7 ) return  19684. ;
  if(n ==  6  && c ==  8 ) return  43800. ;
  if(n ==  6  && c ==  9 ) return  88725. ;
  if(n ==  6  && c ==  10 ) return  166870. ;
  if(n ==  6  && c ==  11 ) return  295526. ;
  if(n ==  6  && c ==  12 ) return  498004. ;
  if(n ==  6  && c ==  13 ) return  804895. ;
  if(n ==  6  && c ==  14 ) return  1255450. ;
  if(n ==  6  && c ==  15 ) return  1899080. ;
  if(n ==  6  && c ==  16 ) return  2796976. ;
  if(n ==  6  && c ==  17 ) return  4023849. ;
  if(n ==  6  && c ==  18 ) return  5669790. ;
  if(n ==  6  && c ==  19 ) return  7842250. ;
  if(n ==  6  && c ==  20 ) return  10668140. ;
  if(n ==  7  && c ==  1 ) return  1. ;
  if(n ==  7  && c ==  2 ) return  20. ;
  if(n ==  7  && c ==  3 ) return  315. ;
  if(n ==  7  && c ==  4 ) return  2344. ;
  if(n ==  7  && c ==  5 ) return  11165. ;
  if(n ==  7  && c ==  6 ) return  39996. ;
  if(n ==  7  && c ==  7 ) return  117655. ;
  if(n ==  7  && c ==  8 ) return  299600. ;
  if(n ==  7  && c ==  9 ) return  683289. ;
  if(n ==  7  && c ==  10 ) return  1428580. ;
  if(n ==  7  && c ==  11 ) return  2783891. ;
  if(n ==  7  && c ==  12 ) return  5118840. ;
  if(n ==  7  && c ==  13 ) return  8964085. ;
  if(n ==  7  && c ==  14 ) return  15059084. ;
  if(n ==  7  && c ==  15 ) return  24408495. ;
  if(n ==  7  && c ==  16 ) return  38347936. ;
  if(n ==  7  && c ==  17 ) return  58619825. ;
  if(n ==  7  && c ==  18 ) return  87460020. ;
  if(n ==  7  && c ==  19 ) return  127695979. ;
  if(n ==  7  && c ==  20 ) return  182857160. ;
  if(n ==  8  && c ==  1 ) return  1. ;
  if(n ==  8  && c ==  2 ) return  36. ;
  if(n ==  8  && c ==  3 ) return  834. ;
  if(n ==  8  && c ==  4 ) return  8230. ;
  if(n ==  8  && c ==  5 ) return  48915. ;
  if(n ==  8  && c ==  6 ) return  210126. ;
  if(n ==  8  && c ==  7 ) return  720916. ;
  if(n ==  8  && c ==  8 ) return  2097684. ;
  if(n ==  8  && c ==  9 ) return  5381685. ;
  if(n ==  8  && c ==  10 ) return  12501280. ;
  if(n ==  8  && c ==  11 ) return  26796726. ;
  if(n ==  8  && c ==  12 ) return  53750346. ;
  if(n ==  8  && c ==  13 ) return  101969959. ;
  if(n ==  8  && c ==  14 ) return  184478490. ;
  if(n ==  8  && c ==  15 ) return  320367720. ;
  if(n ==  8  && c ==  16 ) return  536879176. ;
  if(n ==  8  && c ==  17 ) return  871980201. ;
  if(n ==  8  && c ==  18 ) return  1377508284. ;
  if(n ==  8  && c ==  19 ) return  2122961770. ;
  if(n ==  8  && c ==  20 ) return  3200020110. ;
  if(n ==  9  && c ==  1 ) return  1. ;
  if(n ==  9  && c ==  2 ) return  60. ;
  if(n ==  9  && c ==  3 ) return  2195. ;
  if(n ==  9  && c ==  4 ) return  29144. ;
  if(n ==  9  && c ==  5 ) return  217045. ;
  if(n ==  9  && c ==  6 ) return  1119796. ;
  if(n ==  9  && c ==  7 ) return  4483815. ;
  if(n ==  9  && c ==  8 ) return  14913200. ;
  if(n ==  9  && c ==  9 ) return  43046889. ;
  if(n ==  9  && c ==  10 ) return  111111340. ;
  if(n ==  9  && c ==  11 ) return  261994491. ;
  if(n ==  9  && c ==  12 ) return  573309320. ;
  if(n ==  9  && c ==  13 ) return  1178278205. ;
  if(n ==  9  && c ==  14 ) return  2295672484. ;
  if(n ==  9  && c ==  15 ) return  4271485135. ;
  if(n ==  9  && c ==  16 ) return  7635498336. ;
  if(n ==  9  && c ==  17 ) return  13176431825. ;
  if(n ==  9  && c ==  18 ) return  22039922460. ;
  if(n ==  9  && c ==  19 ) return  35854190179. ;
  if(n ==  9  && c ==  20 ) return  56888890680. ;
  if(n ==  10  && c ==  1 ) return  1. ;
  if(n ==  10  && c ==  2 ) return  108. ;
  if(n ==  10  && c ==  3 ) return  5934. ;
  if(n ==  10  && c ==  4 ) return  104968. ;
  if(n ==  10  && c ==  5 ) return  976887. ;
  if(n ==  10  && c ==  6 ) return  6047412. ;
  if(n ==  10  && c ==  7 ) return  28249228. ;
  if(n ==  10  && c ==  8 ) return  107377488. ;
  if(n ==  10  && c ==  9 ) return  348684381. ;
  if(n ==  10  && c ==  10 ) return  1000010044. ;
  if(n ==  10  && c ==  11 ) return  2593758618. ;
  if(n ==  10  && c ==  12 ) return  6191761368. ;
  if(n ==  10  && c ==  13 ) return  13785886387. ;
  if(n ==  10  && c ==  14 ) return  28925519364. ;
  if(n ==  10  && c ==  15 ) return  57665115096. ;
  if(n ==  10  && c ==  16 ) return  109951267744. ;
  if(n ==  10  && c ==  17 ) return  201599532153. ;
  if(n ==  10  && c ==  18 ) return  357046911756. ;
  if(n ==  10  && c ==  19 ) return  613106873542. ;
  if(n ==  10  && c ==  20 ) return  1024000320168. ;
  if(n ==  11  && c ==  1 ) return  1. ;
  if(n ==  11  && c ==  2 ) return  188. ;
  if(n ==  11  && c ==  3 ) return  16107. ;
  if(n ==  11  && c ==  4 ) return  381304. ;
  if(n ==  11  && c ==  5 ) return  4438925. ;
  if(n ==  11  && c ==  6 ) return  32981556. ;
  if(n ==  11  && c ==  7 ) return  179756983. ;
  if(n ==  11  && c ==  8 ) return  780903152. ;
  if(n ==  11  && c ==  9 ) return  2852823609. ;
  if(n ==  11  && c ==  10 ) return  9090909100. ;
  if(n ==  11  && c ==  11 ) return  25937424611. ;
  if(n ==  11  && c ==  12 ) return  67546215528. ;
  if(n ==  11  && c ==  13 ) return  162923672197. ;
  if(n ==  11  && c ==  14 ) return  368142288164. ;
  if(n ==  11  && c ==  15 ) return  786341441775. ;
  if(n ==  11  && c ==  16 ) return  1599289640416. ;
  if(n ==  11  && c ==  17 ) return  3115626937073. ;
  if(n ==  11  && c ==  18 ) return  5842582734492. ;
  if(n ==  11  && c ==  19 ) return  10590023536219. ;
  if(n ==  11  && c ==  20 ) return  18618181818200. ;
  if(n ==  12  && c ==  1 ) return  1. ;
  if(n ==  12  && c ==  2 ) return  352. ;
  if(n ==  12  && c ==  3 ) return  44368. ;
  if(n ==  12  && c ==  4 ) return  1398500. ;
  if(n ==  12  && c ==  5 ) return  20346485. ;
  if(n ==  12  && c ==  6 ) return  181402676. ;
  if(n ==  12  && c ==  7 ) return  1153450872. ;
  if(n ==  12  && c ==  8 ) return  5726645688. ;
  if(n ==  12  && c ==  9 ) return  23535840225. ;
  if(n ==  12  && c ==  10 ) return  83333418520. ;
  if(n ==  12  && c ==  11 ) return  261535848376. ;
  if(n ==  12  && c ==  12 ) return  743008623292. ;
  if(n ==  12  && c ==  13 ) return  1941507500933. ;
  if(n ==  12  && c ==  14 ) return  4724493332300. ;
  if(n ==  12  && c ==  15 ) return  10812195782480. ;
  if(n ==  12  && c ==  16 ) return  23456249468976. ;
  if(n ==  12  && c ==  17 ) return  48551855128737. ;
  if(n ==  12  && c ==  18 ) return  96402617971728. ;
  if(n ==  12  && c ==  19 ) return  184442913865600. ;
  if(n ==  12  && c ==  20 ) return  341333338694740. ;
  if(n ==  13  && c ==  1 ) return  1. ;
  if(n ==  13  && c ==  2 ) return  632. ;
  if(n ==  13  && c ==  3 ) return  122643. ;
  if(n ==  13  && c ==  4 ) return  5162224. ;
  if(n ==  13  && c ==  5 ) return  93900245. ;
  if(n ==  13  && c ==  6 ) return  1004668776. ;
  if(n ==  13  && c ==  7 ) return  7453000807. ;
  if(n ==  13  && c ==  8 ) return  42288908768. ;
  if(n ==  13  && c ==  9 ) return  195528140649. ;
  if(n ==  13  && c ==  10 ) return  769230769240. ;
  if(n ==  13  && c ==  11 ) return  2655593241851. ;
  if(n ==  13  && c ==  12 ) return  8230246567632. ;
  if(n ==  13  && c ==  13 ) return  23298085122493. ;
  if(n ==  13  && c ==  14 ) return  61054982558024. ;
  if(n ==  13  && c ==  15 ) return  149707312950735. ;
  if(n ==  13  && c ==  16 ) return  346430740566976. ;
  if(n ==  13  && c ==  17 ) return  761890617915857. ;
  if(n ==  13  && c ==  18 ) return  1601766528128568. ;
  if(n ==  13  && c ==  19 ) return  3234844881712099. ;
  if(n ==  13  && c ==  20 ) return  6301538461538480. ;
  if(n ==  14  && c ==  1 ) return  1. ;
  if(n ==  14  && c ==  2 ) return  1182. ;
  if(n ==  14  && c ==  3 ) return  341802. ;
  if(n ==  14  && c ==  4 ) return  19175140. ;
  if(n ==  14  && c ==  5 ) return  435970995. ;
  if(n ==  14  && c ==  6 ) return  5597460306. ;
  if(n ==  14  && c ==  7 ) return  48444564052. ;
  if(n ==  14  && c ==  8 ) return  314146329192. ;
  if(n ==  14  && c ==  9 ) return  1634056945605. ;
  if(n ==  14  && c ==  10 ) return  7142857857190. ;
  if(n ==  14  && c ==  11 ) return  27124989505086. ;
  if(n ==  14  && c ==  12 ) return  91708464312972. ;
  if(n ==  14  && c ==  13 ) return  281241174889207. ;
  if(n ==  14  && c ==  14 ) return  793714780783770. ;
  if(n ==  14  && c ==  15 ) return  2085209014017960. ;
  if(n ==  14  && c ==  16 ) return  5146971021883216. ;
  if(n ==  14  && c ==  17 ) return  12026987640695817. ;
  if(n ==  14  && c ==  18 ) return  26772383442450222. ;
  if(n ==  14  && c ==  19 ) return  57071906191197010. ;
  if(n ==  14  && c ==  20 ) return  117028571520000180. ;
  if(n ==  15  && c ==  1 ) return  1. ;
  if(n ==  15  && c ==  2 ) return  2192. ;
  if(n ==  15  && c ==  3 ) return  956635. ;
  if(n ==  15  && c ==  4 ) return  71582944. ;
  if(n ==  15  && c ==  5 ) return  2034505661. ;
  if(n ==  15  && c ==  6 ) return  31345666736. ;
  if(n ==  15  && c ==  7 ) return  316504102999. ;
  if(n ==  15  && c ==  8 ) return  2345624810432. ;
  if(n ==  15  && c ==  9 ) return  13726075481049. ;
  if(n ==  15  && c ==  10 ) return  66666666680272. ;
  if(n ==  15  && c ==  11 ) return  278483211316211. ;
  if(n ==  15  && c ==  12 ) return  1027134771672736. ;
  if(n ==  15  && c ==  13 ) return  3412392867656149. ;
  if(n ==  15  && c ==  14 ) return  10371206370593264. ;
  if(n ==  15  && c ==  15 ) return  29192926025492783. ;
  if(n ==  15  && c ==  16 ) return  76861433640597376. ;
  if(n ==  15  && c ==  17 ) return  190828203434178353. ;
  if(n ==  15  && c ==  18 ) return  449776041098750736. ;
  if(n ==  15  && c ==  19 ) return  1012075135325318539. ;
  if(n ==  15  && c ==  20 ) return  2184533333333762144. ;
  if(n ==  16  && c ==  1 ) return  1. ;
  if(n ==  16  && c ==  2 ) return  4116. ;
  if(n ==  16  && c ==  3 ) return  2690844. ;
  if(n ==  16  && c ==  4 ) return  268439590. ;
  if(n ==  16  && c ==  5 ) return  9536767665. ;
  if(n ==  16  && c ==  6 ) return  176319474366. ;
  if(n ==  16  && c ==  7 ) return  2077058521216. ;
  if(n ==  16  && c ==  8 ) return  17592187093524. ;
  if(n ==  16  && c ==  9 ) return  115813764494505. ;
  if(n ==  16  && c ==  10 ) return  625000006251280. ;
  if(n ==  16  && c ==  11 ) return  2871858129872556. ;
  if(n ==  16  && c ==  12 ) return  11555266207816266. ;
  if(n ==  16  && c ==  13 ) return  41588538124935529. ;
  if(n ==  16  && c ==  14 ) return  136122083705327370. ;
  if(n ==  16  && c ==  15 ) return  410525522392242720. ;
  if(n ==  16  && c ==  16 ) return  1152921504875290696. ;
  if(n ==  16  && c ==  17 ) return  3041324492665174641. ;
  if(n ==  16  && c ==  18 ) return  7589970694225901484. ;
  if(n ==  16  && c ==  19 ) return  18027588349037812060. ;
  if(n ==  16  && c ==  20 ) return  40960000001600020110. ;
  if(n ==  17  && c ==  1 ) return  1. ;
  if(n ==  17  && c ==  2 ) return  7712. ;
  if(n ==  17  && c ==  3 ) return  7596483. ;
  if(n ==  17  && c ==  4 ) return  1010580544. ;
  if(n ==  17  && c ==  5 ) return  44878791365. ;
  if(n ==  17  && c ==  6 ) return  995685849696. ;
  if(n ==  17  && c ==  7 ) return  13684147881607. ;
  if(n ==  17  && c ==  8 ) return  132458812569728. ;
  if(n ==  17  && c ==  9 ) return  981010688215689. ;
  if(n ==  17  && c ==  10 ) return  5882352941176480. ;
  if(n ==  17  && c ==  11 ) return  29732178147017291. ;
  if(n ==  17  && c ==  12 ) return  130506535690613952. ;
  if(n ==  17  && c ==  13 ) return  508847995257725773. ;
  if(n ==  17  && c ==  14 ) return  1793608631137129184. ;
  if(n ==  17  && c ==  15 ) return  5795654431511374095. ;
  if(n ==  17  && c ==  16 ) return  17361641481138401536. ;
  if(n ==  17  && c ==  17 ) return  48661191875666868497. ;
  if(n ==  17  && c ==  18 ) return  128583032925805678368. ;
  if(n ==  17  && c ==  19 ) return  322375697516753069779. ;
  if(n ==  17  && c ==  20 ) return  771011764705882352960. ;
  if(n ==  18  && c ==  1 ) return  1. ;
  if(n ==  18  && c ==  2 ) return  14602. ;
  if(n ==  18  && c ==  3 ) return  21524542. ;
  if(n ==  18  && c ==  4 ) return  3817763740. ;
  if(n ==  18  && c ==  5 ) return  211927736135. ;
  if(n ==  18  && c ==  6 ) return  5642220380006. ;
  if(n ==  18  && c ==  7 ) return  90467424361132. ;
  if(n ==  18  && c ==  8 ) return  1000799924679192. ;
  if(n ==  18  && c ==  9 ) return  8338590871415805. ;
  if(n ==  18  && c ==  10 ) return  55555555611222370. ;
  if(n ==  18  && c ==  11 ) return  308884295325206986. ;
  if(n ==  18  && c ==  12 ) return  1479074071447277812. ;
  if(n ==  18  && c ==  13 ) return  6247522609031752867. ;
  if(n ==  18  && c ==  14 ) return  23715491901739603070. ;
  if(n ==  18  && c ==  15 ) return  82105104448548141080. ;
  if(n ==  18  && c ==  16 ) return  262353693496577680816. ;
  if(n ==  18  && c ==  17 ) return  781282469565908953017. ;
  if(n ==  18  && c ==  18 ) return  2185911559749720272442. ;
  if(n ==  18  && c ==  19 ) return  5784852794346334629910. ;
  if(n ==  18  && c ==  20 ) return  14563555555584007112140. ;
  if(n ==  19  && c ==  1 ) return  1. ;
  if(n ==  19  && c ==  2 ) return  27596. ;
  if(n ==  19  && c ==  3 ) return  61171659. ;
  if(n ==  19  && c ==  4 ) return  14467258264. ;
  if(n ==  19  && c ==  5 ) return  1003867701485. ;
  if(n ==  19  && c ==  6 ) return  32071565263716. ;
  if(n ==  19  && c ==  7 ) return  599941851861751. ;
  if(n ==  19  && c ==  8 ) return  7585009898729264. ;
  if(n ==  19  && c ==  9 ) return  71097458824894329. ;
  if(n ==  19  && c ==  10 ) return  526315789473684220. ;
  if(n ==  19  && c ==  11 ) return  3218899497284976131. ;
  if(n ==  19  && c ==  12 ) return  16814736808980154056. ;
  if(n ==  19  && c ==  13 ) return  76943173177655058469. ;
  if(n ==  19  && c ==  14 ) return  314542313628890231444. ;
  if(n ==  19  && c ==  15 ) return  1166756747396368729455. ;
  if(n ==  19  && c ==  16 ) return  3976729669784964390496. ;
  if(n ==  19  && c ==  17 ) return  12582759772902701307761. ;
  if(n ==  19  && c ==  18 ) return  37275544492386193492524. ;
  if(n ==  19  && c ==  19 ) return  104127350297911241532859. ;
  if(n ==  19  && c ==  20 ) return  275941052631578947368440. ;
  if(n ==  20  && c ==  1 ) return  1. ;
  if(n ==  20  && c ==  2 ) return  52488. ;
  if(n ==  20  && c ==  3 ) return  174342216. ;
  if(n ==  20  && c ==  4 ) return  54975633976. ;
  if(n ==  20  && c ==  5 ) return  4768372070757. ;
  if(n ==  20  && c ==  6 ) return  182807925027504. ;
  if(n ==  20  && c ==  7 ) return  3989613329006536. ;
  if(n ==  20  && c ==  8 ) return  57646075284033552. ;
  if(n ==  20  && c ==  9 ) return  607883273127192897. ;
  if(n ==  20  && c ==  10 ) return  5000000000500012024. ;
  if(n ==  20  && c ==  11 ) return  33637499747924890752. ;
  if(n ==  20  && c ==  12 ) return  191687999625469653384. ;
  if(n ==  20  && c ==  13 ) return  950248188750932939413. ;
  if(n ==  20  && c ==  14 ) return  4183412771278702872288. ;
  if(n ==  20  && c ==  15 ) return  16626283650427087000176. ;
  if(n ==  20  && c ==  16 ) return  60446290980786434434720. ;
  if(n ==  20  && c ==  17 ) return  203211570332479425973857. ;
  if(n ==  20  && c ==  18 ) return  637411810819982432293224. ;
  if(n ==  20  && c ==  19 ) return  1879498672877604463254424. ;
  if(n ==  20  && c ==  20 ) return  5242880000000512000352088. ;
  return 0;
}
