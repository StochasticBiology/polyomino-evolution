void BreedGen(int *P, float *f, int NPAR, int LEN, float mu);

// core GA routine: takes a population P and vector of fitness values f, performs selection and mutation (rate mu)
void BreedGen(int *P, float *f, int NPAR, int LEN, float mu)
{
  int *tp;
  float total;
  int i, j;
  float r;
  int selected;
  float *cumsum;

  // allocate memory for cumulative sum of fitnesses
  cumsum = (float*)malloc(sizeof(float)*NPAR);
  // tp will store our temporary new population then get copied at the end
  tp = (int*)malloc(sizeof(int)*NPAR*LEN);

  total = 0;
  // construct roulette wheel from fitness values (cumsum = cumulative sum of fitnesses)
  for(i = 0; i < NPAR; i++)
    {
      cumsum[i] = (i == 0 ? f[i] : f[i]+cumsum[i-1]);
      total += f[i];
    }
  if(total == 0)
    {
      for(i = 0; i < NPAR; i++)
	cumsum[i] = i;
      total = i;
    }

  selected = NPAR-1;
  // loop through the members of the next generation
  for(i = 0; i < NPAR; i++)
    {
      // roll roulette ball
      r = drand48()*total;
      for(j = 0; j < NPAR; j++)
	{
	  if(cumsum[j] >= r)
	    {
	      selected = j;
	      break;
	    }
	}
      // copy selected genome with point mutations
      // this isn't an efficient way of applying point mutations: but this isn't the code bottleneck
      for(j = 0; j < LEN; j++)
	{
	  tp[i*LEN+j] = (drand48() < mu ? !P[selected*LEN+j] : P[selected*LEN+j]);
	}
    }
  for(i = 0; i < NPAR*LEN; i++)
    P[i] = tp[i];

  free(tp);
  free(cumsum);
}
