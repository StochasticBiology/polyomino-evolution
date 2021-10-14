// assembly.c
// library for polyominos

void OutputGrid(int *G, int ARR);
void FileOutputGrid(FILE *fp, int *G, int ARR);
void Push(int *data, int *neighbourlist, int *length);
void Pop(int *neighbourlist, int *length, int index);
int Add(int x, int y, int type, int ori, int *Grid, int *sides, int n, int *sidelist, int *neighbourlist, int *length, int allownd, int allowub, int ARR);
int mypow2(int n);
void Convert(int *p, int *q, int NTILE, int NBITCOL);
void Confound(int *s, int *d, int confoundtype, int LEN);
int Grow(int *sides, int n, int *Grid, int *size, int allownd, int allowub, int checksize, int display, int ARR);
  
// outputs a grid containing a polyomino structure to the screen
void OutputGrid(int *G, int ARR)
{
  int x, y;

  for(y = 0; y < ARR; y++)
    {
      for(x = 0; x < ARR; x++)
	printf("%c", (G[y*ARR+x] == -1 ? '.' : G[y*ARR+x]+48));
      printf("\n");
    }
}

// outputs grid to file
void FileOutputGrid(FILE *fp, int *G, int ARR)
{
  int x, y;

  for(y = 0; y < ARR; y++)
    {
      for(x = 0; x < ARR; x++)
	fprintf(fp, "%c", (G[y*ARR+x] == -1 ? '.' : G[y*ARR+x]+48));
      fprintf(fp, "\n");
    }
}

   
/***** this stuff is for basic polyomino assembly ******/

/*********/
// the array "sidelist" stores info on building block structures in the following form:
// type side type side type side ... n_1 times for each "1" patch
// type side type side ... n_2 times for each "2" side
// etc

// so sidelist for 1111 | 2000 would look like
// 0 0 0 1 0 2 0 3 0 0 0 0 0 0 0 0
// 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0

// the length of the array is then 4*2*n*(maxcols)

/**********/
// we'll also maintain an array called "neighbourlist" which stores all possible addition moves we can make
// each element of this array has 4 parts: x, y, block type, orientation 

// add a new possible move (data) to neighbourlist
void Push(int *data, int *neighbourlist, int *length)
{
  int i;
  for(i = 0; i < 4; i++)
    neighbourlist[4*(*length)+i] = data[i];
  (*length)++;
}

// remove a possible move from neighbourlist
void Pop(int *neighbourlist, int *length, int index)
{
  int i;
  (*length)--;
  for(i = index*4; i < (*length)*4; i++)
    neighbourlist[i] = neighbourlist[i+4];
} 

// this function adds a block to the assembly grid
// given the move description x, y, block type, orientation
// and updates neighbourlist appropriately
int Add(int x, int y, int type, int ori, int *Grid, int *sides, int n, int *sidelist, int *neighbourlist, int *length, int allownd, int allowub, int ARR)
{
  int i;
  int bond, partner;
  int data[4];
  int nondet;

  // add this tile to the grid
  Grid[y*ARR+x] = type;
  nondet = 0;
  // this compares our proposed move to all other current possible moves
  // to check for non-determinism and remove all possible moves in the
  // space we've chosen
  for(i = 0; i < *length; i++)
    {
      if(neighbourlist[4*i] == x && neighbourlist[4*i+1] == y)
	{
	  if(allownd == 0 && (neighbourlist[4*i+2] != type || neighbourlist[4*i+3] != ori))
	    {
	      // uncomment this for immediate nondet termination
	      //	      printf("NONDET FOUND: %i,%i : %i %i not %i %i\n", x, y, neighbourlist[4*i+2], neighbourlist[4*i+3], type, ori); 
	      return -1;
	    }
	  Pop(neighbourlist, length, i);
	  i--;
	}
    }

  // the subsequent four sections compute the new possible moves that we've gained on each side of the new block
  // and push these possible moves into neighbourlist
  if(x > 0)
    {
      bond = sides[4*type+(7-ori)%4];
      if(bond != 0 && Grid[y*ARR+x-1] == -1) 
	{
	  partner = (bond%2 == 0 ? bond-1 : bond+1);
	  for(i=0;i<sidelist[partner];i++)
	    {
	      data[0] = x-1; data[1] = y;
	      data[2] = sidelist[4*2*n*partner+2*i]; data[3] = (5-sidelist[4*2*n*partner+2*i+1])%4;
	      Push(data, neighbourlist, length);
	    }
	}
    }

  if(x < ARR-1)
    {
      bond = sides[4*type+(5-ori)%4];
      if(bond != 0 && Grid[y*ARR+x+1] == -1)
	{
	  partner = (bond%2 == 0 ? bond-1 : bond+1);
	  for(i=0;i<sidelist[partner];i++)
	    {
	      data[0] = x+1; data[1] = y;
	      data[2] = sidelist[4*2*n*partner+2*i]; data[3] = (3-sidelist[4*2*n*partner+2*i+1])%4;
	      Push(data, neighbourlist, length);
	    }
	}
    }
  if(y > 0)
    {
      bond = sides[4*type+(4-ori)%4];
      if(bond != 0 && Grid[(y-1)*ARR+x] == -1)
	{
	  partner = (bond%2 == 0 ? bond-1 : bond+1);
	  for(i=0;i<sidelist[partner];i++)
	    {
	      data[0] = x; data[1] = y-1;
	      data[2] = sidelist[4*2*n*partner+2*i]; data[3] = (6-sidelist[4*2*n*partner+2*i+1])%4;
	      Push(data, neighbourlist, length);
	    }
	}
    }
  if(y < ARR-1)
    {
      bond = sides[4*type+(6-ori)%4];
      if(bond != 0 && Grid[(y+1)*ARR+x] == -1) 
	{
	  partner = (bond%2 == 0 ? bond-1 : bond+1);
	  for(i=0;i<sidelist[partner];i++)
	    {
	      data[0] = x; data[1] = y+1;
	      data[2] = sidelist[4*2*n*partner+2*i]; data[3] = (4-sidelist[4*2*n*partner+2*i+1])%4;
	      Push(data, neighbourlist, length);
	    }
	}
    }
  if((x == 0 || x == ARR-1 || y == 0 || y == ARR-1) && allowub == 0)
    return -2;

  return 0;
}

// more efficient than calling pow
int mypow2(int n)
{
  int v = 1;
  int i;

  for(i = 0; i < n; i++)
    v *= 2;
  return v;
}

// convert a bitstring to an integer ruleset
void Convert(int *p, int *q, int NTILE, int NBITCOL)
{
  int i, j, k;

  for(i = 0; i < NTILE; i++)
    {
      for(j = 0; j < 4; j++)
	{
	  q[4*i+j] = 0;
	  for(k = 0; k < NBITCOL; k++)
	    {
	      q[4*i+j] += (p[NBITCOL*4*i+NBITCOL*j+k] == 1 ? mypow2(NBITCOL-1-k) : 0);
	    }
	}
    }
}

// a function to deliberately confuse the GP map
void Confound(int *s, int *d, int confoundtype, int LEN) 
{
  int i;
  int sum = 0;
  float rate;

  rate = confoundtype-1;
  
  if(confoundtype == 1)
    {
      // given a bitstring genome, we first count the number of 1s, then bitshift all elements by that sum
      // mutations thus have structural as well as local effects
      for(i = 0; i < LEN; i++)
	sum += s[i];
      for(i = 0; i < LEN; i++)
	d[(i+sum)%LEN] = s[i];
    }
  else if(confoundtype != 0)
    {
      // simply add random perturbations to the bitstring at some rate
      for(i = 0; i < LEN; i++) 
	d[i] = (drand48() < rate/LEN ? 1-s[i] : s[i]);
    }
}

// grow a polyomino
int Grow(int *sides, int n, int *Grid, int *size, int allownd, int allowub, int checksize, int display, int ARR)
{
  int i, j, k;
  int max = 0;
  int *sidelist;
  int *neighbourlist;
  int length;
  int ref;
  int x, y, type, ori;
  int *done;
  int run, last;

  for(i = 0; i < n*4; i++)
    {
      if(sides[i] > max) max = sides[i];
    }
  
  max++;

  // allocate memory for structures storing the sides of placed blocks and the neighbour list of possible moves we can make at any given time
  sidelist = (int*)malloc(sizeof(int)*4*n*max*2);
  neighbourlist = (int*)calloc(4*ARR*ARR*n, sizeof(int));

  for(i = 0; i < 4*n*max*2; i++)
    sidelist[i] = (i < max ? 0 : -1);

  // horrible looking piece of code to populate list of sides as described above
  // i.e. the (tile,side) coordinates for every instance of colour 1, colour 2, ...
  for(i = 0; i < n; i++)
    {
      for(j = 0; j < 4; j++)
	{
	  if(sides[4*i+j] != 0)
	    {
	      for(k = 0;;k+=2)
		{
		  if(sidelist[sides[4*i+j]*4*2*n+k] == -1)
		    {
		      sidelist[sides[4*i+j]*4*2*n+k] = i;
		      sidelist[sides[4*i+j]*4*2*n+k+1] = j;
		      sidelist[sides[4*i+j]]++;
		      break;
		    }
		}
	    }
	}
    }

  // simulate as many assembly runs as we are using to check determinism
  for(run = 0; run <= checksize; run++)
    {
      // initialise assembly grid with empty space 
      for(i = 0; i < ARR*ARR; i++)
	Grid[i] = -1;

      length = 0;

      // add the first tile; this will update the neighbourlist
      Add(ARR/2, ARR/2, 0, 0, Grid, sides, n, sidelist, neighbourlist, &length, allownd, allowub, ARR);
      *size = 1;

      // as long as the length of neighbourlist is nonzero, pick a random element and use this as the next addition move
      while(length!=0){
	ref = (int)(drand48()*length);
	x = neighbourlist[ref*4]; y = neighbourlist[ref*4+1];
	type = neighbourlist[ref*4+2]; ori = neighbourlist[ref*4+3];
	// add the chosen tile and update neighbourlist appropriately
	ref = Add(x, y, type, ori, Grid, sides, n, sidelist, neighbourlist, &length, allownd, allowub, ARR);
	if(ref < 0) { free(sidelist); free(neighbourlist); return ref; }
	(*size)++;
      }

      // catch nondeterminism
      if(checksize && run > 0 && *size != last)
	{free(sidelist); free(neighbourlist); return -3;}

      last = *size;
    }

  if(display)
    OutputGrid(Grid, ARR);

  free(sidelist); free(neighbourlist);
  return 0;
}


