int SymmLibrary(int *g, int size, int *lib, int *sizelib, int *numlib, int ARR, int *newpheno);

// this function takes a polyomino structure and compares it to a library of structures
// if the structure is found in the library, we return the reference to its entry
// if not, it is added to the library, and the new reference is returned
// we compare all orientations of the structure as we're interested in free polyominoes
int SymmLibrary(int *g, int size, int *lib, int *sizelib, int *numlib, int ARR, int *newpheno)
{
  int T[ARR*ARR];
  int ttry, x, y;
  int success[8];
  int i, j;
  int fail;
  int symm;

  int mini = ARR, minj = ARR;
  int n;

  *newpheno = 0;
  // we're looping through the 8 square lattice symmetry operations
  for(ttry = 0; ttry < 8; ttry++)
    {
      for(x = 0; x < ARR; x++)
	{
	  for(y = 0; y < ARR; y++)
	    {
	      // pop the structure, transformed by the current operation, in array T
	      switch(ttry)
		{
		case 0: T[x+ARR*y] = g[x+ARR*y]; break;
		case 1: T[x+ARR*y] = g[(ARR-1-x)+ARR*y]; break;
		case 2: T[x+ARR*y] = g[ARR-1-x+ARR*(ARR-1-y)]; break;
		case 3: T[x+ARR*y] = g[x+ARR*(ARR-1-y)]; break;
		case 4: T[x+ARR*y] = g[y+ARR*x]; break;
		case 5: T[x+ARR*y] = g[ARR-1-y+ARR*x]; break;
		case 6: T[x+ARR*y] = g[ARR-1-y+ARR*(ARR-1-x)]; break;
		case 7: T[x+ARR*y] = g[y+ARR*(ARR-1-x)]; break;
		}
	    }
	}
      mini = minj = ARR;
      // locate the top left corner of the transformed structure
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

      // loop through all extant structures looking for a match
      for(n = 0; n < *numlib; n++)
	{
	  // first check for a size match: skip if we're the wrong size
	  if(size == sizelib[n])
	    {
	      fail = 0;
	      // compare structures pixel-by-pixel
	      for(i = 0; i < ARR && fail == 0; i++)
		{
		  for(j = 0; j < ARR && fail == 0; j++)
		    {
		      x = i+mini; y = j+minj;
		      if(x < ARR && y < ARR)
			{
			  if((T[y*ARR+x] == -1 && lib[n*ARR*ARR+j*ARR+i] != -1) || (T[y*ARR+x] != -1 && lib[n*ARR*ARR+j*ARR+i] == -1))
			    fail = 1;
			}
		    }
		}
	      // if we've found a match, just return the reference
	      if(fail == 0) return n;
	    }
	}

    }

  // if we made it here, we've got a new structure
  *newpheno = 1;

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

  // add it to the next open slot in the library
  for(i = 0; i < ARR; i++)
    {
      for(j = 0; j < ARR; j++)
	{
	  x = i+mini; y = j+minj;
	  if(x < ARR && y < ARR)
	    lib[(*numlib)*ARR*ARR+j*ARR+i] = g[y*ARR+x];
	  else lib[(*numlib)*ARR*ARR+j*ARR+i] = -1;
	}
    }

  // record size and return new reference
  sizelib[*numlib] = size;
  (*numlib)++;

  return (*numlib)-1;
}

