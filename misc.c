// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
// -                                                                         -
// -           U n i v e r s i t y   o f   A v e i r o   2 0 1 3             -
// -                                                                         -
// -         XS is under GPL,  http://www.gnu.org/licenses/gpl.txt           -
// -                                                                         -
// -                       XS: a FASTQ read simulator                        -
// -                                                                         -
// -             Contacts:  pratas@ua.pt , ap@ua.pt , jmr@ua.pt              -
// -                                                                         -
// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
//////////////////////////////////////////////////////////////////////////////

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <float.h>
#include "mem.h"
#include "misc.h"

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

uint32_t complement(uint32_t base)
  {
  switch(base)
    {
    case 0:  return 3;  break;
    case 1:  return 2;  break;
    case 2:  return 1;  break;
    case 3:  return 0;  break;
    case 4:  return 4;  break;
    default:            break;
    }

  fprintf(stderr, "Error assigning complement: %u\n", base);
  return 0;
  }

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

char numbToChar(uint32_t n, uint32_t type)
  {
  if(type == ABISOLID)
    {
    switch(n)
      {
      case 0:  return '0';  break;
      case 1:  return '1';  break;
      case 2:  return '2';  break;
      case 3:  return '3';  break;
      case 4:  return '.';  break;
      default:              break;
      }
    }
  else
    {
    switch(n)
      {
      case 0:  return 'A';  break;
      case 1:  return 'C';  break;
      case 2:  return 'G';  break;
      case 3:  return 'T';  break;
      case 4:  return 'N';  break;
      default:              break;
      }
    }

  fprintf(stderr, "Error assigning character: %u\n", n);
  return '.';
  }

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

FILE *Fopen(const char *path, const char *mode)
  {
  FILE *file = fopen(path, mode);

  if(file == NULL)
    {
    fprintf(stderr, "Error opening file %s with mode %s.\n", path, mode);
    exit(1);
    }

  return file;
  }

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

uint32_t max(uint32_t a, uint32_t b)
  {
  return a > b ? a : b; 
  }

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

uint32_t min(uint32_t a, uint32_t b)
  { 
  return a < b ? a : b; 
  }

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

void swapVect(uint32_t vector[], uint32_t i)
  { 
  uint32_t tmp = 0; 
  tmp          = vector[i];
  vector[i]    = vector[i+1];
  vector[i+1]  = tmp;
  }

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
 
void sortIntVect(unsigned int vector[], unsigned int size)
  {
  uint32_t i, change;
 
  do
    {
    size--;
    change = 0;
    for(i = 0; i < size; i++)
      if(vector[i] > vector[i+1])
        {
        swapVect(vector, i);
        change = 1;
        }
    }
  while(change);
  }

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

double gaussRand(double mean, double sDev)
  {
  double   v = 0.0, value = 0.0;
  uint32_t div;

  if(!(DBL_EPSILON < sDev))
    {
    fprintf(stderr, "Standard desviation is too small\n");
    return mean;
    }

  v     = (double) rand() / (RAND_MAX + MAX_CUMULATIVE);  
  div   = (v < 0.5) ? 0 : 1;
  v     = fabs(v - 1.0);                        
  v     = sDev * sqrt(-2.0 * log(v));
  value = div ? mean + v : mean - v;
  if(0.0 <= value && value <= 1.0) 
    return value; 

  return OUTBOUND_RAND_VALUE;  
  }

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

void IncreaseTemplateStorage(T *t, int32_t additionalSize)
  {
//  t->position = (int32_t *) (t->size == 0 ? Malloc(additionalSize * 
//  sizeof(int32_t)) : Realloc(t->position, (additionalSize + t->size) * 
//  sizeof(int32_t), additionalSize * sizeof(int32_t)));

  if(t->size == 0)
    t->position = (int32_t *) Malloc(additionalSize * sizeof(int32_t));
  else
    t->position = (int32_t *) Realloc(t->position, (additionalSize + 
      t->size) * sizeof(int32_t), additionalSize * sizeof(int32_t));
  }

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

uint32_t toABISOLID(uint32_t previous, uint32_t current)
  {
  uint32_t symbol;

  //TODO: rewrite with switches
  if((previous == 0 && current == 0) || (previous == 1 && 
    current == 1) || (previous == 2 && current == 2) || 
    (previous == 3 && current == 3))
    symbol = 0;
  else if((previous == 0 && current == 1) || (previous == 1 && 
    current == 0) || (previous == 2 && current == 3) || 
    (previous == 3 && current == 2))
    symbol = 1;
  else if((previous == 2 && current == 0) || (previous == 0 && 
    current == 2) || (previous == 1 && current == 3) || 
    (previous == 3 && current == 1))
    symbol = 2;
  else if((previous == 3 && current == 0) || (previous == 0 && 
    current == 3) || (previous == 1 && current == 2) || 
    (previous == 2 && current == 1))
    symbol = 3;
  else
    symbol = 4;

  return symbol;
  }

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

