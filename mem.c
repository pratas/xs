// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
// -                                                                         -
// -           U n i v e r s i t y   o f   A v e i r o   2 0 1 3             -
// -                                                                         -
// -         XS is under GPL,  http://www.gnu.org/licenses/gpl.txt           -
// -                                                                         -
// -                       XS: a FASTQ read simulator                        -
// -                                                                         -
// -              Contacts:  pratas@ua.pt, ap@ua.pt, jmr@ua.pt               -
// -                                                                         -
// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
//////////////////////////////////////////////////////////////////////////////

#include "mem.h"
#include <stdlib.h>

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

static uint64_t totalMemory = 0;

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

void *Malloc(size_t size)
  {
  void *pointer = malloc(size);

  if(pointer == NULL)
    {
    fprintf(stderr, "Error allocating %zu bytes.\n", size);
    exit(1);
    }

  totalMemory += size;

  return pointer;
  }

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

void *Calloc(size_t nmemb, size_t size)
  {
  void *pointer = calloc(nmemb, size);

  if(pointer == NULL)
    {
    fprintf(stderr, "Error allocating %zu bytes.\n", size);
    exit(1);
    }

  totalMemory += nmemb * size;

  return pointer;
  }

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

void *Realloc(void *ptr, size_t size, size_t additionalSize)
  {
  void *pointer = realloc(ptr, size);
  
  if(pointer == NULL) 
    {
    fprintf(stderr, "Error allocating %zu bytes.\n", size);
    exit(1);
    }

  totalMemory += additionalSize;

  return pointer;
  }

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

uint64_t TotalMemory()
  {
  return totalMemory;
  }

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

