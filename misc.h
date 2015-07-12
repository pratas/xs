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

#ifndef MISC_H_INCLUDED
#define MISC_H_INCLUDED

#include "defs.h"

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

typedef struct
  {
  uint32_t      number;                                   // Number of repeats
  uint32_t      minSize;                                // Minimum repeat size
  uint32_t      maxSize;                                // Maximum repeat size
  double        mutation;                        // Average mutation by repeat
  uint8_t       rcr;                      // Reverse complement repeats on|off
  }
Repeats;

typedef struct
  {
  int32_t       size;
  int32_t       *position;
  int32_t       deepestPosition;
  }
T;

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

uint32_t        complement(uint32_t base);
double          gaussRand(double mean, double sDev);
char            numbToChar(uint32_t n, uint32_t type);
FILE            *Fopen(const char *path, const char *mode);
uint32_t        max(uint32_t a, uint32_t b);
uint32_t        min(uint32_t a, uint32_t b);
void            sortIntVect(uint32_t vector[], uint32_t size);
uint32_t        toABISOLID(uint32_t previous, uint32_t current);
void            IncreaseTemplateStorage(T *t, int additionalSize);

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

#endif

