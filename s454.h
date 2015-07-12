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

#ifndef S454_H_INCLUDED
#define S454_H_INCLUDED

#include "defs.h"

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

typedef struct
  {
  uint32_t year;
  uint32_t month;
  uint32_t day;
  uint32_t hour;
  uint32_t minute;
  uint32_t second;
  }
Date;

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

char *get454UAN(uint32_t e);

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

#endif



