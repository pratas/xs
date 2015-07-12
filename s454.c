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
#include <time.h>
#include "mem.h"
#include "s454.h"

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

void reverse(char *s) 
  {
  int32_t length, n;
  char    *begin, *end, tmp;

  length  = strlen(s);
  begin   = s;
  end     = s;
  for(n = 0 ; n < length - 1 ; ++n)
    end++;

  for(n = 0 ; n < length / 2 ; ++n) 
    {        
    tmp    = *end;
    *end   = *begin;
    *begin = tmp;
    begin++;
    end--;
    }
  }

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

uint32_t getThisDate454TimeStamp(Date d)
  {
  return ((d.year - 2000) * 60 * 60 * 24 * 32 * 13 ) +
          (d.month        * 60 * 60 * 24 * 32      ) +                
          (d.day          * 60 * 60 * 24           ) +
          (d.hour         * 60 * 60                ) +
          (d.minute       * 60                     ) +
          (d.second                                ) ;
  }

uint32_t get454TimeStamp(time_t t)
  {
  return ((localtime(&t)->tm_year - 100) * 60 * 60 * 24 * 32 * 13 ) +
          (localtime(&t)->tm_mon         * 60 * 60 * 24 * 32      ) +                
          (localtime(&t)->tm_mday        * 60 * 60 * 24           ) +
          (localtime(&t)->tm_hour        * 60 * 60                ) +
          (localtime(&t)->tm_min         * 60                     ) +
          (localtime(&t)->tm_sec                                  ) ;
  }

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

char *get454UAN(uint32_t  e)
  {
  time_t  now;
  Date    date;
  int32_t d, n, code, chval = 0;
  char    hash[7], hash2[6], *s, *uan;

  uan = (char *)Malloc(15 * sizeof(char));

  if((now = time(NULL)) == (time_t)-1)
    {
    fprintf(stderr, "Error getting time!");
    //FALLBACK TIME
    date.year   = 2004;
    date.month  = 9;
    date.day    = 22;
    date.hour   = 16;
    date.minute = 59;
    date.second = 10;
    d = getThisDate454TimeStamp(date);
    }
  else
    d = get454TimeStamp(now);  

  for(n = 0 ; n < 6 ; n++)
    {
    code    = d % 36;
    d      /= 36;
    hash[n] = toBase36(code);
    }
  hash[n] = '\0';
  reverse(hash);
  for(n = 0 ; n < 6 ; n++)
    uan[n] = hash[n];

  for(s = hash; *s; s++)
    {
    chval += (int) *s;
    chval %= 31;
    }

  uan[6] = toBase36(chval);
  uan[7] = '0';
  uan[8] = '1';

  for(n = 0 ; n < 5 ; n++)
    {
    code = e % 36;
    e /= 36;
    hash2[n] = toBase36(code);
    }
  hash2[n] = '\0';
  reverse(hash2);
  for(n = 0 ; n < 5 ; n++)
    uan[9+n] = hash2[n];
  uan[14]  = '\0';
    
  return uan;
  }

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

