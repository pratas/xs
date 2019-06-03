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

#ifndef DEFS_H_INCLUDED
#define DEFS_H_INCLUDED

#include <stdint.h>

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

#define VERSION                              "2.1.0"
#define toBase36(x)                          ((x)<26 ? (x)+'A' : ((x)-26)+'0')
#define ABS(a)                               ((a)>0?(a):-(a))
#define MIN(a,b)                             ((a)<(b)?(a):(b))
#define MAX(a,b)                             ((a)>(b)?(a):(b))
#define STOP_RUNNING                         1
#define MAX_PERIOD_RAND                      4294967296
#define DEFAULT_BSIZE                        16384
#define DEFAULT_FSIZE                        1024
#define DEFAULT_LINE_SIZE                    70
#define DEFAULT_NREADS                       1000
#define DEFAULT_SEQUENCING_TYPE              1
#define DEFAULT_QS_NSYMBOLS                  32
#define DEFAULT_NUCLEOTIDE_FREQUENCY         0.225
#define DEFAULT_N_FREQUENCY                  0.100 
#define DEFAULT_REPEATS_NUMBER               0
#define DEFAULT_REPEATS_MAX_SIZE             100
#define DEFAULT_REPEATS_MIN_SIZE             3000
#define DEFAULT_MUTATION_RATE                0.01
#define DEFAULT_REPEATS_VERBOSE_SORTING      1000
#define DEFAULT_READS_VERBOSE_SORTING        100000
#define DEFAULT_COLOR_SPACE_SEPARATOR        5
#define DEFAULT_COLOR_SPACE_MAX_VALUE        10
#define DEFAULT_MINIMUM_LINESIZE             30
#define DEFAULT_MAXIMUM_LINESIZE             40
#define UNIFORM                              1
#define NORMAL                               2
#define MAX_COUNTER_FLOW_CELL_TILE           20
#define MAX_COUNTER_ID                       1100
#define MAX_CUMULATIVE                       1.0
#define MAX_FLOWCELLTILE                     100
#define MAX_XCOR                             2000
#define MAX_X                                2020
#define MAX_Y                                2020
#define MAX_YCOR                             2000
#define MIN_COUNTER_ID                       500
#define MIN_XCOR                             1
#define MIN_YCOR                             1
#define NSYMBOLS_DNA                         5
#define OUTBOUND_RAND_VALUE                  1.0
#define ITORRENT                             4
#define ABISOLID                             3
#define ILLUMINA                             2
#define ROCHE454                             1
#define DEFAULT_ABISOLID_ASCII_INIT          33
#define DEFAULT_ABISOLID_ASCII_END           126
#define DEFAULT_ILLUMINA_ASCII_INIT          33
#define DEFAULT_ILLUMINA_ASCII_END           126
#define DEFAULT_ROCHE454_ASCII_INIT          33
#define DEFAULT_ROCHE454_ASCII_END           126
#define DEFAULT_ITORRENT_ASCII_INIT          33
#define DEFAULT_ITORRENT_ASCII_END           126
#define IT_MIN_XCOR                          0
#define IT_MAX_XCOR                          15456
#define IT_MIN_YCOR                          0
#define IT_MAX_YCOR                          10656

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

#endif
