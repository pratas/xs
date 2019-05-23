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

#include "msg.h"
#include "stdio.h"

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

uint32_t ErrorM(uint32_t id)
  {
  switch(id)
    {
    case 0:
    fprintf(stderr, "Usage: XS   [OPTION]... [FILE] \n\n");
    fprintf(stderr, "System options:\n");
    fprintf(stderr, " -h                       give this help\n");
    fprintf(stderr, " -V                       version number\n");
    fprintf(stderr, " -v                       verbose mode\n\n");
    fprintf(stderr, "Main FASTQ options:\n");
    fprintf(stderr, " -t  <sequencingType>     type: 1=Roche-454, 2=Illumina,"
    " 3=ABI SOLiD, 4=Ion Torrent\n");
    fprintf(stderr, " -hf <headerFormat>       header format: 1=Length append"
    "ix, 2=Pair End\n");
    fprintf(stderr, " -i  n=<instrumentName>   the unique instrument name (us"
    "e n= before name)\n");
    fprintf(stderr, " -o                       use the same header in third l"
    "ine of the read\n");
    fprintf(stderr, " -ls <lineSize>           static line (bases/quality sco"
    "res) size\n");
    fprintf(stderr, " -ld <minSize>:<maxSize>  dynamic line (bases/quality sc"
    "ores) size\n");
    fprintf(stderr, " -n  <numberOfReads>      number of reads per file\n\n");
    fprintf(stderr, "DNA options:\n");
    fprintf(stderr, " -f  <A>,<C>,<G>,<T>,<N>  symbols frequency\n");
    fprintf(stderr, " -rn <numberOfRepeats>    repeats: number (default: 0)"
    "\n");
    fprintf(stderr, " -ri <repeatsMinSize>     repeats: minimum size\n");
    fprintf(stderr, " -ra <repeatsMaxSize>     repeats: maximum size\n");
    fprintf(stderr, " -rm <mutationRate>       repeats: mutation frequency"
    "\n");
    fprintf(stderr, " -rr                      repeats: use reverse complemen"
    "t repeats\n\n");
    fprintf(stderr, "Quality scores options:\n");
    fprintf(stderr, " -qt <assignmentType>     quality scores distribution: 1"
    "=uniform, 2=gaussian\n");
    fprintf(stderr, " -qf <statsFile>          load file: mean, standard devi"
    "ation (when: -qt 2)\n");
    fprintf(stderr, " -qc <template>           custom template ascii alphabet"
    "\n\n");
    fprintf(stderr, "Filtering options:\n");
    fprintf(stderr, " -eh                      excludes the use of headers fr"
    "om output\n");
    fprintf(stderr, " -eo                      excludes the use of optional h"
    "eaders (+) from output\n");
    fprintf(stderr, " -ed                      excludes the use of DNA bases "
    "from output\n");
    fprintf(stderr, " -edb                     excludes '\\n' when DNA bases "
    "line size is reached\n");
    fprintf(stderr, " -es                      excludes the use of quality sc"
    "ores from output\n\n");
    fprintf(stderr, "Stochastic options:\n");
    fprintf(stderr, " -s  <seed>               generation seed\n\n");
    fprintf(stderr, "<genFile>                 simulated output file\n\n");
    fprintf(stderr, "Common usage:\n");
    fprintf(stderr, " ./XS -v -t 1 -i n=MySeq -ld 30:80 -n 20000 -qt=1 -qc 33"
    ",36,39:43 File\n");
    fprintf(stderr, " ./XS -v -ls 100 -n 10000 -eh -eo -es -edb -f 0.3,0.2,0."
    "2,0.3,0.0 -rn 50 -ri 300 -ra 3000 -rm 0.1 File\n\n");
    fprintf(stderr, "Report bugs to {pratas,ap,jmr}@ua.pt\n");
    break;

    case 1:
    fprintf(stderr, "For help, type: XS -h\n");
    break;

    case 2:
    fprintf(stderr, "Error: can't handle non-causal contexts\n");
    break;

    case 3:
    fprintf(stderr, "Error: can't handle type of data\n");
    break;

    case 4:
    fprintf(stderr, "Error: unknown sequencing type\n");
    break;

    case 5:
    fprintf(stderr, "Error: unmatching size of stats\n");
    break;

    case 6:
    fprintf(stderr, "Error: invalid minimum and maximum repeats size\n");
    break;

    default:
    fprintf(stderr, "Error: unknown error\n");
    break;
    }

  return STOP_RUNNING;    
  }

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

