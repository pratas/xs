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
#include <ctype.h>
#include <math.h>
#include <time.h>
#include "defs.h"
#include "misc.h"
#include "msg.h"
#include "mem.h"
#include "init.h"
#include "s454.h"

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

void ParseListOfPositions(char *str, T *t)
  {
  int n, low, high;
  t->size = 0;
  while(1)
    {
    switch(sscanf(str, "%d:%d", &low, &high))
      {
      case 1: 
      if(low == 0) return;
      if(low < 0) exit(ErrorM(2));
      IncreaseTemplateStorage(t, 1);
      t->position[t->size++] = low;
      break;

      case 2: 
      if(low <= 0) exit(ErrorM(2));
      if(high - low + 1 > 0)
        {
        IncreaseTemplateStorage(t, high - low + 1);
        for(n = 0 ; n < high - low + 1 ; ++n)
          t->position[t->size++] = low + n;
        }
      break;

      default: return;
      }

    while(*str != '\0' && *str++ != ',')
      ;
    }
  }

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

int32_t main(int32_t argc, char *argv[])
  {
  FILE        *OUT = NULL, *STATS = NULL;
  uint32_t    n, m, k, z, numberOfReads, lineSize, sequencingType, flowCell, 
              flowCellTile, xCor, yCor, counterId, counterFlowCellTile, id, 
              indexBase, repeatsMinSize, repeatsMaxSize, numberOfRepeats, 
              rSize, *imagePositions, *objectPositions, *orderPos, maxVect,
              *memoryBases, minVect, *tmp1, *tmp2, usedRepeats, nEntries, 
              asciiOut[4][2], scoresSize, previousCS, minLine, maxLine, 
              minLineSize, maxLineSize, currentLineSize, headerFormat;
  uint8_t     help, verbose, optionalHeader, templateOn, edb;
  uint64_t    randSeed;
  char        instrumentName[DEFAULT_FSIZE], usageInstrumentName, symbol, 
              fNameGen[DEFAULT_FSIZE], buffer[DEFAULT_BSIZE], *bases, 
              *scores, rcr, doingRepeat, rcrActive, QStype, setSeed,
              dynamicLineSize, printHeaders, printOptionalH, printScores, 
              printDNA, tmpString[DEFAULT_FSIZE], fNameStats[DEFAULT_FSIZE],
              xStr[10], yStr[10]; 
  double      cpuTimeUsed, freqs[NSYMBOLS_DNA], cumulative[NSYMBOLS_DNA], tmp,
              totalFrequency, mutationRate, *QSfreqs, *QScumulative, v1, v2;
  clock_t     tic, start;
  Repeats     repeats;
  T           *t;

  start               = clock();
  k                   =  0;
  n                   =  0;
  m                   =  0;
  z                   =  0;
  id                  =  1;
  v1                  =  0.0;
  v2                  =  0.0;
  xCor                =  1;
  yCor                =  1;
  edb                 =  'n';
  rcr                 =  'n';
  help                =  'n';
  symbol              =  0;
  rSize               =  100;
  minVect             =  1;
  maxVect             =  2;
  minLine             =  10;
  maxLine             =  20;
  nEntries            =  0;
  randSeed            =  0;
  flowCell            =  1;
  orderPos            =  0;
  setSeed             =  'n';
  verbose             =  'n';
  previousCS          =  0;
  printDNA            =  'y';
  memoryBases         =  0; 
  rcrActive           =  'n';
  usedRepeats         =  1;
  templateOn          =  'n';
  flowCellTile        =  1;
  headerFormat        =  1;
  QStype              =  UNIFORM;
  printScores         =  'y';
  doingRepeat         =  'n';
  printHeaders        =  'y';
  imagePositions      =  0; 
  objectPositions     =  0;
  currentLineSize     =  30;
  optionalHeader      =  'n';
  printOptionalH      =  'y';
  dynamicLineSize     =  'n';
  usageInstrumentName =  'n';
  lineSize            =  DEFAULT_LINE_SIZE;
  numberOfReads       =  DEFAULT_NREADS;
  mutationRate        =  DEFAULT_MUTATION_RATE;
  minLineSize         =  DEFAULT_MINIMUM_LINESIZE;
  maxLineSize         =  DEFAULT_MAXIMUM_LINESIZE;
  sequencingType      =  DEFAULT_SEQUENCING_TYPE;
  numberOfRepeats     =  DEFAULT_REPEATS_NUMBER;
  repeatsMinSize      =  DEFAULT_REPEATS_MAX_SIZE; 
  repeatsMaxSize      =  DEFAULT_REPEATS_MIN_SIZE;
  asciiOut[3][1]      =  DEFAULT_ITORRENT_ASCII_END;
  asciiOut[2][1]      =  DEFAULT_ABISOLID_ASCII_END;
  asciiOut[1][1]      =  DEFAULT_ILLUMINA_ASCII_END;
  asciiOut[0][1]      =  DEFAULT_ROCHE454_ASCII_END;
  asciiOut[0][0]      =  DEFAULT_ROCHE454_ASCII_INIT;
  asciiOut[1][0]      =  DEFAULT_ILLUMINA_ASCII_INIT;
  asciiOut[2][0]      =  DEFAULT_ABISOLID_ASCII_INIT;
  asciiOut[3][0]      =  DEFAULT_ITORRENT_ASCII_INIT;

  for(k = 0 ; k < NSYMBOLS_DNA-1 ; k++)
    freqs[k] = DEFAULT_NUCLEOTIDE_FREQUENCY;

  freqs[NSYMBOLS_DNA-1] = DEFAULT_N_FREQUENCY;

  for(n = 1 ; n < argc ; n++)
    if(!strcmp("-h", argv[n]))
      { help = 'y'; break; }

  if(help == 'y')
    exit(ErrorM(0));
  else if(argc < 2)
    exit(ErrorM(1));

  for(n = 1 ; n < argc ; n++)
    if(!strcmp("-v", argv[n]))
      { verbose = 'y'; break;  }

  for(n = 1 ; n < argc ; n++)
    if(!strcmp("-t", argv[n]))
      { sequencingType = atoi(argv[n+1]); break; }

  for(n = 1 ; n < argc ; n++)
    if(!strcmp("-hf", argv[n]))
      { headerFormat = atoi(argv[n+1]); break; }

  if(sequencingType < 1 || sequencingType > 4)
    exit(ErrorM(4));

  for(n = 1 ; n < argc ; n++)
    if(!strcmp("-i", argv[n]))
      { 
      if(sscanf(argv[n+1], "n=%s", tmpString) == 1)
        usageInstrumentName = 'y'; 
      break;
      }

  if(usageInstrumentName == 'y')
    strcpy(instrumentName, tmpString); 
  else
    switch(sequencingType)
      {
      case ABISOLID: 
      strcpy(instrumentName, "VAB_BARB_20130115_2_Broad_3b"); 
      break;
      case ILLUMINA: 
      strcpy(instrumentName, "HWI-EAS100R");
      break;
      case ROCHE454: 
      strcpy(instrumentName, "E70AHBIO");
      break;
      case ITORRENT: 
      strcpy(instrumentName, "QWRK0");
      break;
      default:       
      strcpy(instrumentName, "XS-UA");
      }
  
  for(n = 1 ; n < argc ; n++)
    if(!strcmp("-o", argv[n]))
      { optionalHeader = 'y'; break; }

  for(n = 1 ; n < argc ; n++)
    if(!strcmp("-ls", argv[n]))
      { lineSize = atoi(argv[n+1]); break; }

  for(n = 1 ; n < argc ; n++)
    if(!strcmp("-ld", argv[n]))
      { 
      if(sscanf(argv[n+1], "%u:%u", &minLine, &maxLine) == 2)
        {
        minLineSize     = minLine;
        maxLineSize     = maxLine;
        dynamicLineSize = 'y';
        }
      break; 
      }

  for(n = 1 ; n < argc ; n++)
    if(!strcmp("-n", argv[n]))
      { numberOfReads = atoi(argv[n+1]); break; }

  for(n = 1 ; n < argc ; n++)
    if(!strcmp("-f", argv[n]))
      {
      if(sscanf(argv[n+1], "%lf,%lf,%lf,%lf,%lf", &freqs[0], 
        &freqs[1], &freqs[2], &freqs[3], &freqs[4]) != 5)
        {
        for(k = 0 ; k < NSYMBOLS_DNA-1 ; k++)
          freqs[k] = DEFAULT_NUCLEOTIDE_FREQUENCY;
        freqs[NSYMBOLS_DNA-1] = DEFAULT_N_FREQUENCY;
        }
      totalFrequency = 0.0;
      for(k = 0 ; k < NSYMBOLS_DNA ; k++)
        totalFrequency += freqs[k];

      for(k = 0 ; k < NSYMBOLS_DNA ; k++)
        freqs[k] /= totalFrequency;
      break;
      }

  for(n = 1 ; n < argc ; n++)
    if(!strcmp("-rn", argv[n]))
      { numberOfRepeats = atoi(argv[n+1]); break; }

  for(n = 1 ; n < argc ; n++)
    if(!strcmp("-ri", argv[n]))
      { repeatsMinSize  = atoi(argv[n+1]); break; }

  for(n = 1 ; n < argc ; n++)
    if(!strcmp("-ra", argv[n]))
      { repeatsMaxSize  = atoi(argv[n+1]); break; }

  for(n = 1 ; n < argc ; n++)
    if(!strcmp("-rm", argv[n]))
      { mutationRate    = atof(argv[n+1]); break; }

  for(n = 1 ; n < argc ; n++)
    if(!strcmp("-rr", argv[n]))
      { rcr = 'y'; break; }

  cumulative[0] = freqs[0];
  for(k = 1 ; k < NSYMBOLS_DNA-1 ; k++)
    cumulative[k] = cumulative[k-1] + freqs[k];
  cumulative[NSYMBOLS_DNA-1] = MAX_CUMULATIVE;

  for(n = 1 ; n < argc ; n++) 
    if(!strcmp("-s", argv[n]))
      { randSeed = atol(argv[n+1]); setSeed = 'y'; break; }

  for(n = 1 ; n < argc ; n++) 
    if(!strcmp("-qt", argv[n]))
      { QStype = atoi(argv[n+1]); break; }

  t = (T *) Calloc(1, sizeof(T));

  for(n = 1 ; n < argc ; n++) 
    if(!strcmp("-qc", argv[n]))
      { 
      ParseListOfPositions(argv[n+1], &t[0]);
      templateOn = 'y';
      break;
      }
  
  for(n = 1 ; n < argc ; n++) 
    if(!strcmp("-qf", argv[n]))
      { STATS = Fopen(strcpy(fNameStats, argv[n+1]), "r"); break; }

  for(n = 1 ; n < argc ; n++)
    if(!strcmp("-eh", argv[n]))
      { printHeaders = 'n'; break; }

  for(n = 1 ; n < argc ; n++)
    if(!strcmp("-eo", argv[n]))
      { printOptionalH = 'n'; break; }

  for(n = 1 ; n < argc ; n++)
    if(!strcmp("-ed", argv[n]))
      { printDNA = 'n'; break; }

  for(n = 1 ; n < argc ; n++)
    if(!strcmp("-edb", argv[n]))
      { edb = 'y'; break; }

  for(n = 1 ; n < argc ; n++)
    if(!strcmp("-es", argv[n]))
      { printScores = 'n'; break; }

  scoresSize = (templateOn == 'n' ? asciiOut[sequencingType-1][1] - 
  asciiOut[sequencingType-1][0] + 1 : t[0].size);
  scores = (char *) Malloc((scoresSize + 1) * sizeof(char));
  for(n = 0 ; n < scoresSize ; ++n)
    scores[n] = (templateOn == 'n' ? asciiOut[sequencingType-1][0] + n : 
    t[0].position[n]);

  nEntries = 0;
  if(STATS)
    {
    while(nEntries < scoresSize && fscanf(STATS, "%lf:%lf", &v1, &v2) == 2)
      {
      mean[nEntries] = v1;
      sDev[nEntries] = v2;
      nEntries++;
      }
    if(nEntries != scoresSize - 1)
      exit(ErrorM(5));
    }
    
  if(setSeed == 'n')
    randSeed = time(NULL);
  srand(randSeed);

  repeats.number   = (printDNA == 'y' ? numberOfRepeats : 0);
  repeats.minSize  = repeatsMinSize;
  repeats.maxSize  = repeatsMaxSize;
  repeats.mutation = mutationRate;
  repeats.rcr      = rcr;

  if(repeatsMinSize >= repeatsMaxSize)
    exit(ErrorM(6));

  if(verbose == 'y')
    {
    printf("Sequencing type: ");
    switch(sequencingType)
      {
      case ABISOLID: printf("ABI SOLiD (Working on color space)\n"); break;
      case ILLUMINA: printf("Illumina\n");                           break;
      case ROCHE454: printf("Roche-454\n");                          break;
      case ITORRENT: printf("Ion Torrent\n");                        break;
      default:       printf("XS-UA\n");                              break;
      }

    if(printHeaders == 'y' && optionalHeader == 'y') 
      printf("Using optional copy of header\n");
    dynamicLineSize == 'n' ? printf("Static line size: %u\n", lineSize) : 
      printf("Dynamic line size, min: %u, max: %u\n", minLineSize, maxLineSize);      
    printf("Number of reads: %u\n", numberOfReads);

    if(printDNA == 'y')
      {
      printf("Frequencies: A=%.3f, C=%.3f, G=%.3f, T=%.3f, N=%.3f\n", 
      freqs[0], freqs[1], freqs[2], freqs[3], freqs[4]);
      printf("Number of repeats: %u   \n", repeats.number);
      printf("Repeats minimum size: %u\n", repeats.minSize);
      printf("Repeats maximum size: %u\n", repeats.maxSize);
      printf("Mutation frequency: %.3f\n", repeats.mutation);
      if(repeats.rcr == 'y')
        printf("Using reverse complement repeats\n");
      }

    if(printScores == 'y')
      printf("Scores template size: %d\n", scoresSize);

    printf("Using seed: %lu\n", randSeed);
    }

  OUT = Fopen(strcpy(fNameGen, argv[argc-1]), "w");

  counterFlowCellTile = rand() % MAX_COUNTER_FLOW_CELL_TILE + 1;

  if(dynamicLineSize == 'y')
    lineSize = maxLineSize;

  bases        = (char *)   Malloc((lineSize   + 1) * sizeof(char));
  QSfreqs      = (double *) Malloc((scoresSize + 1) * sizeof(double));
  QScumulative = (double *) Malloc((scoresSize + 1) * sizeof(double));
  
  //Uniform Mask for stochastic interval selection
  for(k = 0 ; k < scoresSize-1 ; k++)
    QSfreqs[k] = 1.0 / scoresSize;
  QSfreqs[scoresSize-1] = DEFAULT_N_FREQUENCY;
  QScumulative[0] = QSfreqs[0];
  for(k = 1 ; k < scoresSize-1 ; k++)
    QScumulative[k] = QScumulative[k-1] + QSfreqs[k];
  QScumulative[scoresSize-1] = MAX_CUMULATIVE;

  if(repeats.number > 0)
    {
    memoryBases     = (uint32_t *) Malloc(((lineSize * numberOfReads) + 1) * 
    sizeof(uint32_t));
    imagePositions  = (uint32_t *) Malloc(repeats.number * sizeof(uint32_t));
    objectPositions = (uint32_t *) Malloc(repeats.number * sizeof(uint32_t));
    orderPos        = (uint32_t *) Malloc(repeats.number * sizeof(uint32_t));
    tmp1            = (uint32_t *) Malloc(repeats.number * sizeof(uint32_t));
    tmp2            = (uint32_t *) Malloc(repeats.number * sizeof(uint32_t));
  
    for(n = 0 ; n < repeats.number ; n++)
      {
      imagePositions[n]  = rand() % (lineSize * numberOfReads) + 1;
      objectPositions[n] = rand() % (lineSize * numberOfReads) + 1;
      }
  
    if(verbose == 'y' && repeats.number > DEFAULT_REPEATS_VERBOSE_SORTING && 
    numberOfReads > DEFAULT_READS_VERBOSE_SORTING)
      printf("Sorting repeats...\n");

    sortIntVect(objectPositions, repeats.number);
    for(n = 0 ; n < repeats.number ; n++)
      {
      tmp1[n] = max(imagePositions[n], objectPositions[n]);
      tmp2[n] = tmp1[n];
      }
    sortIntVect(tmp1, repeats.number);
    for(n = 0 ; n < repeats.number ; n++)
      {
      for(m = 0 ; m < repeats.number ; m++)
        if(tmp1[n] == tmp2[m])
          {
          orderPos[n] = m;
          break;
          }
      if(verbose == 'y' && repeats.number > DEFAULT_REPEATS_VERBOSE_SORTING && 
      numberOfReads > DEFAULT_READS_VERBOSE_SORTING && n % (repeats.number/100) == 0)
        fprintf(stderr, "%3d%%\r", (int) (n/(repeats.number/100)));
      }

    z           = 0;
    maxVect     = max(imagePositions[orderPos[z]], objectPositions[orderPos[z]]);
    minVect     = min(imagePositions[orderPos[z]], objectPositions[orderPos[z]]);
    usedRepeats = 1;
    z++;
    if(repeats.rcr == 'y')
      rcrActive = rand() % 1;

    if(verbose == 'y' && repeats.number > DEFAULT_REPEATS_VERBOSE_SORTING && 
    numberOfReads > DEFAULT_READS_VERBOSE_SORTING)
      printf("Done!\n");
    }

  printf("Starting simulation...\n");

  currentLineSize = lineSize;
  indexBase = 0;
  for(m = 0 ; m < numberOfReads ; ++m)                 //TODO: WRITE W/ BUFFER
    {
    if(dynamicLineSize == 'y')
      currentLineSize = rand() % (maxLineSize + 1 -minLineSize) + minLineSize;

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
// - [1] Header

    if(counterFlowCellTile-- < 1) 
      {
      flowCellTile++;
      counterFlowCellTile = rand() % MAX_COUNTER_FLOW_CELL_TILE + 1;
      }

    if(flowCellTile == MAX_FLOWCELLTILE)
      {
      flowCellTile = 1;
      flowCell++;
      counterFlowCellTile = rand() % MAX_COUNTER_FLOW_CELL_TILE + 1;
      }

    counterId = rand() % (MAX_COUNTER_ID+1-MIN_COUNTER_ID) + MIN_COUNTER_ID;
    id       += counterId;
    xCor      = rand() % (MAX_XCOR + 1 - MIN_XCOR) + MIN_XCOR;
    yCor      = rand() % (MAX_YCOR + 1 - MIN_YCOR) + MIN_YCOR;

    if(printHeaders == 'y')
      {
      switch(sequencingType)
        {
        case ABISOLID:
        headerFormat == 1 ? sprintf(buffer, "%s.%u %s_%u_%u_%u_%u length=%u", 
        fNameGen, id, instrumentName, flowCell, flowCellTile, xCor, yCor, 
        currentLineSize) : sprintf(buffer, "%s.%u %s_%u_%u_%u_%u/%u", 
        fNameGen, id, instrumentName, flowCell, flowCellTile, xCor, yCor, 
        rand() % 2 + 1);
        //@OUTFILENAME.2286 VAB_BARB_20080515_2_Broad_3b_150_1176_133_2030 length=35
        //@OUTFILENAME.2286 VAB_solid0019_20100107_QR1004_50x2_QR1000417_21_640/2
        break; 
        
        case ILLUMINA:
        //NCBI Sequence Read Archive 1
        //@OUTFILENAME.2286 HWUSI-EAS100R:1:1:210:308 length=30
        //Illumina headers 2
        //@OUTFILENAME.2286 HWUSI-EAS100R:6:73:941:1973/1
        headerFormat == 1 ? sprintf(buffer, "%s.%u %s:%u:%u:%u:%u length=%u", 
        fNameGen, id, instrumentName, flowCell, flowCellTile, xCor, yCor, 
        currentLineSize) : sprintf(buffer, "%s.%u %s:%u:%u:%u:%u/%u", 
        fNameGen, id, instrumentName, flowCell, flowCellTile, xCor, yCor, 
        rand() % 2 + 1);
        break;
    
        case ROCHE454:
        sprintf(buffer, "%s.%u %s length=%u", fNameGen, id, get454UAN(((rand()
        % MAX_X + 1) * 4096) + (rand() % MAX_Y + 1)), currentLineSize);
        break;

        case ITORRENT:
        xCor = rand() % (IT_MAX_XCOR + 1 - IT_MIN_XCOR) + IT_MIN_XCOR;
        yCor = rand() % (IT_MAX_YCOR + 1 - IT_MIN_YCOR) + IT_MIN_YCOR;

        if(xCor < 10)
          sprintf(xStr, "0000%d", xCor);
        else if(xCor < 100)
          sprintf(xStr, "000%d", xCor);
        else if(xCor < 1000)
          sprintf(xStr, "00%d", xCor);
        else if(xCor < 10000)
          sprintf(xStr, "0%d", xCor);
        else
          sprintf(xStr, "%d", xCor);

        if(yCor < 10)
          sprintf(yStr, "0000%d", yCor);
        else if(yCor < 100)
          sprintf(yStr, "000%d", yCor);
        else if(yCor < 1000)
          sprintf(yStr, "00%d", yCor);
        else if(yCor < 10000)
          sprintf(yStr, "0%d", yCor);
        else 
          sprintf(yStr, "%d", yCor);

        sprintf(buffer, "%s:%s:%s", instrumentName, xStr, yStr);
        break;
        }
      fprintf(OUT, "@%s\n", buffer);
      }

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
// - [2] DNA Bases

    if(sequencingType == ABISOLID)
      {
      if(DEFAULT_COLOR_SPACE_SEPARATOR > rand() % 
      DEFAULT_COLOR_SPACE_MAX_VALUE + 1)
        {
        fprintf(OUT, "T");
        previousCS = 3;
        }
      else
        {
        fprintf(OUT, "G");
        previousCS = 2;
        }
      }
 
    for(n = 0 ; n < lineSize ; n++)
      {
      if(indexBase == maxVect && repeats.number > 0)
        {
        rSize = rand() % (repeats.maxSize-repeats.minSize) + repeats.minSize;
        doingRepeat = 'y';
        }
    
      if(doingRepeat == 'y')
        {
        if(repeats.rcr == 'y' && rcrActive > repeats.number / 2)
          {
          if(rand() / (RAND_MAX + 1.0) < mutationRate)
            {
            k = rand() % NSYMBOLS_DNA;
            minVect--;
            }
          else
            {
            k = complement(memoryBases[rSize+minVect--]);
            }            
          }
        else
          {
          if(rand() / (RAND_MAX + 1.0) < mutationRate)
            {
            k = rand() % NSYMBOLS_DNA;
            minVect++;
            }
          else
            {
            k = memoryBases[minVect++];
            }
          }
    
        if(indexBase == maxVect + rSize)
          {
          usedRepeats++;
          doingRepeat = 'n';
        
          if(repeats.number > z)
            {
            do            
              {
              maxVect = max(imagePositions[orderPos[z]], 
                objectPositions[orderPos[z]]);
              minVect = min(imagePositions[orderPos[z]], 
                objectPositions[orderPos[z]]);
              z++;
              }
            while(maxVect <= indexBase && indexBase < lineSize * numberOfReads
            && z < repeats.number);

            if(repeats.rcr == 'y')
              rcrActive = rand() % repeats.number + 1;
            }
          }
        }
      else
        {
        tmp = rand() / (RAND_MAX + MAX_CUMULATIVE);
        for(k = 0 ; k < NSYMBOLS_DNA ; ++k)
          if(tmp <= cumulative[k])
            break;
        }

      symbol = (sequencingType == ABISOLID ? numbToChar(toABISOLID(previousCS, 
      k), sequencingType) : numbToChar(k, sequencingType));

      bases[n] = symbol;
      if(repeats.number > 0)
        memoryBases[indexBase] = k;
      indexBase++;

      if(sequencingType == ABISOLID)
        previousCS = k;

      if(printDNA == 'y' && n < currentLineSize)
        fprintf(OUT, "%c", symbol);
      }

    if(printDNA == 'y' && edb == 'n')
      fprintf(OUT, "\n");

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
// - [3] Optional headers

    if(printOptionalH == 'y')
      {
      fprintf(OUT, "+");
      if(optionalHeader == 'y')
        fprintf(OUT, "%s", buffer);
      fprintf(OUT, "\n");
      }

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
// - [4] Quality scores

    if(sequencingType == ABISOLID)
      fprintf(OUT, "%c", scores[scoresSize - 1]);

    for(n = 0 ; n < lineSize ; n++)
      {
      if(QStype == UNIFORM || QStype == NORMAL)
        {
        if((bases[n] == 'N' && sequencingType != ABISOLID) || (bases[n] == '.'
        && sequencingType == ABISOLID))
          k = 0;
        else
          {
          tmp = (QStype == NORMAL ? gaussRand(mean[n], sDev[n]) : rand() / 
          (RAND_MAX + 1.0));
          for(k = 0 ; k < scoresSize ; k++)
            if(tmp <= QScumulative[k])
              {
              k = scoresSize - k - 1;
              break;
              }
          }
        }
      if(printScores == 'y' && n < currentLineSize)
        fprintf(OUT, "%c", scores[k]);
      }

    if(printScores == 'y')
      fprintf(OUT, "\n");

    if(numberOfReads > 100 && m % (numberOfReads / 100) == 0)
      fprintf(stderr, "%3d%%\r", (int) (m / (numberOfReads / 100)));
    }

  tic = clock();
  cpuTimeUsed = ((double) (tic-start)) / CLOCKS_PER_SEC;

  printf("Done!\n");
  if(verbose == 'y')
    printf("Used %u repeats from %u\n", usedRepeats-1, repeats.number);
  printf("Spent %g s and %.1f MB.\n", cpuTimeUsed, TotalMemory() / 1024. / 
  1024.);

  if(STATS) fclose(STATS);
  fclose(OUT);

  return EXIT_SUCCESS;
  }

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

