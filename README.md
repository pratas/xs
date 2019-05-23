<p align="center"><img src="/xs.png" 
alt="XS" width="200" height="200" border="0" /></p>
XS is a skilled FASTQ read simulation tool, flexible, portable (does not need a reference sequence) and tunable in terms of sequence complexity. XS handles Ion Torrent, Roche-454, Illumina and ABI-SOLiD simulation sequencing types. It has several running modes, depending on the time and memory available, and is aimed at testing computing infrastructures, namely cloud computing of large-scale projects, and testing FASTQ compression algorithms. Moreover, XS offers the possibility of simulating the three main FASTQ components individually (headers, DNA sequences and quality-scores). Quality-scores can be simulated using uniform and Gaussian distributions.

## INSTALLATION ##

<pre>
git clone https://github.com/pratas/XS.git
cd XS/
make
</pre>

### Run XS

Run XS:

<pre>
./XS [parameters] [outFile]
</pre>

## PARAMETERS

To see the possible options type
<pre>
./XS -h
</pre>

For additional help:

```
                                  PARAMETERS
------------------------------------------------------------------------------

System options:

 -h                       give this help -> show possible parameters;

 -v                       verbose mode -> show more information;

------------------------------------------------------------------------------

Main FASTQ options:

 -t  <sequencingType>     type: 1=Roche-454, 2=Illumina, 3=ABI SOLiD, 4=Ion 
                          Torrent -> the four main types of sequencing 
                          supported by XS, the default is option 1;

 -hf <headerFormat>       header format: 1=Length appendix, 2=Pair End -> Uses
                          format with pair ends or with appendix reporting the
                          size of the read, such as: "length=80" (only for 
                          sequencing types 2 and 3);

 -i  n=<instrumentName>   the unique instrument name (use n= before name) -> 
                          this name appears in the read;

 -o                       use the same header in third line of the read ->
                          some FASTQ files have an optional header which is a 
                          copy of the real header;

 -ls <lineSize>           static line (bases/quality scores) size -> the size
                          of each DNA base and quality scores size. This 
                          option assume static size;

 -ld <minSize>:<maxSize>  dynamic line (bases/quality scores) size -> This 
                          sets a dynamic line size ranging sizes between a 
                          minimum and a maximum, such as "-ls 35:80";

 -n  <numberOfReads>      number of reads per file -> example: "-n 10000";

------------------------------------------------------------------------------

DNA options:

 -f  <A>,<C>,<G>,<T>,<N>  symbols frequency -> this sets the nucleotide 
                          distribution. A for Adenine, C for Cytosine, G for 
                          Guanine, T for Thymine and N for raw (any symbol)
                          Example for the overall human DNA:
                          "-f 0.29,0.19,0.19,0.29,0.04";

 -rn <numberOfRepeats>    repeats: number (default: 0) -> number of copies
                          (read article for more detailed explanation);

 -ri <repeatsMinSize>     repeats: minimum size -> minimum size of the 
                          repeats, such as: "-ri 300" (average minimum size of
                          transposable elements, TEs);

 -ra <repeatsMaxSize>     repeats: maximum size -> maximum size of the 
                          repeats. Example; "-ra 3000" (average maximum size 
                          of transposable elements, TEs);
                    
 -rm <mutationRate>       repeats: mutation frequency -> mutation rate.
                          Example: "-rm 0.01" mean that 1 base in 100 have an
                          uniform chance of be mutated;

 -rr                      repeats: use reverse complement repeats (also known 
                          as inverted repeats). Use: "-rr". For example: the 
                          string "ACGTA" will be reverse complemented to 
                          "TACGT";

------------------------------------------------------------------------------

Quality scores options:

 -qt <assignmentType>     quality scores distribution: 1=uniform, 2=gaussian. 
                          There are two possible quality scores destributions:
                          uniform (each symbol has got equal probability) and
                          normal (gaussian). Example: "-qt 1" will use an 
                          uniform distribution;

 -qf <statsFile>          load file: mean, standard deviation (when: -qt 2)
                          For a normal distribution there is the option of 
                          load a file with mean and standard deviation:
                          "-qf FILE", vi FILE:
                          0.4:0.10
                          0.3:0.15
                          ...
                          0.1:0.01
                          If this option is not used an array with both will
                          be used (mean and stdev in init.h); 

 -qc <template>           custom template ascii alphabet; It is possible to 
                          use custom quality scores, such as: 
                          "-qc 33:36,55,57:59" wil use the ascii values:
                          33,34,35,36,55,57,58,59;

------------------------------------------------------------------------------

Filtering options:

 -eh                      excludes the use of headers from output -> The 
                          output FASTQ file will omit this data source;

 -eo                      excludes the use of optional headers (+) from output

 -ed                      excludes the use of DNA bases from output -> The
                          output FASTQ file will omit this data source;

 -edb                     excludes '\n' when DNA bases line size is reached ->
                          this is normally used to facilitate the output
                          analysis, leaving the necessity of post-processing;

 -es                      excludes the use of quality scores from output ->
                          The output FASTQ file will omit this data source;

------------------------------------------------------------------------------

Stochastic options:

 -s  <seed>               generation seed -> a fixed seed with fixed 
                          parameters will generate the same output data.
                          However, a different seed with the same parameters
                          will generate different output data, but with the
                          same distribution;

------------------------------------------------------------------------------

<genFile>                 simulated output file -> FASTQ simulated file;

------------------------------------------------------------------------------

Common usage:

 ./XS -v -t 1 -i n=MySeq -ld 30:80 -n 20000 -qt=1 -qc 33,36,39:43 File
 (Common FASTQ file simulation)

 ./XS -v -ls 100 -n 10000 -eh -eo -es -edb -f 0.3,0.2,0.2,0.3,0.0 -rn 50 
-ri 300 -ra 3000 -rm 0.1 File
 (Transposable Elements simulation example)
```

## CITATION ##

On using this software please cite:

Pratas, D., Pinho, A. J., & Rodrigues, J. M. R. (2014). XS: a FASTQ read simulator. BMC research notes, 7(1), 40.

## ISSUES ##

For any issue let us know at [issues link](https://github.com/pratas/XS/issues).

## LICENSE ##

GPL v3.

For more information:
<pre>http://www.gnu.org/licenses/gpl-3.0.html</pre>
                                                    

