<p align="center"><img src="/xs.png" 
alt="XS" width="200" height="200" border="0" /></p>
XS is a skilled FASTQ read simulation tool, flexible, portable (does not need a reference sequence) and tunable in terms of sequence complexity. XS handles Ion Torrent, Roche-454, Illumina and ABI-SOLiD simulation sequencing types. It has several running modes, depending on the time and memory available, and is aimed at testing computing infrastructures, namely cloud computing of large-scale projects, and testing FASTQ compression algorithms. Moreover, XS offers the possibility of simulating the three main FASTQ components individually (headers, DNA sequences and quality-scores).

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

## CITATION ##

On using this software/method please cite:

Pratas, D., Pinho, A. J., & Rodrigues, J. M. R. (2014). XS: a FASTQ read simulator. BMC research notes, 7(1), 40.

## ISSUES ##

For any issue let us know at [issues link](https://github.com/pratas/XS/issues).

## LICENSE ##

GPL v2.

For more information:
<pre>http://www.gnu.org/licenses/gpl-2.0.html</pre>

                                                    

