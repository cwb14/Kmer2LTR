#### Pipeline to estimate LTR-RT insertion time.   
(1) Use kmers to identify LTRs from a library of LTR-RTs.  
(2) MAFFT and trimal to clean alignments.  
(2) Force global alignment gapless LTRs.   
(3) Date with p-dist, JC69, and K2P.  

```
python Kmer2LTR/Kmer2LTR.py -h
usage: Kmer2LTR.py [-h] [-k] [-v] [-d DIST] [-l KMIN] [-U KMAX] [-u MUTATION_RATE] [-f STD_FACTOR] [-e EXTENSION] [-t TEMP_DIR] [-o OUTFILE]
                   [-p THREADS]
                   input_fasta

Process multi-seq LTR-RT FASTA to extract and align LTRs.

positional arguments:
  input_fasta       Path to multi-sequence LTR-RT FASTA file.

options:
  -h, --help        show this help message and exit
  -k                Keep temp directory after processing.
  -v                Verbose mode; print each command before executing.
  -d DIST           Minimum distance between kmer pairs (default: 80).
  -l KMIN           Minimum kmer length (default: 8).
  -U KMAX           Maximum kmer length (default: 12).
  -u MUTATION_RATE  Mutation rate μ (default: 3e-8).
  -f STD_FACTOR     Standard deviation factor for kmer filtering.
  -e EXTENSION      Extension length for LTR extraction (default: 65).
  -t TEMP_DIR       Temporary directory name (default: ./temp).
  -o OUTFILE        Output filename (default: ./LTRs.alns.results).
  -p THREADS        Number of parallel threads (default: 20).
  -D DOMAINS_TSV, --domains DOMAINS_TSV
                        Optional TSV (name\tLTR_len). If provided, skip LTR discovery and go stright to alignment
```


```
python Kmer2LTR/Kmer2LTR.py LTR-RT.fa
```

Creates `LTRs.alns.results`.
Output format:
```
<LTR-RT>  <LTR_LEN>  <ALN_LEN>  <substitions> <transitions>  <transversions>  <p-dist> <p-time> <JC69-dist> <JC69-time> <K2P-dist>  <K2P-time>
Gypsy1#LTR_Ty3	584  574	131	109	22	0.228223	3803717	0.272125	4535411	0.290682	4844705
Gypsy2#LTR_Ty3  260  260	55	47	8	0.211538	3525641	0.248518	4141964	0.264922	4415361
Gypsy3#LTR_Ty3  741  744	180	137	43	0.241935	4032258	0.292099	4868310	0.308338	5138959
```
`LTR_LEN` is whats discovered via kmer.  
`ALN_LEN` is the length of aligned bases.  
`LTR_LEN` may be less than or greater than `ALN_LEN`.  
`LTR_LEN` > `ALN_LEN` if gaps in LTR.   
`LTR_LEN` < `ALN_LEN` if extension expands the LTR boundry discovered by kmer.   
















# Developers note.
TESS-PrinTE makes `lib_clean.fa` with LTR length in the header.   
All LTRs are unmutated.   

PIPELINE OUTLINE: (1) Use kmers to identify LTR boundry from LTR-RTs. (2) Force global alignment of LTRs (gaps penalized to promote subsitution).  

For benchmarking the pipeline, I can:
```
python lib_mutator.py -i lib_clean.fa -o lib_clean.15mp.titv2.fa -mp 15 -TiTv 2 --seed 51
```
LTRs are 15% mutated. TiTv ratio = 2. 
Internal sequence is not mutated. 

This new kmer2ltr pipeline has issues with ':' and '/' in header. 
I will fix if the pipeline is adopted. 
But for now...
```
cat lib_clean.15mp.titv2.fa | sed 's/\//_/g' | sed 's/LTRlen:/LTRlen/g' > lib_clean.15mp.titv2.clean.fa
```

Run the pipeline.
```
nohup /usr/bin/time -v python Kmer2LTR.py -p 100 lib_clean.15mp.titv2.clean.fa &
```
Default parameters seem fairly optimized. 

How many does it match in full length?
```
awk -F'\t' '{ split($1, a, /[a-zA-Z_~#-]+/); n=a[length(a)]; if (n == $2) c++ } END { print c }' LTRs.alns.results
```
83%

How many does it get within 20bp?
```
awk -F'\t' '{ split($1,a,/[a-zA-Z_~#-]+/); n=a[length(a)]+0; if (n>=($2-20) && n<=($2+20)) c++ } END{print c}' LTRs.alns.results
```
89%

Not bad for 30% divergence (15% * 2).

Performance.
P-distance (raw div).
```
awk '{a[NR]=$6; sum+=$6} END {n=asort(a); median=(n%2 ? a[(n+1)/2] : (a[n/2]+a[n/2+1])/2); print "Mean:", sum/n, "Median:", median}' LTRs.alns.results
```
Mean: 0.24198 Median: 0.244444

JC69.
```
awk '{a[NR]=$8; sum+=$8} END {n=asort(a); median=(n%2 ? a[(n+1)/2] : (a[n/2]+a[n/2+1])/2); print "Mean:", sum/n, "Median:", median}'
LTRs.alns.results
```
Mean: 0.293474 Median: 0.295789

K2P.
```
awk '{a[NR]=$10; sum+=$10} END {n=asort(a); median=(n%2 ? a[(n+1)/2] : (a[n/2]+a[n/2+1])/2); print "Mean:", sum/n, "Median:", median}' LTRs.alns.results
```
Mean: 0.295765 Median: 0.297855


Runtime.
With 100 threads, it processes 20,571 LTR-RTs in 27m:7s.

Idea: Calculate and report p-dist, JC69-dist, and K2P-dist for all LTR-RTs as a single composite measurment using the sum of all lengths, sum of all transversions, etc. this might be better than taking the mean of all LTR-RT divergencews. 
