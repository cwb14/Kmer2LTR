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
`LTR_LEN` is whats discovered via kmer (or provided in the domains file).  
`ALN_LEN` is the length of aligned bases.  
`LTR_LEN` may be less than or greater than `ALN_LEN`.  
`LTR_LEN` > `ALN_LEN` if gaps in LTR.   
`LTR_LEN` < `ALN_LEN` if extension expands the LTR boundry discovered by kmer.   


`OUTFILE` can be directly used as `DOMAINS_TSV`.













# Developers note.
TESS-PrinTE makes `lib_clean.fa` with LTR length in the header. All LTRs are unmutated. 
Other libraries in TESS are also helpful: `maizeTE02052020.ltr.full`, `rice7.0.0.liban.ltr.full`, `ltr-db.fa`, `athrep.updated.fasta_062024.ltr.full`.

For benchmarking the pipeline, I can:
```
python lib_mutator.py -i lib_clean.fa -o lib_clean.15mp.titv2.fa -mp 15 -TiTv 2 --seed 51
```
LTRs are 15% mutated. TiTv ratio = 2. 
Internal sequence is not mutated. 

Runtime.
With 100 threads, it processes 20,703 LTR-RTs in 25m:22s.
With DOMAINS_TSV provided, and 100 threads, it processes those 20,703 LTR-RTs in 3m:28s.
