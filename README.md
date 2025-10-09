#### Pipeline to estimate LTR-RT insertion time.   
(1) Use kmers to identify LTRs from a library of LTR-RTs.  
(2) MAFFT and trimal to clean alignments.  
(2) Force global alignment gapless LTRs.   
(3) Date with p-dist, JC69, and K2P.  

```diff
Process multi-seq LTR-RT FASTA(s) to extract and align LTRs. 

options:
  -h, --help            show this help message and exit.
  -k                Keep temp directory after processing.
  -v                    Verbose mode; print each command before executing.
  -d DIST               Minimum distance between kmer pairs (default: 80).
  -l KMIN               Minimum kmer length (default: 8).
  -U KMAX               Maximum kmer length (default: 12).
+ -u MUTATION_RATE      Mutation rate μ (default: 3e-8).
  -f STD_FACTOR         Standard deviation factor for kmer filtering.
  -e EXTENSION          Extension length for LTR extraction (default: 65).
  -t TEMP_DIR           Temporary directory name (single input only; ignored with multiple inputs).
+ -o OUTFILE            Output filename (single input only; ignored with multiple inputs).
+ -p THREADS            Number of parallel threads per input (default: 20).
+ -D [DOMAINS_TSVS ...], --domains [DOMAINS_TSVS ...]
                        Optional domains TSV file(s) (format: name\tLTR_len). With a single input, provide one TSV. With multiple inputs, you may
                        supply multiple TSVs; each TSV is matched to an input by prefix. Matching uses the TSV filename up to the first '.'; if
                        it ends with '_domains', that suffix is ignored for matching.
  --max-win-overdisp MAX_WIN_OVERDISP
                        If set, exclude alignments with window overdispersion (win_overdisp) greater than this value.
  --min-retained-fraction MIN_RETAINED_FRACTION
                        Minimum fraction of ungapped columns retained after trimming required to proceed (default: 0.6).
  --assume-duplicate-same-ltr
                        Override duplicate-header safety fallback in fast path. Assumes all duplicate header tokens share the same LTR length
                        from the domains TSV. Use with caution.
  --no-plot             Disable plotting of results into kmer2ltr_density.pdf.
+ -i INPUT_FASTAS [INPUT_FASTAS ...], --input-fastas INPUT_FASTAS [INPUT_FASTAS ...]
                        Path(s) to multi-sequence LTR-RT FASTA file(s).
```

Most users probably only need the flags that are highlighted in green.

Single-file input:
```
python Kmer2LTR/Kmer2LTR.py -i test.fa -u 1.8e-8
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

`LTRs.alns.results.summary` shows the cummulative numbers.
```
total_length	22408270
total_transitions	532951
total_transversions	217017
raw_d	0.033468
JC69_d	0.034238
K2P_d	0.034368
```


Multi-file input with 50 threads:
```
python Kmer2LTR/Kmer2LTR.py -i species1.ltr.fa species2.ltr.fa species3.ltr.fa -p 50
```

We could save time if we already know the lengths of the LTR:
```
ls *.domains
species1.ltr.domains  species2.ltr.domains  species3.ltr.domains


head -1 species1.ltr.domains
CMHA_chr1:90368..96317	172


# Tells the code that theres an LTR-RT in 'species1.ltr.fa' named 'CMHA_chr1:90368..96317'.
# Its LTRs are each 172bp in length.

# With domain input, we can save time by skipping LTR identification.
python Kmer2LTR/Kmer2LTR.py -i species1.ltr.fa species2.ltr.fa species3.ltr.fa -p 50 -D *.domains
```

- `--max-win-overdisp` and `--min-retained-fraction` are filtration flags to remove candidate LTR-RTs if they appear dubious.    
  - `--max-win-overdisp` is used to exclude LTR-RTs if their LTRs are inconsistent diveregence. Eg, the 5' half of the LTRs have low divergence, but the 3' half appears quite divergence. This can happen if the LTR boundaries are not accurate. Start with `--max-win-overdisp 6` and lower for additional stringency.  
  - `--min-retained-fraction` is used to exclude LTR-RTs if the mafft/trimal-cleaned LTRs are shorter than threshold of the kmer-defined LTR boundary. Start with `--min-retained-fraction 0.5` or `--min-retained-fraction 0.6`.  

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

Convert pass list to domains file. 
```
for f in *.pass.list; do out="${f%.pass.list}.domains"; python pass_list_domians.py "$f" > "$out"; done
```

Convert pass list and genome to LTR-RT file:
```
for tsv in *.pass.list; do
    prefix="${tsv%.pass.list}"
    fasta="${prefix}.fa" 
    out="${tsv}.fa" 
    python pass_list_fa_extractor.py -fa "$fasta" -tsv "$tsv" > "$out" &
done
wait
```
