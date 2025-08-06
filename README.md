TESS-PrinTE makes `lib_clean.fa` with LTR length in the header.   
All LTRs are unmutated.   

PIPELINE OUTLINE: (1) Use kmers to identify LTR boundry from LTR-RTs. (2) Force global alignment of LTRs (gaps penalized to promote subsitution).  

For benchmarking the pipeline, I can:
```
python lib_mutator.py -i lib_clean.fa -o lib_clean.15mp.titv2.fa -mp 15 -TiTv 2 --seed 51
```
LTRs are 15% mutated. TiTv ratio = 2. 
Internal sequence is not mutated. 

This new seq_divergence pipeline has issues with ':' and '/' in header. 
I will fix if the pipeline is adopted. 
But for now...
```
cat lib_clean.15mp.titv2.fa | sed 's/\//_/g' | sed 's/LTRlen:/LTRlen/g' > lib_clean.15mp.titv2.clean.fa
```

Run the pipeline.
```
nohup /usr/bin/time -v python seq_divergence.py -p 100 lib_clean.15mp.titv2.clean.fa &
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
