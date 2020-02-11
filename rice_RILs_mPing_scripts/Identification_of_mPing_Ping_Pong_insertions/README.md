## Identification of mPing, Ping, and Pong insertions in RILs

+ Run RelocaTE2 with mPing, Ping as query

```
python Run_RelocaTE2.py --input fastq --repeat repeat.fa

optional arguments:
  -h, --help            show this help message and exit
  -i INPUT, --input input folder with fastq files
  -o OUTPUT, --output output folder for RelocaTE2 results
  -g GENOME, --genome reference genome (MSU7)
  -r REPEAT, --repeat fasta sequence of mPing, Ping, Pong

```

+ Distinguish Ping from mPing in RelocaTE2 results

```
python CallPing.py --input fastq_RelocaTE2/RIL1_RelocaTE2

optional arguments:
  -h, --help            show this help message and exit
  -i INPUT, --input INPUT
  -o OUTPUT, --output OUTPUT
  --SNP
  -v

```


+ Denovo parental and denovo mPing loci

```
GFF files of mPing loci:
parental (n=466): RILs_ALL.parental_mPing.loci.gff
denovo shared (n=1914): RILs_ALL.denovo_shared_mPing.loci.gff
denovo unique heterzygous (n=4007): RILs_ALL.denovo_unique_mPing_het.loci.gff  
denovo unique homozygous (n=10527): RILs_ALL.denovo_unique_mPing_hom.loci.gff
```
