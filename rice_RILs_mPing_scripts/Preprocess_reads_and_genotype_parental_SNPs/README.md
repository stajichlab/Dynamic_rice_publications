## Preprocess Illumina paired-end reads and genotyping parental SNPs 

+ Generate a Makefile to run preprocessing Illumina reads and genotying parental SNPs

```
perl make_genotype_makefile.pl -r /PATH_TO_FILE/Illumina -o /PATH_TO_FILE/genotypes -f . -i genotype_sample.txt -q batch -l "nodes=1:ppn=24,mem=100g,walltime=50:00:00"

Assistants in creating the config file for makefile to process reads for genotyping.

usage: $0 -o output_dir -f fastq_dir -i barcode_sample_id_file -q qsub_options -h this message

  o|output_dir:s        => name of directory for all outputed directories and files
  f|fastq_dir:s         => name of directory of fastq files to be processed
  r|rename_fastq_dir:s  => name of directory of in which the renamed fastq files will be placed
  i|barcode_sample:s    => name of file with Sample_names and illumina barcode ids
  g|genome:s            => name of genome file 
  d|dbSNPs:s            => name of dbSNP file 
  q|qsub_q:s            => qsub -q options, if a queue is available. ex: -q highmem
  l|qsub_l:s            => qsub -l options, if a queue is available. ex: -l nodes=1:ppn=8
  h|help                => this message

```

+ Excute the file "Makefile" to run the preproess and genotyping.

```
make all
```
