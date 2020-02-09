GENOME=/PATH_TO_FILE/MSU_r7.fa
HOME=/PATH_TO_FILE/fastq
DESTINATION=/PATH_TO_FILE/genotypes
RENAME_DIR=/PATH_TO_FILE/Illumina
OUTDIR=$(DESTINATION)/MSU_r7
SCRIPTS=/PATH_TO_FILE/fastq
FINAL_RENAME=$(HOME)/final_genotype_renamed.txt
dbSNP=/PATH_TO_FILE/dbSNPs_VQSR.vcf

all: genotype.todo.txt
clean: 
	rm $(FINAL_RENAME)
	rm $(HOME)/genotype.todo.txt 
	
genotype.todo.txt: $(FINAL_RENAME)
	for i in `cat $^` ; do\
	 SHORT=`echo $$i|cut -d ',' -f1`;\
	 FQ_1=`echo $$i|cut -d ',' -f2`;\
	 FQ_2=`echo $$i|cut -d ',' -f3`;\
	 printf "#!/bin/bash\n\
	 date\n\
	 if [ ! -d $(OUTDIR) ]; then\n\
	  mkdir -p $(OUTDIR)\n\
	 fi\n\
	 $(SCRIPTS)/process_reads_genotype.sh $(GENOME) $$FQ_1 $$FQ_2 $$SHORT $(dbSNP) $(OUTDIR)\n\
	 date\n" > $(HOME)/run.genotype.$$SHORT.sh;\
	 JOB=`sbatch -p intel --nodes=1 --ntasks=8 --mem=60G --time=40:00:00 --workdir=./ $(HOME)/run.genotype.$$SHORT.sh`;\
	 echo "`date` $$JOB run.genotype.$$SHORT.sh" >> $@;\
	done;\
        printf "\nThe following jobs have been submitted to the queue:\n";\
	cat $@

$(FINAL_RENAME): $(SCRIPTS)/rename_files.pl $(RENAME_DIR) $(HOME) /PATH_TO_FILE/genotype_sample.txt
	$^ > $@

.PHONY: all clean

