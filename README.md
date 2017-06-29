# SeqSero2 alpha-test version
Salmonella serotyping from genome sequencing data


# Introduction 
SeqSero2 is a pipeline for Salmonella serotype determination from raw sequencing reads or genome assemblies. This is a alpha test version. For now, it can only accept separated paired-end reads. A web app will be available soon.

# Dependencies 
SeqSero depends on:

1. Python 2.7 and [Biopython 1.65](http://biopython.org/wiki/Download); 

2. [Burrows-Wheeler Aligner](http://sourceforge.net/projects/bio-bwa/files/); 

3. [Samtools](http://sourceforge.net/projects/samtools/files/samtools/);

4. [NCBI BLAST](https://blast.ncbi.nlm.nih.gov/Blast.cgi?PAGE_TYPE=BlastDocs&DOC_TYPE=Download);

5. [SRA Toolkit](http://www.ncbi.nlm.nih.gov/Traces/sra/sra.cgi?cmd=show&f=software&m=software&s=software);

6. [SPAdes](http://bioinf.spbau.ru/spades). 

# Executing the code 
    Usage: SeqSero2.py 

    -p <int> (number of threads, if p >4, only 4 threads will be used for assembly since the size of extracted reads is small, default=1)

    -i <file> (/path/to/input/file; for now, SeqSero2 just accepts separated paired-end reads ) 

    -b <string> (algorithms for bwa mapping; 'mem' for mem, 'sam' for samse/sampe; default=mem; optional; for now SeqSero2 is only optimized for "mem" mode) 

# Output 
Upon executing the command, a directory named 'SeqSero_result_<time_you_run_SeqSero>' will be created. Your result will be stored in 'Seqsero_result.txt' in that directory. And the assemblied alleles can also be found in the directory.

# Citation
Zhang S, Yin Y, Jones MB, Zhang Z, Deatherage Kaiser BL, Dinsmore BA, Fitzgerald C, Fields PI, Deng X.  
Salmonella serotype determination utilizing high-throughput genome sequencing data.  
**J Clin Microbiol.** 2015 May;53(5):1685-92.[PMID:25762776](http://jcm.asm.org/content/early/2015/03/05/JCM.00323-15)

