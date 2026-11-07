# SeqSero2
Salmonella serotype prediction from genome sequencing data.

Online version: http://www.denglab.info/SeqSero2

# Introduction 
SeqSero2 is a pipeline for Salmonella serotype prediction from raw sequencing reads or genome assemblies

# Dependencies and installation

Python >=3.8 and Biopython >=1.75 are required. Biopython 1.74 has a historical BLAST XML parser problem; it is excluded by the minimum version. The validation report records tested versions rather than implying that every version above a minimum was tested.

| Workflow | Additional programs |
|---|---|
| Raw-read allele microassembly (`-m a`, default) | BWA >=0.7.12, samtools >=1.8, BLAST+ >=2.2.28 (`blastn`, `makeblastdb`), SPAdes >=3.15.5, SalmID 0.1.23 (`SalmID.py`) |
| Raw-read k-mer (`-m k -t 1/2/3`) | No external tool for FASTQ; SRA Toolkit >=2.8.0 only when supplying `.sra` archives |
| Assembly k-mer (`-m k -t 4`) | No external tool for FASTA |
| Nanopore k-mer (`-m k -t 5`) | seqtk >=1.3 when converting FASTQ to FASTA |

Bedtools is no longer called: current read extraction uses `samtools bam2fq`. `--check` checks the union of supported workflow executables and returns a nonzero status when any are missing; a missing optional program does not prevent workflows that do not use it.

Create an environment and install this checkout:

```sh
conda env create -f environment.yml
conda activate seqsero2
python -m pip install .
SeqSero2_package.py --check
```

On Apple Silicon, create a native environment with `conda env create --platform osx-arm64 -f environment.yml`. Existing published Bioconda builds may retain their old dependency pins; installing this checkout does not update those external builds. `packaging/conda-meta.proposed.yaml` is a proposed local-source recipe, not a published package.

# Executing the code 
Make sure all SeqSero2 and its dependency executables are added to your path (e.g. to ~/.bashrc). Then type SeqSero2_package.py to get detailed instructions.

    Usage: SeqSero2_package.py 

    -m <string> (which workflow to apply, 'a'(raw reads allele micro-assembly), 'k'(raw reads and genome assembly k-mer), default=a)

    -t <string> (input data type, '1' for interleaved paired-end reads, '2' for separated paired-end reads, '3' for single reads, '4' for genome assembly, '5' for nanopore reads (fasta/fastq))

    -i <file> (/path/to/input/file)

    -p <int> (number of threads for allele mode, if p >4, only 4 threads will be used for assembly since the amount of extracted reads is small, default=1) 

    -b <string> (algorithms for bwa mapping for allele mode; 'mem' for mem, 'sam' for samse/sampe; default=mem; optional; for now we only optimized for default "mem" mode)
 
    -d <string> (output directory name, if not set, the output directory would be 'SeqSero_result_'+time stamp+one random number)
	
    -c <flag> (if '-c' was flagged, SeqSero2 will only output serotype prediction without the directory containing log files)
    
    -n <string> (optional, to specify a sample name in the report output)
    
    -s <flag> (if '-s' was flagged, SeqSero2 will not output header in SeqSero_result.tsv)
		    
    --check <flag> (use '--check' flag to check the required dependencies)
    
    -v, --version (show program's version number and exit)
	

# Examples
Allele mode:

    # Allele workflow ("-m a", default), for separated paired-end raw reads ("-t 2"), use 10 threads in mapping and assembly ("-p 10")
    SeqSero2_package.py -p 10 -t 2 -i R1.fastq.gz R2.fastq.gz
	
K-mer mode:

    # Raw reads k-mer ("-m k"), for separated paired-end raw reads ("-t 2")
    SeqSero2_package.py -m k -t 2 -i R1.fastq.gz R2.fastq.gz

    # Genome assembly k-mer ("-t 4", genome assemblies only predicted by the k-mer workflow, "-m k")
    SeqSero2_package.py -m k -t 4 -i assembly.fasta
	
# Output 
Upon executing the command, a directory named 'SeqSero_result_Time_your_run' will be created. Your result will be stored in 'SeqSero_result.txt' in that directory. And the assembled alleles can also be found in the directory if using "-m a" (allele mode).


# Citation
Zhang S, Den-Bakker HC, Li S, Dinsmore BA, Lane C, Lauer AC, Fields PI, Deng X. 
SeqSero2: rapid and improved Salmonella serotype determination using whole genome sequencing data.
**Appl Environ Microbiology. 2019 Sep; 85(23):e01746-19.** [PMID: 31540993](https://aem.asm.org/content/early/2019/09/17/AEM.01746-19.long) 

Zhang S, Yin Y, Jones MB, Zhang Z, Deatherage Kaiser BL, Dinsmore BA, Fitzgerald C, Fields PI, Deng X.  
Salmonella serotype determination utilizing high-throughput genome sequencing data.  
**J Clin Microbiol. 2015 May;53(5):1685-92.** [PMID: 25762776](http://jcm.asm.org/content/early/2015/03/05/JCM.00323-15)

# Workflow and input notes

The allele workflow maps reads to the determinant reference, extracts mapped reads, assembles them with SPAdes, scores BLAST alignments, assigns O/H determinants and then applies the existing serotype interpretation table. Species/subspecies identification in this workflow is delegated to SalmID, using the forward input. A tool failure returns a nonzero exit and is not converted into a successful antigen-negative report.

The raw-read k-mer workflow uses selected 27-mers from the forward input for speed. The paired-end report lists both submitted paths as provenance; this does not mean both mates were sampled. The assembly workflow compares reference 27-mers with assembly sequence directly. It does not remap reads, assemble, or run BLAST. These are different scoring methods and can disagree on incomplete or divergent determinants.

K-mer subspecies identification uses the bundled invA marker sets; it is not an abundance estimate or a validated mixed-subspecies classifier. Retain taxonomic/QC evidence when results disagree. Updating subspecies taxonomy or training references requires a separate reference validation.

Use `-t 2 -i forward.fastq.gz reverse.fastq.gz` for separated mates, `-t 1` for a single interleaved FASTQ, `-t 3` for single reads, `-m k -t 4 -i assembly.fna` for assemblies, and `-m k -t 5` for nanopore FASTA/FASTQ. FASTA starts with `>name`; FASTQ has a `@name`, sequence, `+`, and matching-length quality line. Supply readable file paths. For Docker, mount the host directory and use paths inside the container, for example `-v /absolute/reads:/data`, then `-t 2 -i /data/R1.fastq.gz /data/R2.fastq.gz`.

Output directories must be new, including with `-c`; existing files are never overwritten or removed. Quote paths containing whitespace or shell metacharacters. The program passes literal arguments to tools and uses simple temporary paths for SPAdes/BLAST.

Paratyphi B marker scoring uses the strongest positive and negative evidence across the original and added context probes. It is independent of MLST. Missing or tied state evidence is reported as an unresolved tartrate variant. The additional negative context probe is a designed SNP counterpart, not a naturally sampled ST28 negative isolate. See the issue resolution report for tested samples and limits.

# Combining results

```sh
SeqSero2_summarize.py results/*/SeqSero_result.tsv -o combined.tsv
```

The helper joins columns by header name so allele and k-mer reports can be combined. Use reports generated without `-s` (headers are required). Existing summary files are not overwritten. For many genomes, loop over one assembly per invocation and include `-m k -t 4`; shell globs that exceed the OS argument limit must be processed in a loop.
