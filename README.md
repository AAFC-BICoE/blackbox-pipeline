[Parallel ITSx]:https://github.com/AAFC-MBB/parallel_itsx 


Blackbox Pipeline for genome assembly
==============

## Background
### Introduction

This is an automated genome assembly pipeline designed to optimally run inside a Docker image.  The purpose of this pipeline is to provide researchers with a standard and simplified workflow for microbial genome assembly while keeping track of executed commands and their parameters, program version numbers and associated sample metadata, which is important for ensuring experimental reproducibility and traceability.

This pipeline uses SPAdes as the genome assembler and therefore would work best with bacterial or smaller fungal genomes (< 100 Mb).  MiSeq raw reads and metadata files are required as input.  Alternatively, archived files from BaseSpace (e.g. analysis_14348334_fastq.zip) can also be used as input.  In addition to the forward and reverse reads file, the following files are required:

1. GenerateFASTQRunStatistics.xml
2. RunInfo.xml
3. SampleSheet.csv

These files are located within the appropriate subfolder (e.g. 140922_M02466_0030_000000000-AARWU - the naming of this folder consists of the date (140922), the MiSeq designation (M02466), and the flowcell number (000000000-AARWU)) of the MiSeqOutput
directory in the MiSeq onboard computer. The fastq files are located in the ../MiSeqOutput/140922_M02466_0030_000000000-AARWU/Data/Intensities/BaseCalls
folder.

Copy the forward and reverse fastq files, GenerateFASTQRunStatistics.xml, RunInfo.xml and SampleSheet.csv to a different working location (e.g. ../Sequencing/user_name/project_name/genome_assembly/140922).

### Contents 

[MBBspades]: bin/MBBspades
[SPAdes]: http://spades.bioinf.spbau.ru 

This pipeline includes a main script ([MBBSpades]) that executes the following helper modules:

* [Metadata Reader](blackbox/runMetadata.py): Extracts metadata from sequencing run reports/files 
* [Fastq Mover](blackbox/fastqmover.py): Moves and/or extracts archived files 
* [SPAdes Run](blackbox/spadesRun.py): Runs SPAdes assembler with error correction enabled by default as well as estimates the insert size of input libraries.
* [QUAST Parser](blackbox/quastParser.py): Determines assembly quality metrics using QUAST 
* [Bowtie2 Wrapper](blackbox/bowtie.py): Subclass of [Biopython's 
  AbstractCommandline](http://biopython.org/DIST/docs/api/Bio.Application.AbstractCommandline-class.html) for 
  [Bowtie2](http://bowtie-bio.sourceforge.net/bowtie2/index.shtml) with support for piping to SAMtools for sorted bam output 
* [Metadata Writer](blackbox/metadataprinter.py): Creates a JSON report of all metadata for each sequenced strain 
* [QualiMap Parser](blackbox/qualimapR.py): Wrapper and parser for [qualimap](http://qualimap.bioinfo.cipf.es/) using the bowtie2 wrapper that maps the corrected reads from SPAdes 
* [ITSx Parser](blackbox/its.py): Utilizes [Parallel ITSx] and parses the results. This will not run if `--clade` is set to `bacteria`

## Installation

### Requirements 

#### Docker 

* This is the ideal solution.  See the [parent project](https://github.com/AAFC-MBB/docker-assembly) for details on how to install and set up.

## Outputs

1. Assembled contigs are collected in the 'BestAssemblies' folder
2. Reports in JSON format are located in the genome folder with the suffix _\_metadata.json_

## Usage

```
usage: MBBSpades [-h] [-v] [-n numreads] [-t threads] [-o] [-F]
                                   [-d destinationfastq] [-m miSeqPath] [-f miseqfolder]
                                   [-r1 readLengthForward] [-r2 readLengthReverse]
                                   [-r referenceFilePath] [-k kmerRange]
                                   [-c customSampleSheet] [-b] [--clade CLADE]
                                   [--itsx ITSX] [--trimoff]
                                   path

Assemble genomes from Illumina fastq files

positional arguments:
  path                  Specify path

optional arguments:
  -h, --help            show this help message and exit
  -v, --version         show program's version number and exit
  -n numreads           Specify the number of reads. Paired-reads: 2,
                        unpaired-reads: 1. Default is paired-end
  -t threads            Number of threads. Default is the number of cores in
                        the system
  -o, --offHours        Optionally run the off-hours module that will search
                        for MiSeq runs in progress, wait until the run is
                        complete, and assemble the run
  -F, --FastqCreation   Optionally run the fastq creation modulethat will
                        search for MiSeq runs in progress, run bcl2fastq to
                        create fastq files, and assemble the run
  -d destinationfastq   Optional folder path to store .fastq files created
                        using the fastqCreation module. Defaults to
                        path/miseqfolder
  -m miSeqPath          Path of the folder containing MiSeq run data folder
  -f miseqfolder        Name of the folder containing MiSeq run data
  -r1 readLengthForward
                        Length of forward reads to use. Can specify "full" to
                        take the full length of forward reads specified on the
                        SampleSheet. Defaults to full
  -r2 readLengthReverse
                        Length of reverse reads to use. Can specify "full" to
                        take the full length of reverse reads specified on the
                        SampleSheet. Defaults to full
  -r referenceFilePath  Provide the location of the folder containing the
                        pipeline accessory files (reference genomes, MLST
                        data, etc.
  -k kmerRange          The range of kmers used in SPAdes assembly. Default is
                        21,33,55,77,99,127
  -c customSampleSheet  Path of folder containing a custom sample sheet and
                        name of sample sheet file e.g.
                        /home/name/folder/BackupSampleSheet.csv. Note that
                        this sheet must still have the same format of Illumina
                        SampleSheet.csv files
  -b, --basicAssembly   Performs a basic de novo assembly, and does not
                        collect metadata
  --clade CLADE         Specifiy HMM database for BUSCO
  --itsx ITSX           Specifiy comma-seperated HMM database for ITSx
  --trimoff             Turn off trimming with bbduk
```
