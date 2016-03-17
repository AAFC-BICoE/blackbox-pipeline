[Parallel ITSx]:https://github.com/AAFC-MBB/parallel_itsx 


SPAdesPipeline
==============
##### Table of Contents

<ul>
<li><a href="#background">Background</a>
<ul>
<li><a href="#introduction">Introduction</a></li>
<li><a href="#contents">Contents</a></li>
</ul></li>
<li><a href="#installation">Installation</a><ul>
<li><a href="#untestednon-dockerinstallation">Untested non-Docker installation</a></li>
<li><a href="#requirements">Requirements</a><ul>
<li><a href="#docker">Docker</a></li>
<li><a href="#non-dockerinstallation">non-Docker installation</a></li>
</ul></li>
</ul></li>
<li><a href="#outputs">Outputs</a></li>
<li><a href="#usage">Usage</a></li></ul>

## Background
### Introduction

This pipeline is designed to be able to assemble and type raw fastq data generated using an Illumina MiSeq.

It has been designed to run using either the archived files obtained from BaseSpace (something similar to analysis_14348334_fastq.zip),
alternatively, it is possible to use the individual fastq files taken directly from the MiSeq. This is probably better, as,
in addition to the fastq files, three other files from the MiSeq are required for the pipeline to function:

1. GenerateFASTQRunStatistics.xml
2. RunInfo.xml
3. SampleSheet.csv

These files are located within the appropriate subfolder (e.g. 140922_M02466_0030_000000000-AARWU - this folder consists
of the date (140922), the MiSeq designation (M02466), and the flowcell number (000000000-AARWU)) of the MiSeqOutput
directory in the MiSeq onboard computer. The fastq files are located in the ../MiSeqOutput/140922_M02466_0030_000000000-AARWU/Data/Intensities/BaseCalls
folder.

Copy all necessary files to a properly named folder in an easy to remember location (e.g. ../Sequencing/140922).

### Contents 

[MBBspades]: bin/MBBspades
[SPAdes]: http://spades.bioinf.spbau.ru 

This pipeline includes a main script ([MBBSpades]), and which relies on helper modules located in helper scripts, 
including: 

* [Metadata Reader](blackbox/runMetadata.py): Pulling metadata from sequencing run reports/files 
* [Fastq Mover](blackbox/fastqmover.py): Moving and/or extracting archived files 
* [SPAdes Run](blackbox/spadesRun.py) Running SPAdes assembler with error correction enabled by default 
    * Estimating the size of the library fragments
* [Quast Parser](blackbox/quastParser.py): Determining assembly quality metrics using quast 
* [Bowtie2 Wrapper](blackbox/bowtie.py): Subclass of [Biopython's 
  AbstractCommandline](http://biopython.org/DIST/docs/api/Bio.Application.AbstractCommandline-class.html) for 
  [Bowtie2](http://bowtie-bio.sourceforge.net/bowtie2/index.shtml) with support for piping to SAMtools for sorted bam output 
* [Busco Parser](blackbox/BuscoParser.py): Wrapper and parser of [BUSCO](http://busco.ezlab.org) 
* [Metadata Writer](blackbox/metadataprinter.py): Creating a JSON report of all collected metadata for each sequenced strain 
* [QualiMap Parser](blackbox/qualimapR.py): Wrapper and parser for [qualimap](http://qualimap.bioinfo.cipf.es/) with use 
  of bowtie2 wrapper using corrected reads from SPAdes 
* [ITSx Parser](blackbox/its.py): Utilizes [Parallel ITSx] and parses result 

## Installation
### Untested non-Docker installation
After cloning the git: 

```commandline
git clone https://github.com/AAFC-MBB/blackbox-pipeline
```

Install the python package:

```commandline
cd blackbox-pipeline
python setup.py install
```

[MBBSpades] will now be in your `$PATH`

### Requirements 

#### Docker 

* This is the ideal solution see the [parent project](https://github.com/AAFC-MBB/docker-assembly) 
* Everything else should be contained within the docker container, and is ready to run. 

#### non-Docker installation

* [Parallel ITSx]
* Quast and SPAdes apart of the `$PYTHONPATH`
* All tools contained in [accessoryfiles 
  script](https://github.com/AFFC-MBB/docker-assembly/blob/master/accessoryfiles.sh) added to the `$PATH` environment 
* Requirements in [Dockerfile](https://github.com/AFFC-MBB/docker-assembly/blob/master/Dockerfile) in `apt-get`

## Outputs
This pipeline generates multiple outputs.

1. Assembled contigs - these are collected in the 'BestAssemblies' folder
2. JSON reports - these are located in the genome folder with the suffix _\_metadata.json_ 

Additionally, within the individual strain subfolders, a .pdf output of plotted insert sizes is included in the 'insertSizes' folder.
Detailed reports can be found in the 'quast_results' folder, and the reference genome file is located in 'referenceGenome'

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