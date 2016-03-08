#!/usr/bin/env python
import os
from threading import Thread
from Queue import Queue
import time
from accessoryFunctions import *

__author__ = 'adamkoziol,mikeknowles'


class Quality(object):
    def fastqcthreader(self, level):
        printtime('Running quality control on {} fastq files'.format(level), self.start)
        version = get_version(['fastqc', '-v']).rstrip()
        for sample in self.metadata:
            if type(sample.general.fastqfiles) is list:
                # Create and start threads for each fasta file in the list
                # Send the threads to bbduker. :args is empty as I'm using
                threads = Thread(target=self.fastqc, args=())
                # Set the daemon to true - something to do with thread management
                threads.setDaemon(True)
                # Start the threading
                threads.start()
        # Iterate through strains with fastq files to set variables to add to the multithreading queue
        for sample in self.metadata:
            fastqccall = ""
            # Check to see if the fastq files exist
            if level == 'Trimmed':
                # Try except loop to allow for missing samples
                try:
                    fastqfiles = sample.general.trimmedfastqfiles
                except KeyError:
                    fastqfiles = ""
                    pass
            else:
                fastqfiles = sample.general.fastqfiles
            # As the metadata can be populated with 'NA' (string) if there are no fastq files, only process if
            # :fastqfiles is a list
            if type(fastqfiles) is list:
                # Set the output directory location
                outdir = '{}/fastqc/fastqc{}'.format(os.path.split(fastqfiles[0])[0], level)
                # Separate system calls for paired and unpaired fastq files
                if len(fastqfiles) == 2:
                    # Call fastqc with -q (quiet), -o (output directory), -d (where to store temp files) flags, and
                    # -t (number of threads) flags
                    fastqccall = "fastqc {} {} -q -o {} -t 12".format(fastqfiles[0], fastqfiles[1], outdir)
                elif len(fastqfiles) == 1:
                    fastqccall = "fastqc {} -q -o {} -t 12".format(fastqfiles[0], outdir)
                # Record FastQC commands
                sample.commands.FastQC = fastqccall
                # Record FastQC version
                sample.software.FastQC = version
                # Add the arguments to the queue
                self.qcqueue.put((fastqccall, outdir))
        # Wait on the trimqueue until everything has been processed
        self.qcqueue.join()
        self.qcqueue = Queue()

    def fastqc(self):
        """Run fastqc system calls"""
        from accessoryFunctions import make_path
        while True:  # while daemon
            # Unpack the variables from the queue
            (systemcall, outputdir) = self.qcqueue.get()
            # Check to see if the output directory already exists
            if not os.path.isdir(outputdir):
                # Make the output directory
                make_path(outputdir)
                # Run the system call
                execute(systemcall)
                # call(systemcall, shell=True, stdout=devnull, stderr=devnull)
            # Signal to qcqueue that job is done
            self.qcqueue.task_done()

    def trimquality(self):
        """Uses bbduk from the bbmap tool suite to quality and adapter trim"""
        version = get_version(['bbduk.sh', '-version']).split('\n')[-3].split()[-1]
        from glob import glob
        print "\r[{:}] Trimming fastq files".format(time.strftime("%H:%M:%S"))
        # Create and start threads for each strain with fastq files
        for sample in self.metadata:
            if type(sample.general.fastqfiles) is list:
                # Create and start threads for each fasta file in the list
                # Send the threads to bbduker. :args is empty as I'm using
                threads = Thread(target=self.bbduker, args=())
                # Set the daemon to true - something to do with thread management
                threads.setDaemon(True)
                # Start the threading
                threads.start()
        # Iterate through strains with fastq files to set variables to add to the multithreading queue
        for sample in self.metadata:
            # As the metadata can be populated with 'NA' (string) if there are no fastq files, only process if
            # :fastqfiles is a list
            if type(sample.general.fastqfiles) is list:
                # Check to see if the fastq files exist
                fastqfiles = sorted(sample.general.fastqfiles)
                # Define the output directory
                outputdir = sample.general.outputdirectory
                # Define the name of the forward trimmed fastq file
                fastqfiles.append('{}/{}_R1_trimmed.fastq'.format(outputdir, sample.name))
                # Separate system calls for paired and unpaired fastq files
                # TODO minlen=number - incorporate read length
                bbdukcall = "bbduk.sh -Xmx1g qtrim=w trimq=20 ktrim=l minlength=50" \
                            " k=25 mink=11 ref={}/resources/adapters.fa hdist=1".format(self.bbduklocation)
                # http://seqanswers.com/forums/showthread.php?t=42776
                if len(fastqfiles) == 2:
                    fastqfiles.append('{}/{}_R2_trimmed.fastq'.format(outputdir, sample.name))
                    bbdukcall += "in1={} in2={} out1={} out2={} tpe tbo".format(*fastqfiles)
                elif len(fastqfiles) == 1:
                    bbdukcall += "bbduk.sh in={} out={}".format(*fastqfiles)
                else:
                    bbdukcall = ""
                # Record bbMap commands
                sample.commands.bbduk = bbdukcall if bbdukcall else "NA"
                # Record FastQC version
                sample.software.bbduk = version
                # Add the arguments to the queue
                self.trimqueue.put((bbdukcall, fastqfiles[2]))
        # Wait on the trimqueue until everything has been processed
        self.trimqueue.join()
        # Add all the trimmed files to the metadata
        for sample in self.metadata:
            # Define the output directory
            outputdir = sample.general.outputdirectory
            # Add the trimmed fastq files to a list
            trimmedfastqfiles = glob('{}/*trimmed.fastq'.format(outputdir, sample.name))
            # Populate the metadata if the files exist
            sample.general.trimmedfastqfiles = trimmedfastqfiles if trimmedfastqfiles else 'NA'
        print "\r[{:}] Fastq files trimmed".format(time.strftime("%H:%M:%S"))
        self.fastqcthreader('Trimmed')

    def bbduker(self):
        """Run bbduk system calls"""
        while True:  # while daemon
            # Unpack the variables from the queue
            (systemcall, forwardname) = self.trimqueue.get()
            # Check to see if the forward file already exists
            if not os.path.isfile(forwardname):
                execute(systemcall)
            # Signal to trimqueue that job is done
            self.trimqueue.task_done()

    def __init__(self, inputobject):
        self.metadata = inputobject.runmetadata.samples
        self.qcqueue = Queue()
        self.trimqueue = Queue()
        self.correctqueue = Queue()
        self.start = inputobject.starttime
        # Find the location of the bbduk.sh script. This will be used in finding the adapter file
        self.bbduklocation = os.path.dirname(which('bbduk.sh'))
