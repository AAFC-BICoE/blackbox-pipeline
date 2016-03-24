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
                sample.general.FastQCoutput = os.path.join(sample.general.outputdirectory, 'fastqc')
                outdir = '{}/fastqc/fastqc{}'.format(os.path.split(fastqfiles[0])[0], level)
                # Separate system calls for paired and unpaired fastq files
                if len(fastqfiles) == 2:
                    # Call fastqc with -q (quiet), -o (output directory), -d (where to store temp files) flags, and
                    # -t (number of threads) flags
                    fastqccall = "fastqc {} {} -q -o {} -t 12 --extract".format(fastqfiles[0], fastqfiles[1], outdir)
                elif len(fastqfiles) == 1:
                    fastqccall = "fastqc {} -q -o {} -t 12 --extract".format(fastqfiles[0], outdir)
                # Record FastQC commands
                sample.commands.FastQC = fastqccall
                # Record FastQC version
                sample.software.FastQC = version
                # Add the arguments to the queue
                fastqfile = os.path.basename(fastqfiles[0]).split('.')[0]
                fastqcfolder = os.path.join(outdir, fastqfile + '_fastqc')
                self.qcqueue.put((fastqccall, fastqcfolder, outdir))
        # Wait on the trimqueue until everything has been processed
        self.qcqueue.join()
        self.qcqueue = Queue()

    def fastqc(self):
        """Run fastqc system calls"""
        from accessoryFunctions import make_path
        while True:  # while daemon
            # Unpack the variables from the queue
            (systemcall, fastqcfolder, outputdir) = self.qcqueue.get()
            # Check to see if the output directory already exists
            if not os.path.isdir(fastqcfolder):
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
                bbdukcall = "bbduk.sh -Xmx1g qtrim=w trimq=25 ktrim=r minlength=50 ftl=10" \
                            " k=25 mink=11 ref={}/resources/adapters.fa hdist=1 ".format(self.bbduklocation)
                # http://seqanswers.com/forums/showthread.php?t=42776
                if len(sample.general.fastqfiles) == 2:
                    fastqfiles.append('{}/{}_R2_trimmed.fastq'.format(outputdir, sample.name))
                    bbdukcall += "in1={} in2={} out1={} out2={} tpe tbo".format(*fastqfiles)
                elif len(sample.general.fastqfiles) == 1:
                    bbdukcall += "in={} out={}".format(*fastqfiles)
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
        import shutil
        printtime('Performing secondary trimming', self.start)
        print "\r[{:}] Trimming fastq files".format(time.strftime("%H:%M:%S"))
        for sample in self.metadata:
            m = int(sample.run.forwardlength) - 50
            if sample.general.trimmedfastqfiles != 'NA':
                bbcall = sample.commands.bbduk
                for fastq in sample.general.trimmedfastqfiles:
                    folder = os.path.join(sample.general.FastQCoutput,
                                          "fastqcTrimmed",
                                          os.path.splitext(os.path.basename(fastq))[0] + '_fastqc')
                    with open(os.path.join(folder, 'fastqc_data.txt')) as summary:
                        newm = maxcut(summary, int(sample.run.forwardlength)) + 10
                    m = newm if newm < m else m
                    if m < (int(sample.run.forwardlength) - 50) and not os.path.isdir("{0:s}_fail".format(folder)):
                        shutil.move(folder, "{0:s}_fail".format(folder))
                        if 'ftr' not in sample.commands.bbduk:
                            bbcall = bbcall.replace('ftl=10', 'ftl=10 ftr=' + str(m))
                        else:
                            import re
                            bbcall = re.sub('ftr=\d+', 'ftr=' + str(m), bbcall)
                if sample.commands.bbduk != bbcall:
                    bbcall = bbcall.replace('hdist=1', 'hdist=2')
                    sample.commands.bbduk = bbcall
                    map(os.unlink, sample.general.trimmedfastqfiles)
                    self.trimqueue.put((sample.commands.bbduk, sample.general.trimmedfastqfiles[0]))
            else:
                self.trimqueue.put(("NA", "NA"))
        # Wait on the trimqueue until everything has been processed
        self.trimqueue.join()

        # Add all the trimmed files to the metadata
        self.fastqcthreader('Trimmed')



    def bbduker(self):
        """Run bbduk system calls"""
        while True:  # while daemon
            # Unpack the variables from the queue
            (systemcall, forwardname) = self.trimqueue.get()
            # Check to see if the forward file already exists
            if not os.path.isfile(forwardname) and systemcall != "NA" and systemcall:
                execute(systemcall, shell=True)
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


def maxcut(summary, full):
    # Set flag to find file location from fastqc
    perbase = False
    variance = [[0.0, 0.0], [0.0, 0.0], [0.0, 0.0], [0.0, 0.0]]
    m = int(full)
    # Parse the summary file and look for the correct line
    for line in summary:
        if "Per base sequence content" in line:
            # Find the line, set the flag then skip past the header
            perbase = True
            summary.next()
        elif perbase:
            # Reset the flag once we find the last line of the important data
            if "END_MODULE" in line:
                perbase = False
            else:
                # extract teh variables, python3 supports partial unpack so we could integrate this in the future
                location, g, a, t, c = line.rstrip().split()
                if '-' in location:
                    # most of the lines are a range from m to n seperated by a dash
                    m, n = map(int, location.split('-'))
                else:
                    # First lines are a single digit
                    n = int(location)
                    m = n
                if n >= 54:
                    # Only start after the 54 to calculate the stdev thereby weighting the errors on the end
                    values = map(float, [g, a, t, c])
                    for v in range(4):
                        # online variance calculation
                        mean, m2 = variance[v]
                        delta = values[v] - mean
                        mean += delta / (n - 53)
                        m2 += delta * (values[v] - mean)
                        s = m2 / (n - 53)
                        if s > 0.0035 and mean and n > 50:
                            # we found that the above value acts for a 1% variation
                            return m
                        else:
                            # add the variance to the appropraite list
                            variance[v] = [mean, m2]
    return m
