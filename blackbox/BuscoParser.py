#!/usr/bin/env python
from accessoryFunctions import *
import os
import quast

__author__ = 'mikeknowles,akoziol'


class Busco(object):

    def buscoprocess(self):
        from threading import Thread
        os.chdir(self.path)
        # Find the fasta files for each sample
        # Only make as many threads are there are samples with fasta files
        for i in range(len([sample.general for sample in self.metadata if sample.general.bestassemblyfile])):
            # Send the threads to the merge method. :args is empty as I'm using
            threads = Thread(target=self.analyze, args=())
            # Set the daemon to true - something to do with thread management
            threads.setDaemon(True)
            # Start the threading
            threads.start()
        for sample in self.metadata:
            # Save augustus, blast and BUSCO versions
            sample.software.BUSCO, sample.software.Blastn, sample.software.Augustus =\
                self.version, self.blast, self.augustus
            if sample.general.bestassemblyfile:
                sample.general.buscoresults = '{}/busco_results'.format(sample.general.outputdirectory)
                sample.commands.BUSCO = "python3 {} -in {} -o {} -l {} -m genome".\
                    format(self.executable, sample.general.bestassemblyfile, sample.name, self.lineage)
                if os.path.isdir("{0:s}/referencegenome".format(self.path)):
                    from glob import glob
                    referencegenome = glob("{0:s}/referencegenome/*".format(self.path))
                    sample.commands.Quast = "quast.py -R {0:s} --gage {1:s} -o {2:s}".\
                        format(referencegenome[0], sample.general.bestassemblyfile, sample.general.quastresults)
                else:
                    sample.commands.Quast = "quast.py {0:s} -o {1:s}".\
                        format(sample.general.bestassemblyfile, sample.general.quastresults)
                self.qqueue.put(sample)
            else:
                sample.commands.QuastCommand = "NA"
        self.qqueue.join()

    def analyze(self):
        """Run the quast command in a multi-threaded fashion"""
        while True:
            sample = self.qqueue.get()
            if sample.general.bestassemblyfile != 'NA' \
                    and not os.path.isfile('{}/report.tsv'.format(sample.general.quastresults)):
                log = os.path.join(sample.general.quastresults, 'stdout.log')
                # sys.stdout, sys.stderr = out, err
                execute(sample.commands.Quast)
            if os.path.isfile('{}/report.tsv'.format(sample.general.quastresults)):
                self.metaparse(sample)
            # Signal to the queue that the job is done
            self.qqueue.task_done()

    def metaparse(self, sample):
        repls = ('>=', 'Over'), ('000 Bp', 'kbp'), ('#', 'Num'), \
                ("'", ''), ('(', ''), (')', ''), (' ', ''), ('>', 'Less'), ('Gc%', 'GC%')
        if not os.path.isfile("%s/report.tsv" % sample.general.quastresults):
            print "There was an issue getting the metadata from {0:s}".format(sample.name)
        else:
            quast = dict()
            resfile = "{0:s}/gage_report.tsv".format(sample.general.quastresults) \
                if os.path.isfile("{0:s}/gage_report.tsv".format(sample.general.quastresults)) \
                else "{0:s}/report.tsv".format(sample.general.quastresults)
            with open(resfile) as report:
                report.next()
                for line in report:
                    # Use headings in report as keys for the GenObject supplied from generator and replace incrementally
                    # with reduce and lambda function below
                    k, v = [reduce(lambda a, kv: a.title().replace(*kv), repls, s) for s in line.rstrip().split('\t')]
                    quast[k] = v
            sample.assembly = GenObject(quast)
            sample.assembly.kmers = self.kmers

    def __init__(self, inputobject):
        from Queue import Queue
        from Bio.Blast.Applications import NcbiblastnCommandline
        from distutils import spawn
        # Find blastn and augustus version
        self.version = "v1.1b1"
        self.augustus = " ".join(Popen(['augustus', '--version'], stderr=PIPE).stderr.read().split()[:2])
        self.blast = NcbiblastnCommandline(version=True)()[0].replace('\n', ' ').rstrip()
        self.metadata = inputobject.runmetadata.samples
        self.executable = os.path.abspath(spawn.find_executable("BUSCO_{}.py".format(self.version)))
        self.start = inputobject.starttime
        self.threads = inputobject.cpus
        self.path = inputobject.path
        self.qqueue = Queue()
        printtime('Running BUSCO {} for gene discovery metrics'.format(self.version.split(",")[0]), self.start)
        # Testing with bacterial HMMs
        self.lineage = 'bacteri'
        self.buscoprocess()


