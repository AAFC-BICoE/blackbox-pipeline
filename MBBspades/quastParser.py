#!/usr/bin/env python
from accessoryFunctions import *
import os


__author__ = 'mikeknowles,akoziol'


class Quast(object):

    def quastprocess(self):
        from threading import Thread
        # Find the fasta files for each sample
        # Only make as many threads are there are samples with fasta files
        for i in range(len([sample.general for sample in self.metadata if type(sample.general.fastqfiles) is list])):
            # Send the threads to the merge method. :args is empty as I'm using
            threads = Thread(target=self.analyze, args=())
            # Set the daemon to true - something to do with thread management
            threads.setDaemon(True)
            # Start the threading
            threads.start()
        for sample in self.metadata:
            # Initialise the quast command
            sample.software.Quast = self.version
            if sample.general.bestassemblyfile:
                sample.general.quastresults = '{}/quast_results'.format(sample.general.outputdirectory)
                make_path(sample.general.quastresults)
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

    def analyze(self):
        """Run the quast command in a multi-threaded fashion"""
        while True:
            sample = self.qqueue.get()
            if sample.commands.Quast and not os.path.isfile('{}/report.tsv'.format(sample.general.quastresults)):
                execute(sample.commands.Quast, sample.general.quastresults)
                self.metadata(sample)
            # Signal to the queue that the job is done
            self.qqueue.task_done()

    @staticmethod
    def metadata(self, sample):
        repls = ('>=', 'Over'), ('000 Bp', 'kbp'), ('#', 'Num'), ("'", ''), ('(', ''), (')', ''), (' ', '')
        if not os.path.isfile("%s/report.tsv" % sample.general.quastresults):
            print "There was an issue getting the metadata from {0:s}".format(sample.name)
        else:
            quast = dict()
            with open("{0:s}/report.tsv".format(sample.general.quastresults)) as report:
                report.next()
                for line in report:
                    # Use headings in report as keys for the GenObject
                    k, v = [reduce(lambda a, kv: a.title().replace(*kv), repls, s) for s in line.rstrip().split('\t')]
                    quast[k] = v
            sample.assembly = GenObject(quast)
            sample.assembly.kmers = self.kmers

    def __init__(self, inputobject):
        from Queue import Queue
        # Find quast version
        from libs.qconfig import quast_version
        self.version = quast_version()
        self.metadata = inputobject.runmetadata.samples
        self.start = inputobject.starttime
        self.kmers = inputobject.kmers
        self.threads = inputobject.cpus
        self.path = inputobject.path
        self.qqueue = Queue()
        self.kmers = inputobject.kmers
        self.quastprocess()
        printtime('Running Quast {} for assembly metrics', self.start)
