#!/usr/bin/env python
from accessoryFunctions import *
from bowtie import *
__author__ = 'mike knowles'


class QualiMap(object):

    def bowtie(self):
        from threading import Thread
        for i in range(len([sample.general for sample in self.metadata if sample.general.bestassemblyfile])):
            # Send the threads to the merge method. :args is empty as I'm using
            threads = Thread(target=self.build(), args=())
            # Set the daemon to true - something to do with thread management
            threads.setDaemon(True)
            # Start the threading
            threads.start()
        for sample in self.metadata:
            # Initialise the bowtie command and version
            sample.software.Bowtie2 = self.version
            sagen = sample.general
            if sagen.bestassemblyfile:
                sagen.qualimapresults = '{}/qualimap_results'.format(sagen.outputdirectory)
                sagen.bowtie2results = os.path.join(sagen.qualimapresults, sample.name)
                sample.commands.Bowtie2Build = Bowtie2BuildCommandLine(reference=sagen.bestassemblyfile,
                                                                       bt2=sagen.bowtie2results)
                if len(sagen.fastqfiles) == 2:
                    indict = dict(("m"+str(x+1), sagen.fastqfiles[x]) for x in range(2))
                else:
                    indict = dict(("U", ",".join(sagen.fastqfiles)))
                sample.commands.Bowtie2Align = Bowtie2CommandLine(bt2=sagen.bowtie2results, **indict)
                self.qqueue.put((sample.commands.Bowtie2Build, sample.commands.Bowtie2Align))
            else:
                sample.commands.Quast = "NA"
        self.qqueue.join()

    def align(self):
        while True:
            index, bwa = self.qqueue
            index()
            bwa()
            # Signal to the queue that the job is done
            self.qqueue.task_done()


    def build(self):
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

    def __init__(self, inputobject):
        from Queue import Queue
        self.bowversion = Bowtie2CommandLine(version=True)()[0].split('\n')[0].split()[-1]
        self.metadata = inputobject.runmetadata.samples
        self.start = inputobject.starttime
        self.threads = inputobject.cpus
        self.path = inputobject.path
        self.qqueue = Queue()
        printtime('Aligning reads with Bowtie2 {} for qualimap'.format(self.version.split(",")[0]), self.start)


if __name__ == '__main__':
    pass