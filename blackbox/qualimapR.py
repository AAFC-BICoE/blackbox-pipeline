#!/usr/bin/env python
from accessoryFunctions import *
from bowtie import *
from Bio.Sequencing.Applications import SamtoolsViewCommandline
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
            sample.software.Bowtie2 = self.bowversion
            sagen = sample.general
            if sagen.bestassemblyfile:
                sagen.QualimapResults = '{}/qualimap_results'.format(sagen.outputdirectory)
                make_path(sagen.QualimapResults)
                sagen.bowtie2results = os.path.join(sagen.QualimapResults, sample.name)
                sample.commands.Bowtie2Build = Bowtie2BuildCommandLine(reference=sagen.bestassemblyfile,
                                                                       bt2=sagen.bowtie2results)
                if len(sagen.fastqfiles) == 2:
                    indict = dict(("m"+str(x+1), sagen.fastqfiles[x]) for x in range(2))
                else:
                    indict = dict(("U", ",".join(sagen.fastqfiles)))
                sample.commands.Bowtie2Align = Bowtie2CommandLine(bt2=sagen.bowtie2results, **indict)
                self.bowqueue.put((sample.commands.Bowtie2Build, sample.commands.Bowtie2Align))
            else:
                sample.commands.Quast = "NA"
        self.bowqueue.join()

    def align(self):
        while True:
            sample = self.bowqueue.get()
            for func in sample.commands.Bowtie2Build, sample.commands.Bowtie2Align:
                stdout, stderr = func()
                if stderr:
                    # Write the standard error to log, bowtie2 puts alignmentsummary here
                    with open(os.path.join(sample.general.QualimapResults, "bowtie.log"), "ab+") as log:
                        log.write("{}\n{}\n{}".format(func, stderr, '-'*60))
            # Signal to the queue that the job is done
            self.bowqueue.task_done()


    def build(self):
        """Run the quast command in a multi-threaded fashion"""
        while True:
            sample = self.bowqueue.get()
            if sample.general.bestassemblyfile != 'NA' \
                    and not os.path.isfile('{}/report.tsv'.format(sample.general.quastresults)):
                log = os.path.join(sample.general.quastresults, 'stdout.log')
                # sys.stdout, sys.stderr = out, err
                execute(sample.commands.Quast)
            if os.path.isfile('{}/report.tsv'.format(sample.general.quastresults)):
                self.metaparse(sample)
            # Signal to the queue that the job is done
            self.bowqueue.task_done()

    def __init__(self, inputobject):
        from Queue import Queue
        self.bowversion = Bowtie2CommandLine(version=True)()[0].split('\n')[0].split()[-1]
        self.metadata = inputobject.runmetadata.samples
        self.start = inputobject.starttime
        self.threads = inputobject.cpus
        self.path = inputobject.path
        self.bowqueue = Queue()
        printtime('Aligning reads with Bowtie2 {} for qualimap'.format(self.bowversion.split(",")[0]), self.start)


if __name__ == '__main__':
    pass