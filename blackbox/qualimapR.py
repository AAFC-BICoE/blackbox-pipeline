#!/usr/bin/env python
from accessoryFunctions import *
from bowtie import *
from Bio.Sequencing.Applications import SamtoolsViewCommandline, SamtoolsSortCommandline
from cStringIO import StringIO
__author__ = 'mike knowles'


class QualiMap(object):

    def __init__(self, inputobject):
        from Queue import Queue
        self.bowversion = Bowtie2CommandLine(version=True)()[0].split('\n')[0].split()[-1]
        self.samversion = get_version(['samtools', '--version']).split('\n')[0].split()[1]
        self.version = get_version(['qualimap', '--help']).split('\n')[3].split()[1]
        self.metadata = inputobject.runmetadata.samples
        self.start = inputobject.starttime
        self.threads = inputobject.cpus
        self.path = inputobject.path
        self.bowqueue = Queue()
        self.qqueue = Queue()
        printtime('Aligning reads with Bowtie2 {} for Qualimap'.format(self.bowversion.split(",")[0]), self.start)
        self.bowtie()

    def bowtie(self):
        from threading import Thread
        for i in range(len([sample.general for sample in self.metadata if sample.general.bestassemblyfile])):
            # Send the threads to the merge method. :args is empty as I'm using
            threads = Thread(target=self.align, args=())
            # Set the daemon to true - something to do with thread management
            threads.setDaemon(True)
            # Start the threading
            threads.start()
        for sample in self.metadata:
            # Initialise the bowtie command and version
            sample.software.Bowtie2 = self.bowversion
            sample.software.SAMtools = self.samversion
            sagen = sample.general
            if sagen.bestassemblyfile:
                sagen.QualimapResults = '{}/qualimap_results'.format(sagen.outputdirectory)
                make_path(sagen.QualimapResults)
                sagen.bowtie2results = os.path.join(sagen.QualimapResults, sample.name)
                # Use fancy new bowtie2 wrapper
                sample.commands.Bowtie2Build = Bowtie2BuildCommandLine(reference=sagen.bestassemblyfile,
                                                                       bt2=sagen.bowtie2results)
                sbam = sagen.bowtie2results + ".sorted.bam"
                sample.commands.SAMtools = "samtools view -bS - | samtools sort -o {} -".format(sbam)
                if len(sagen.fastqfiles) == 2:
                    indict = dict(("m"+str(x+1), sagen.fastqfiles[x]) for x in range(2))
                else:
                    indict = dict(("U", ",".join(sagen.fastqfiles)))
                sample.commands.Bowtie2Align = Bowtie2CommandLine(bt2=sagen.bowtie2results,
                                                                  threads=self.threads,
                                                                  **indict)
                self.bowqueue.put((sample, sbam))
            else:
                sample.commands.Bowtie2Align = "NA"
                sample.commands.Bowtie2Build = "NA"
                sample.commands.SAMtools = "NA"
        self.bowqueue.join()

    def align(self):
        while True:
            sample, sbam = self.bowqueue.get()
            if not os.path.isfile(sbam) and sample.general.bestassemblyfile != 'NA':
                stdout = StringIO()
                for func in sample.commands.Bowtie2Build, sample.commands.Bowtie2Align:
                    stdout.close()
                    # Use cStringIO streams to handle bowtie output
                    stdout, stderr = map(StringIO, func())
                    if stderr:
                        # Write the standard error to log, bowtie2 puts alignmentsummary here
                        with open(os.path.join(sample.general.QualimapResults, "bowtie.log"), "ab+") as log:
                            log.write("{}\n{}\n{}".format(func, stderr.getvalue(), '-'*60))
                    stderr.close()
                # stdout will be the SAM file from alignment
                if stdout:
                    # PIPE stdout to stdin of samtools view then sort (only outputing sorted bam)
                    # SAMtools sort v1.3 has different run parameters
                    if self.samversion < "1.3":
                        samsort = SamtoolsSortCommandline(input_bam="-", out_prefix=sbam[:-4])
                    else:
                        samsort = SamtoolsSortCommandline(input_bam=sbam, o=True, out_prefix="-")
                    # Use cStringIO streams to handle bowtie output
                    for func in [SamtoolsViewCommandline(b=True, S=True, input_file="-"), samsort]:
                        stdout, stderr = map(StringIO, func(stdin=stdout.getvalue()))
                        stdout.close()
                        # Write the standard error to log
                        with open(os.path.join(sample.general.QualimapResults, "samtools.log"), "ab+") as log:
                            log.write("{}\n{}\n{}".format(func, stderr.getvalue, '-'*60))
                        stderr.close()
            # Signal to the queue that the job is done
            self.bowqueue.task_done()

    def __call__(self, *args, **kwargs):
        """Exectute Qualimap on call"""
        from threading import Thread
        printtime('Reading BAM file for Qualimap output', self.start)
        for i in range(len([sample.general for sample in self.metadata if sample.general.bestassemblyfile])):
            # Send the threads to the merge method. :args is empty as I'm using
            threads = Thread(target=self.mapper, args=())
            # Set the daemon to true - something to do with thread management
            threads.setDaemon(True)
            # Start the threading
            threads.start()
        for sample in self.metadata:
            sbam = os.path.join(sample.general.QualimapResults, sample.name + ".sorted.bam")
            sample.software.Qualimap = self.version
            if sample.general.bestassemblyfile:
                sample.commands.Qualimap = 'qualimap bamqc -bam {} -outdir {}'.\
                    format(sbam, sample.general.QualimapResults)
                self.qqueue.put(sample)
            else:
                sample.commands.Qualimap = "NA"
        self.qqueue.join()

    def mapper(self):
        while True:
            sample = self.qqueue.get()
            log = os.path.join(sample.general.QualimapResults, "qualimap.log")
            reportfile = os.path.join(sample.general.QualimapResults, 'genome_results.txt')
            if not os.path.isfile(reportfile):
                execute(sample.commands.Qualimap, log)
            qdict = dict()
            if os.path.isfile(reportfile):
                with open(reportfile) as report:
                    for line in report:
                        key, value = self.analyze(line)
                        if (key,value) != (None, None):
                            qdict[key] = value
            if qdict:
                # Make new category for Qualimap results
                sample.mapping = GenObject(qdict)
            self.qqueue.task_done()

    def analyze(self, line):
        if ' = ' in line:
            key, value = line.split(' = ')
            key = key.replace('number of ', "").replace("'", "").title().replace(" ", "")
            # Should we keep comma seperation?
            value = value.replace(",", "").replace(" ", "").rstrip()
        else:
            key, value = None, None
        return key, value

if __name__ == '__main__':
    from metadataReader import MetadataReader
    from multiprocessing import cpu_count
    from time import time
    import json
    metadata = MetadataObject()
    metadata.path = "/data/"
    metadata.samples = [GenObject({'name': '2015-SEQ-1283'})]
    # metadata.cpus = cpu_count()
    metadata.cpus = 4
    metadata.starttime = time()
    metadata.runmetadata.samples = MetadataReader(metadata).samples
    QualiMap(metadata)()
    print json.dumps(dict(metadata.runmetadata.samples[0]), indent=4, sort_keys=True)
