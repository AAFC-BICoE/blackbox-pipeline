#!/usr/bin/env python
from threading import Thread
from accessoryFunctions import *
from itsx.parallel import ITSx
import os

__author__ = 'mike knowles'


class ITS(object):
    def __init__(self, inputobject):
        from Queue import Queue
        self.metadata = [sample for sample in inputobject.runmetadata.samples
                         if os.path.isfile(sample.general.bestassemblyfile)]
        self.start = inputobject.starttime
        self.threads = int(inputobject.cpus)
        self.path = inputobject.path
        self.smallqueue = Queue()
        self.hmm = inputobject.hmm
        with open(which("ITSx")) as f:
            for line in f:
                if line.startswith("$app_version"):
                    self.version = line.split('=')[1].strip()[1:-2]
                    break
        printtime('Performing ITSx {} analysis'.format(self.version), self.start)
        for _ in range(len(self.metadata)):
            # Send the threads to the merge method. :args is empty as I'm using
            threads = Thread(target=self.run, args=())
            # Set the daemon to true - something to do with thread management
            threads.setDaemon(True)
            # Start the threading
            threads.start()
        for i in self.metadata:
            i.general.ITSxresults = '{}/ITSx_results'.format(i.general.outputdirectory)
            make_path(i.general.ITSxresults)
            self.smallqueue.put(i)
        self.smallqueue.join()

    @staticmethod
    def parse(sample, pos):
        sample.ITS = GenObject()
        main = lambda k, v: getattr(*k).append(v) if hasattr(*k) else setattr(*k + ([v],))
        for line in pos:
            line = line.strip()
            lst = line.split('\t')
            contig = lst[0]
            for ele in lst[2:]:
                if "No" not in ele and ': ' in ele:
                    k, v = ele.split(": ")
                    main((sample.ITS, k), "{}[{}]".format(contig, v.replace('-', ':')))

    def run(self):
        while True:
            sample = self.smallqueue.get()
            sample.software.ITSx = self.version
            positions, summary = [os.path.join(sample.general.ITSxresults, sample.name + f)
                                  for f in ['.positions.txt', '.summary.txt']]
            sample.commands.ITSx = ITSx(o=sample.general.ITSxresults,
                                        i=sample.general.bestassemblyfile,
                                        cpu=self.threads,
                                        N=2, t=self.hmm,
                                        detailed_results="T",
                                        preserve="T")
            print repr(sample.commands.ITSx)
            if not all(map(os.path.isfile, [positions, summary])):
                sample.commands.ITSx(name=sample.name, total=int(sample.assembly.TotalLength))
            if all(map(os.path.isfile, [positions, summary])):
                if os.stat([positions, summary][0]).st_size:
                    with open([positions, summary][0]) as pos:
                        self.parse(sample, pos)
            else:
                printtime("ERROR: No output generated for " + sample.name, self.start)
            self.smallqueue.task_done()


if __name__ == '__main__':
    from metadataReader import MetadataReader
    # from multiprocessing import cpu_count
    from time import time
    import json

    metadata = MetadataObject()
    metadata.path = "/data/"
    metadata.samples = [GenObject({'name': '2015-SEQ-1283'})]
    # metadata.cpus = cpu_count()
    metadata.cpus = 4
    metadata.hmm = "F,O"
    metadata.starttime = time()
    metadata.runmetadata.samples = MetadataReader(metadata).samples
    ITS(metadata)
    # metadata.runmetadata.samples[0].commands.Blast = NcbiblastnCommandline(help=True)
    # print json.dumps(dict(metadata.runmetadata.samples[0]), indent=4, sort_keys=True)
