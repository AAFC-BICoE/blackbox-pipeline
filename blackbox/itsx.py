#!/usr/bin/env python
from Bio import SeqIO
from threading import Thread
from collections import defaultdict
from accessoryFunctions import *
import os
import shutil
__author__ = 'mike knowles'


class Bin:
    """Bin for holding sequences"""

    def __init__(self, capacity, contents=None):
        self.capacity = capacity
        self.contents = contents if contents else list()
        self.sum = sum(map(len, list()))

    def add(self, x):
        self.contents.append(x)
        self.sum += len(x)

    def __iter__(self):
        for i in self.contents:
            yield i

    def free_capacity(self):
        return self.capacity - self.sum


class BinPacker:
    """ Uses first fit algorithm to solve the bin-packing problem
    """
    def __init__(self, record, cap):
        # Extra is ideal for small contigs
        self.bins = [Bin(cap)]
        for seq_record in record:
            # Add the item to the first bin that can hold it
            # If no bin can hold it, make a new bin
            item = len(seq_record)
            for xBin in self.bins:
                if item > cap:
                    if xBin.sum == 0:
                        xBin.add(seq_record)
                        xBin.capacity = item
                        break
                    # Sometimes, large contigs are bigger than cap
                if xBin.free_capacity() >= item:
                    xBin.add(seq_record)
                    break
                if self.bins.index(xBin) == len(self.bins) - 1:
                    self.bins.append(Bin(cap))

    def __iter__(self):
        for xBin in self.bins:
            yield xBin


class ITSx(object):

    def __init__(self, inputobject):
        from Queue import Queue
        self.metadata = [sample for sample in inputobject.runmetadata.samples
                         if os.path.isfile(sample.general.bestassemblyfile)]
        self.start = inputobject.starttime
        self.kmers = inputobject.kmers
        self.threads = inputobject.cpus
        self.path = inputobject.path
        self.smallqueue = Queue()
        self.itsxqueue = Queue()

        printtime('Performing ITSx analysis', self.start)
        with open(which("ITSx")) as f:
            for line in f:
                if line.startswith("$app_version"):
                    self.version = line.split('=')[1].strip()[1:-2]
                    break
        self.pool = list()
        for _ in range(len(self.metadata)):
            # Send the threads to the merge method. :args is empty as I'm using
            threads = Thread(target=self.run, args=())
            # Set the daemon to true - something to do with thread management
            threads.setDaemon(True)
            # Start the threading
            threads.start()
            self.pool.append(threads)
            for i in range(int(self.threads)):
                itsx = Thread(target=self.parallel, args=())
                itsx.setDaemon(True)
                itsx.start()
                self.pool.append(itsx)
        for i in self.metadata:
            i.general.ITSxresults = '{}/ITSx_results'.format(i.general.outputdirectory)
            make_path(i.general.ITSxresults)
            self.smallqueue.put(i)
        self.smallqueue.join()
        for i in self.pool:
            print i.is_alive(), i.getName()
        self.


    def parallel(self):
        while True:
            execute("ITSx -i {} -o {} -t {} --debug --cpu {}".format(*self.itsxqueue.get()))
            self.itsxqueue.task_done()

    @staticmethod
    def parse(sample, pos):
        itsdict = defaultdict()
        for line in pos:
            line = line.strip()
            lst = line.split('\t')
            contig = lst.pop()
            for ele in lst[1:]:
                if "Not found" not in ele:
                    k, v = ele.split(": ")
                    itsdict[k].append("{}[{}]".format(contig, v))
        if itsdict:
            sample.ITS = GenObject(itsdict)

    def run(self):
        import math
        while True:
            sample = self.smallqueue.get()
            sample.software.ITSx = self.version
            cap = int(math.ceil(float(sample.assembly.TotalLength)/self.threads))
            cap += int(cap/2e4)
            baselist = []
            files = ['.positions.txt', '.ITS1.fasta', '.ITS2.fasta']
            finalfiles = [os.path.join(sample.general.ITSxresults, sample.name + f) for f in files]
            if not all(map(os.path.isfile, finalfiles)):
                with open(sample.general.bestassemblyfile) as bestassemblyfile:
                    record = SeqIO.parse(bestassemblyfile, "fasta")
                    for i, batch in enumerate(BinPacker(record, cap)):
                        base = os.path.join(sample.general.ITSxresults, str(i+1))
                        output = os.path.join(base, "{0:s}_group_{1:d}".format(sample.name, i + 1))
                        intermediate = [output + f for f in files + [".summary.txt"]]
                        if not all(map(os.path.isfile, intermediate)):
                            make_path(base)
                            filename = output + ".fasta"
                            with open(filename, "w") as handle:
                                SeqIO.write(list(batch), handle, "fasta")
                            output = os.path.splitext(filename)[0]
                            self.itsxqueue.put((filename, output, "F,O", self.threads))
                        baselist.append((intermediate, base))
            self.itsxqueue.join()
            for output in finalfiles:
                # Low level file i/o operation to quickly append files without significant overhead
                if hasattr(os, 'O_BINARY'):
                    o_binary = getattr(os, 'O_BINARY')
                else:
                    o_binary = 0
                output_file = os.open(output, os.O_WRONLY | o_binary | os.O_CREAT)
                for intermediate, folder in baselist:
                    input_filename = intermediate[finalfiles.index(output)]
                    input_file = os.open(input_filename, os.O_RDONLY | o_binary)
                    while True:
                        input_block = os.read(input_file, 1024 * 1024)
                        if not input_block:
                            break
                        os.write(output_file, input_block)
                    os.close(input_file)
                os.close(output_file)
            if all(map(os.path.isfile, finalfiles)):
                if os.stat(finalfiles[0]).st_size:
                    with open(finalfiles[0]) as pos:
                        self.parse(sample, pos)
                print finalfiles
                summarylines = list()
                for intermediate, folder in baselist:
                    with open(intermediate[3]) as summary:
                        if summarylines:
                            for idx, line in enumerate(summary):
                                if '\t' in line:
                                    line.split('\t')
                        else:
                            summarylines = summary.readlines()
                        for line in
                    with open(os.path.join(sample.general.ITSxresults, sample.name + '.summary.txt'), 'w+') as full:
                    # shutil.rmtree(folder)
            else:
                print_function("ERROR: No output generated for " + sample.name, self.start)
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
    metadata.starttime = time()
    metadata.runmetadata.samples = MetadataReader(metadata).samples
    ITSx(metadata)
    # metadata.runmetadata.samples[0].commands.Blast = NcbiblastnCommandline(help=True)
    # print json.dumps(dict(metadata.runmetadata.samples[0]), indent=4, sort_keys=True)
