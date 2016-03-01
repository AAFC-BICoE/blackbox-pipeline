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

        printtime('Assembling sequences', self.start)
        with open(which("ITSx")) as f:
            for line in f:
                if f.startswith("$app_version"):
                    self.version = line.split('=')[1].strip()[1:-2]
                    break
        for i in self.metadata:
            # Send the threads to the merge method. :args is empty as I'm using
            threads = Thread(target=self, args=())
            itsx = Thread(target=self.parallel, args=())
            # Set the daemon to true - something to do with thread management
            threads.setDaemon(True)
            itsx.setDaemon(True)
            # Start the threading
            threads.start()
            itsx.start()
        for i in self.metadata:
            i.general.ITSxresults = '{}/ITSx_results'.format(i.general.outputdirectory)
            make_path(i.general.ITSxresults)
            self.smallqueue.put(i)
        self.smallqueue.join()


    def parallel(self):
        while True:
            execute("ITSx -i {} -o {} -t {} --debug --cpu {}".format(*self.itsxqueue.get()))
            self.itsxqueue.task_done()

    @staticmethod
    def parse(sample):
        itsdict = defaultdict()
        with open(os.path.join(sample.general.ITSxresults, sample.name + ".positions.txt")) as pos:
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

    def __call__(self):
        import math
        while True:
            sample = self.smallqueue.get()
            sample.software.ITSx = self.version
            cap = int(math.ceil(float(sample.assembly.TotalLength)/self.threads))
            cap /= 2e4
            baselist = []
            files = ['.positions.txt', '.ITS1.fasta', '.ITS2.fasta']
            with open(sample.general.bestassemblyfile) as record:
                for i, batch in enumerate(BinPacker(record, cap)):
                    base = os.path.join(sample.general.ITSxresults, i)
                    make_path(base)
                    filename = os.path.join(base, "%s_group_%i.fasta" % (sample.name ,i+1))
                    handle = open(filename, "w")
                    SeqIO.write(list(batch), handle, "fasta")
                    handle.close()
                    output = os.path.splitext(filename)[0]
                    baselist.append((output, base))
                    self.itsxqueue.put((filename, output, "F,O", self.threads))
            self.itsxqueue.join()
            for f in files:
                output = os.path.join(sample.general.ITSxresults, sample.name + f)
                if hasattr(os, 'O_BINARY'):
                    o_binary = getattr(os, 'O_BINARY')
                else:
                    o_binary = 0
                output_file = os.open(output, os.O_WRONLY | o_binary)
                for input_filename, folder in [[x+f, fold] for x, fold in baselist]:
                    input_file = os.open(input_filename, os.O_RDONLY | o_binary)
                    while True:
                        input_block = os.read(input_file, 1024 * 1024)
                        if not input_block:
                            break
                        os.write(output_file, input_block)
                    os.close(input_file)
                    shutil.rmtree(folder)
                os.close(output_file)
            self.parse(sample)
            self.smallqueue.task_done()


if __name__ == '__main__':
    pass
    # import math
    # from Bio import SeqIO
    # # from time import clock
    # # # record = list(SeqIO.parse("/Users/mike/database/testset/BestAssemblies/WG11-0646B-BC06.fasta", "fasta"))
    # record = SeqIO.parse("/Users/mike/sshfs/test/BestAssemblies/2015-SEQ-1283.fasta", "fasta")
    # # record.close()
    # # t1 = clock()
    # # items = 0
    # #
    # # # Extra is ideal for small contigs
    # # i =0
    # cap = int(math.ceil(float(2977468)/5)) + int(math.ceil(float(2977468)/1e5))
    # base = "/Users/mike/test"
    # for i, batch in enumerate(BinPacker(record, cap)):
    #     # print i, batch.free_capacity(), batch.capacity
    #     filename = os.path.join(base, "group_%i.fasta" % (i+1))
    #     handle = open(filename, "w")
    #     count = SeqIO.write(list(batch), handle, "fasta")
    #     handle.close()
    #     print "Wrote %i records to %s" % (count, filename)
    # print i
    # bins = []
    # bins.append(Bin(cap))
    # lst = list()
    # for seq_record in record:
    #     # Add the item to the first bin that can hold it
    #     # If no bin can hold it, make a new bin
    #     item = len(seq_record)
    #     items += item
    #     # cap = item if item > cap else cap
    #     for xBin in bins:
    #         if item > cap:
    #             if xBin.sum == 0:
    #                 print item
    #                 xBin.add(seq_record)
    #                 xBin.capacity = item
    #                 print xBin.free_capacity()
    #                 break
    #             # Sometimes, large contigs are bigger than cap
    #         if xBin.free_capacity() >= item:
    #             xBin.add(seq_record)
    #             break
    #         if bins.index(xBin) == len(bins) - 1:
    #             bins.append(Bin(cap))
    # remaining = [xBin.free_capacity() for xBin in bins]
    # capacity_left = sum(remaining)
    # capacity_total = sum([xBin.capacity for xBin in bins])
    # capacity_used = capacity_total - capacity_left
    # average = capacity_left / float(len(bins))
    # sd = (sum((x-average)**2 for x in remaining)/len(remaining))**0.5
    # print "Algorithm runtime (s):", clock()-t1
    # print "First-fit Decreasing algorithm for", items, "with capacity", cap, "used", len(bins), "bins"
    # # print "The configuration was", bins
    # print "Capacity remaining per bin: ", remaining
    # print "Average capacity remaining per bin: ", average, sd
    # print "Total capacity used:", capacity_used
    # print "Total capacity remaining: ", capacity_left
    # print "Efficiency (%): ", 100 * (1 - capacity_left/float(capacity_total))
