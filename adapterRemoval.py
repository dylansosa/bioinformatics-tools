# -*- coding: utf-8 -*-
# !/usr/bin/env python2.7
# @author dylan sosa
# December 8 2017
# BCB 5200
# Dr. Ahn
import os, sys, re
from Bio.Seq import Seq
from Bio import SeqIO


def trim_adapters(records, start_adapter, end_adapter, min_length):
    for record in records:
        if len(record) < min_length:
            continue
        elif record.seq.startswith(start_adapter) and record.seq.endswith(end_adapter):
            index = len(record) - (len(start_adapter) + len(end_adapter))
            if(index >= min_length):
                yield record[len(start_adapter):-len(end_adapter)]
        elif record.seq.startswith(start_adapter): # trimming primer
            index = len(record) - len(start_adapter)
            if(index >= min_length):
                yield record[len(start_adapter):]
        elif record.seq.endswith(end_adapter): # trimming primer
            index = len(record) - len(end_adapter)
            if(index >= min_length):
                yield record[:-len(end_adapter)]
        else:
            yield record

infile = open(sys.argv[1],"r")
outfile = open(sys.argv[2], "w")
start_adapter = sys.argv[3]
end_adapter = sys.argv[4]
min_length = int(sys.argv[5])

original_reads = SeqIO.parse(infile, "fastq")
primer_trimmed_reads = trim_adapters(original_reads, start_adapter,end_adapter,min_length)
count = SeqIO.write(primer_trimmed_reads, outfile, "fastq")
print("Saved %i reads" % count)
