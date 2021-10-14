#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
import re
import sys
import gzip
import logging
import argparse

LOG = logging.getLogger(__name__)

__version__ = "1.1.0"
__author__ = ("Xingguo Zhang",)
__email__ = "invicoun@foxmail.com"
__all__ = []


def read_paf(file, coverage, score):

    for line in open(file):
        line = line.strip()

        if not line or line.startswith("#"):
            continue

        line = line.split('\t')
        qlen = int(line[1])
        match = int(line[9])
        blen = int(line[10])
        match = min(blen, match)
        
        if ((blen*100.0/qlen)>=coverage and (match*100.0/blen)>=coverage) or int(line[11])>score:
            yield line


def read_fasta(file):
    '''Read fasta file'''

    if file.endswith(".gz"):
        fp = gzip.open(file)
    elif file.endswith(".fasta") or file.endswith(".fa"):
        fp = open(file)
    else:
        raise Exception("%r file format error" % file)

    seq = ''

    for line in fp:
        if isinstance(line, bytes):
            line = line.decode('utf-8')
        line = line.strip()

        if not line:
            continue
        if not seq:
            seq += "%s\n" % line.strip(">").split()[0]
            continue
        if line.startswith(">"):
            line = line.strip(">").split()[0]
            seq = seq.split('\n')

            yield [seq[0], seq[1]]
            seq = ''
            seq += "%s\n" % line
        else:
            seq += line

    seq = seq.split('\n')
    if len(seq)==2:
        yield [seq[0], seq[1]]
    fp.close()


def read_fastq(file):
    '''Read fastq file'''

    if file.endswith(".gz"):
        fp = gzip.open(file)
    elif file.endswith(".fastq") or file.endswith(".fq"):
        fp = open(file)
    else:
        raise Exception("%r file format error" % file)

    seq = []

    for line in fp:
        if isinstance(line, bytes):
            line = line.decode('utf-8')
        line = line.strip()

        if not line:
            continue
        if line.startswith('@') and (len(seq)==0 or len(seq)>=5):
            seq.append(line.strip("@").split()[0])
            continue
        if line.startswith('@') and len(seq)==4:
            yield seq
            seq = []
            seq.append(line.strip("@").split()[0])
        else:
            seq.append(line)

    if len(seq)==4:
        yield seq
    fp.close()


def stat_gc(string):

    gc = string.upper().count('G') + string.upper().count('C')
    at = string.upper().count('A') + string.upper().count('T')

    return gc*100.0/(gc+at)


def choose_reads(file, gc, readid, name, rn, base, note):

    if base=="all":
        nbase = "all"
        base = 0
    else:
        nbase = ""
        base = int(base)/2

    n = 0
    output = open("%s.choose_r%s.fq" % (name, rn), "w")
    
    if file.endswith(".fastq") or file.endswith(".fq") or file.endswith(".fastq.gz") or file.endswith(".fq.gz"):
        fh = read_fastq(file)
    else:
        raise Exception("%r file format error" % file)

    for line in fh:
        if rn==1 and stat_gc(line[1])<gc:
            readid.add(line[0])
        if line[0] not in readid:
            continue
        if n>=base and nbase!="all":
            break
        output.write('@%s %s\n%s\n%s\n%s\n' % (line[0], note, line[1], line[2], line[3]))
        n+=len(line[1])
    output.close()


def choose_map_reads(read1, read2, pafile, platform, name, base, gc, coverage, score):

    read = set()

    for line in read_paf(pafile, coverage, score):
        if platform=="mgi":
            ids = line[0].split('/')[0]
            read.add("%s/1" % ids)
            read.add("%s/2" % ids)
        else:
            read.add(line[0])
    note = "1:N:0"
    choose_reads(read1, gc, read, name, 1, base, note)
    note = "2:N:0"
    choose_reads(read2, gc, read, name, 2, base, note)


def add_hlep_args(parser):

    parser.add_argument("-i", "--input", metavar='FILE', type=str, required=True,
        help="Input files in paf")
    parser.add_argument("-r1", "--read1", metavar='FILE', type=str, required=True,
        help="Input files in fastq format")
    parser.add_argument("-r2", "--read2", metavar='FILE', type=str, required=True,
        help="Input files in fastq format")
    parser.add_argument("-p", "--platform", choices=["illumina", "mgi"], default="illumina",
        help="Input the platform for second-generation sequencing, default=illumina")
    parser.add_argument("-b", "--base", metavar='STR', type=str, default='all',
        help="Selected data volume, default=all")
    parser.add_argument("--gc", metavar='INT', type=int, default=30,
        help="Minimum GC content of reads, default=30 ")
    parser.add_argument("-c", "--coverage", metavar='INT', type=int, default=60,
        help="Match coverage, default=60 ")
    parser.add_argument("-s", "--score", metavar='INT', type=int, default=0,
        help="Match score, default=0")
    parser.add_argument("-n", "--name", metavar='STR', type=str, default='out',
        help="The name of the output file")
    
    return parser


def main():

    logging.basicConfig(
        stream=sys.stderr,
        level=logging.INFO,
        format="[%(levelname)s] %(message)s"
    )
    parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter,
        description='''
name:
    choose_map_ngs.py  Select no reads on the map

attention:
    choose_map_ngs.py -i paf -r1 R1.fastq -r2 R2.fastq  -n unmap

version: %s
contact:  %s <%s>\
        ''' % (__version__, ' '.join(__author__), __email__))

    args = add_hlep_args(parser).parse_args()

    choose_map_reads(args.read1, args.read2, args.input, args.platform, args.name, args.base, args.gc, args.coverage, args.score)


if __name__ == "__main__":

    main()
