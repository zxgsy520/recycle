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


def read_tsv(file, sep=None):

    LOG.info("Reading message from %r" % file)

    for line in open(file):
        line = line.strip()

        if not line or line.startswith("#"):
            continue

        yield line.split(sep)


def read_stdin(sep=None):

    for line in sys.stdin:
        line = line.strip()

        if not line or line.startswith("#"):
            continue

        yield line.split(sep)


def filter_paf(file, score=0, platform="mgi"):

    reads = set()

    if file:
        fh = read_tsv(file, '\t')
    else:
        fh = read_stdin('\t')

    for line in fh:
        if int(line[11]) <= score:
            continue
        seqid = line[0]
        if platform=="mgi":
            seqid = seqid.split('/')[0]
            reads.add("%s/1" % seqid)
        else:
            reads.add(seqid)

    return reads


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
        elif line.startswith('@') and len(seq)==0:
            seq.append(line.strip("@").split()[0])
        else:
            seq.append(line)

        if len(seq)==4:
            yield seq
            seq = []

    fp.close()


def stat_gc(string):

    gc = string.upper().count('G') + string.upper().count('C')
    at = string.upper().count('A') + string.upper().count('T')

    return gc*100.0/(gc+at)


def choose_read1(file, gc, readid, base, prefix="out", platform="mgi"):

    reads = set()
    n = 0

    fo = open("%s.choose_r1.fq" % prefix, "w")
    for line in read_fastq(file):
        if (line[0] not in readid) and (stat_gc(line[1]) >= gc):
            continue
        if platform == "mgi":
            seqid = line[0].split('/')[0]
            reads.add("%s/2" % seqid)
        else:
            reads.add(line[0])
        fo.write("@%s\n" % "\n".join(line))
        n += len(line[1])

        if n >= base:
            break
    fo.close()

    return reads


def gmk2pb(string):

    string = string.lower().strip()

    if string.endswith('g') or string.endswith('gb'):
        base = string.split('g')[0]
        base = float(base)*1e+09
    elif string.endswith('m') or string.endswith('mb'):
        base = string.split('m')[0]
        base = float(base)*1e+06
    elif string.endswith('k') or string.endswith('kb'):
        base = string.split('k')[0]
        base = float(base)*1e+03
    else:
        base = float(string)

    return int(base)


def choose_ml_map(file, read1, read2, base="all", prefix="out",
                     gc=0, score=0, platform="mgi"):

    readids = filter_paf(file, score, platform)

    if base == "all":
        base = float("inf")
    else:
        base = gmk2pb(base)*1.0/2

    readids = choose_read1(read1, gc, readids, base, prefix, platform)

    if read2:
        fo = open("%s.choose_r2.fq" % prefix, "w")
        for line in read_fastq(read2):
            if line[0] not in readids:
                continue
            fo.write("@%s\n" % "\n".join(line))

        fo.close()
    else:
        LOG.info("Did you enter third-generation sequencing?")

    return 0


def add_hlep_args(parser):

    parser.add_argument("-i", "--input", metavar='FILE', type=str, default="",
        help="Input files in paf.")
    parser.add_argument("-r1", "--read1", metavar='FILE', type=str, required=True,
        help="Input sequencing data R1(fastq).")
    parser.add_argument("-r2", "--read2", metavar='FILE', type=str, default="",
        help="Input sequencing data R2(fastq).")
    parser.add_argument("--platform", choices=["illumina", "mgi"], default="illumina",
        help="Input the platform for second-generation sequencing, default=illumina")
    parser.add_argument("-b", "--base", metavar='STR', type=str, default='all',
        help="Selected data volume, default=all")
    parser.add_argument("--gc", metavar='INT', type=int, default=30,
        help="Minimum GC content of reads, default=30 ")
    parser.add_argument("-s", "--score", metavar='INT', type=int, default=0,
        help="Match score, default=0")
    parser.add_argument("-p", "--prefix", metavar='STR', type=str, default='out',
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
    choose_ml_map Extract possible mitochondrial sequences based on alignment quality and GC.

attention:
    choose_ml_map -i input.paf -r1 R1.fastq -r2 R2.fastq  -p mat
    choose_ml_map -i input.paf -r1 R1.fastq -r2 R2.fastq --platform  mgi -p mat
    choose_ml_map -i input.paf -r1 R1.fastq -p mat
    minimap2 -t 10 -x sr ref.fa R1.fastq R2.fastq | choose_ml_map -r1 R1.fastq -r2 R2.fastq  -p mat 

version: %s
contact:  %s <%s>\
        ''' % (__version__, ' '.join(__author__), __email__))

    args = add_hlep_args(parser).parse_args()

    choose_ml_map(args.input, args.read1, args.read2, args.base, args.prefix,
                     args.gc, args.score, args.platform)


if __name__ == "__main__":

    main()
