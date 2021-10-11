#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
import re
import sys
import logging
import argparse

LOG = logging.getLogger(__name__)

__version__ = "1.0.0"
__author__ = ("Xingguo Zhang",)
__email__ = "invicoun@foxmail.com"
__all__ = []


def read_tsv(file, sep=None):

    LOG.info("reading message from %r" % file)

    for line in open(file):
        line = line.strip()

        if not line or line.startswith("#"):
            continue

        yield line.split(sep)


def read_gfa(file):

    segment = {}
    link = {}

    for line in read_tsv(file, "\t"):
        if line[0] =='S':
            depth = "%.2f" % float(line[4].split(":")[-1])
            segment[line[1]] = (line[2], depth)
        if line[0] =='L':
            if line[1] not in link:
                link[line[1]] = []
            if line[3] not in link:
                link[line[3]] = []

            values = (line[2], line[3], line[4])
            if values not in link[line[1]]:
                link[line[1]].append(values)

            qs = "+"
            rs = "+"
            if line[4]=="+":
                qs = "-"
            if line[2]=="+":
                rs = "-"
            values = (qs, line[1], rs)
            if values not in link[line[3]]:
                link[line[3]].append(values)
    return segment, link


def get_id(file):

    r = []

    for line in read_tsv(file, "\t"):
        r.append(line[0])

    return r


def get_gfa(gfa, seqids):

    ids = get_id(seqids)
    segment, link = read_gfa(gfa)

    seqid = set()
    for i in link:
        if i not in ids:
            continue
        seqid.add(i)
        for qs, rid, rs in link[i]:
            seqid.add(rid)
            print("L\t%s\t%s\t%s\t%s\t0M" % (i, qs, rid, rs))

    for i in segment:
        if i not in seqid:
            continue
        seq, depth = segment[i]
        print("S\t%s\t%s\tLN:i:%s\tdp:f:%s" % (i, seq, len(seq), depth))

#    for i in link:
#        if i not in seqid:
#            continue
#        for qs, rid, rs in link[i]:
#            print("L\t%s\t%s\t%s\t%s\t0M" % (i, qs, rid, rs))

    return 0


def add_help(parser):

    parser.add_argument('gfa', metavar='FILE', type=str,
        help='Input the assembled gfa file.')
    parser.add_argument('-i', '--seqids', metavar='FILE', type=str, required=True,
        help='Input the id file of the extracted sequence.')

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
    get_gfa.py Extract gfa file
attention:
    get_gfa.py assembly.gfa -i id >assembly_new.gfa
version: %s
contact:  %s <%s>\
        ''' % (__version__, ' '.join(__author__), __email__))

    args = add_help(parser).parse_args()

    get_gfa(args.gfa, args.seqids)


if __name__ == "__main__":

    main()
