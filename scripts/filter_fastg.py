#!/usr/bin/env python
#coding:utf-8

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


def read_fastg(file):

    '''Read fastg file'''
    if file.endswith(".gz"):
        fp = gzip.open(file)
    else:
        fp = open(file)

    seq = []
    for line in fp:
        if isinstance(line, bytes):
            line = line.decode('utf-8')
        line = line.strip()

        if not line:
            continue
        if line.startswith(">"):
            line = line.strip(">")
            if len(seq) == 2:
                yield seq
            seq = []
            seq.append(line.split()[0])
            continue
        if len(seq) == 2:
            seq[1] += line
        else:
            seq.append(line)

    if len(seq) == 2:
        yield seq
    fp.close()


def split_tag(string):

    r = {}
    key = ""
    value = ""
    for i in string.split("_"):
        try:
            key = float(i)
        except:
            value = i
        r[value] = key

    return r


def read_head(string):

    r = {}

    for ann in string.split(";"):
        for head in ann.split(":"):
            head = split_tag(head.strip("'"))
            for i in head:
                if i not in r:
                    r[i] = [head[i]]
                    continue
                r[i].append(head[i])
    return r


def filter_fastg(file, minlen=500, mincov=0.9):

    for seqid,seq in read_fastg(file):
        tags = read_head(seqid)
        if "length" in tags:
            length = max(tags["length"])
        else:
            length = len(seq)
        if "cov" in tags:
            cov = max(tags["cov"])
        else:
            mincov = 0
            cov = 1

        if length <= minlen:
            continue
        if cov <= mincov:
            continue
        print(">%s\n%s" % (seqid, seq))

    return 0


def add_hlep_args(parser):

    parser.add_argument('fastg', metavar='FILE', type=str,
        help='Input fastg file.')
    parser.add_argument('-ml', '--minlen', metavar='FLOAT', type=float, default=50,
        help='Set the minimum length of the filter, default=50.')
    parser.add_argument('-mc', '--mincov', metavar='FLOAT', type=float, default=0.5,
        help='Set the minimum coverage for filter, default=0.5.')

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
    filter_fastg.py: Filter fastg files
attention:
    filter_fastg.py assembly_graph.fastg >new_graph.fastg
version: %s
contact:  %s <%s>\
        ''' % (__version__, ' '.join(__author__), __email__))

    args = add_hlep_args(parser).parse_args()

    filter_fastg(args.fastg, args.minlen, args.mincov)


if __name__ == "__main__":

    main()
