#!/usr/bin/env python

# This script takes some mutations annotated with annovar, a table of af for a bunch of samples
# and generates a tab delimited files with these fields:
# case, AF, myid, gene, class, is_cosmic

import argparse as ap
import sys

parser = ap.ArgumentParser(description="From wide AF table and annovar annotation to long tsv with useful info to filter muts and compute models")
parser.add_argument('annovar', type=str)
parser.add_argument('af_table', type=str)
args = parser.parse_args()

annotations = {}
with open(args.annovar) as anno:
    header = anno.readline()
    for line in anno:
        line = line.rstrip()
        fields = line.split("\t")
        fields[0] = 'chr' + fields[0]
        id = ':'.join(fields[0:2]+fields[3:5])
        if not id in annotations.keys():
            annotations[id] = fields[6] + "\t" + fields[5] + "\t" + fields[8] + "\t" + fields[11]
        else:
            sys.exit('duplicated annotations!')


with open(args.af_table) as aftable:
    header = aftable.readline()
    header = header.rstrip()
    samples = header.split("\t")[1:]
    for line in aftable:
        line = line.rstrip()
        fields = line.split("\t")
        id = fields.pop(0)
        for i in range(0, len(fields)):
            if float(fields[i]) != 0:
                annot = ["\t\t.\t"]
                if id in annotations.keys():
                    annot = [annotations.get(id)]
                res = [samples[i],fields[i],id] + annot
                print("\t".join(res))
