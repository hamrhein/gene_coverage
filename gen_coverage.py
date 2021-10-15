#!/usr/bin/env python3

import argparse
import pysam
import sys

import numpy as np
import hamrhein.file_parsers.gtfparse as gtfparse


def build_parser():
    parser = argparse.ArgumentParser(prog="bam_coverage.py", description="")

    parser.add_argument(
        "-mn",
        "--min_gene_len",
        dest="min_gene_len",
        type=int,
        default=100,
        help="Set the minimum gene length",
    )

    parser.add_argument(
        "-mx",
        "--max_gene_len",
        dest="max_gene_len",
        type=int,
        help="Set the maximum gene length",
    )

    parser.add_argument(
        "-s",
        "--source",
        dest="source",
        help="Set the source filter for the GTF file",
    )

    parser.add_argument(
        "-n",
        "--normalize",
        dest="normalize",
        action="store_true",
        help="Normalize the reads",
    )

    parser.add_argument(
        "-sm",
        "--single_model",
        dest="single",
        action="store_true",
        help="Single model genes",
    )

    parser.add_argument(
        "-p",
        "--print_list",
        dest="catout",
        action="store_true",
        help="Print the output genes",
    )

    parser.add_argument(
        "-u",
        "--unique",
        dest="unique",
        action="store_true",
        help="Only use uniquely-mapped genes",
    )

    parser.add_argument(
        "-v",
        "--verbose",
        dest="verbose",
        action="store_true",
        help="Print progress and diagnostic information",
    )

    parser.add_argument("bamfile", help="Input BAM file")

    parser.add_argument("gtffile", help="Annotation file in GTF format")

    parser.add_argument("outfile", help="Output file")

    return parser


def build_gene_locus_dict(gtf_file, source=None):
    ifs = open(gtf_file, "r")
    gtf = gtfparse.GtfParse(ifs)

    gdict = dict()
    gdict_index = dict()

    for gi in gtf:
        if gi["feature"] != "exon":
            continue

        if source:
            if gi["source"] != source:
                continue

        chrom = gi["seqname"]
        left = int(gi["start"])
        right = int(gi["end"])

        gene_id = gi.gene_id
        trans_id = gi.transcript_id

        if gene_id not in gdict:
            gdict[gene_id] = dict()

        if trans_id in gdict[gene_id]:
            if chrom != gdict[gene_id][trans_id][0][0]:
                continue
        else:
            gdict[gene_id][trans_id] = list()

        gdict[gene_id][trans_id].append((chrom, left, right, gi["strand"]))

        if chrom not in gdict_index:
            gdict_index[chrom] = set()

        gdict_index[chrom].add(gene_id)

    ifs.close()

    return (gdict, gdict_index)


def build_coverage_array(bfd, reference, length, unique=False):
    arr = np.zeros(length)

    aiter = bfd.fetch(reference)

    for algn in aiter:
        if algn.is_unmapped:
            continue

        tags = dict(algn.get_tags())

        if unique and 'NH' not in tags:
            print('NH tags missing from BAM alignments')
            return

        if unique and tags['NH'] > 1:
            continue

        for b in algn.blocks:
            arr[b[0]:b[1]] += 1

    return arr


def main(args):
    if args.min_gene_len < 100:
        print("Minimum gene length set manually < 100...Setting to 100.")
        args.min_gene_len == 100

    if args.verbose:
        print("Reading and indexing annotation file.")

    gdict, gdict_idx = build_gene_locus_dict(args.gtffile, args.source)

    bfd = pysam.AlignmentFile(args.bamfile, "rb")

    name_lengths = dict(zip(bfd.references, bfd.lengths))

    gene_number = len(gdict)

    output_vector = np.zeros(100)

    for ref, length in name_lengths.items():

        if ref not in gdict_idx:
            continue

        if args.verbose:
            print("Processing {}: size {}".format(ref, length))

        coverage_arr = build_coverage_array(bfd, ref, length, args.unique)

        if coverage_arr is None:
            return 1

        for gid in gdict_idx[ref]:

            gmodel = gdict[gid]

            if args.single and len(gmodel) > 1:
                continue

            nl = set()

            for tid in gmodel:
                for (chrom, left, right, strand) in gmodel[tid]:
                    for i in range(left, right):
                        nl.add(i)
                        
            nl = sorted(nl)

            if strand == '-' or strand == 'R':
                nl.reverse()

            gene_len = len(nl)

            if gene_len < args.min_gene_len:
                gene_number -= 1
                continue

            if args.max_gene_len is not None and gene_len > args.max_gene_len:
                gene_number -= 1
                continue

            final_vector = np.zeros(100)
            bins = np.linspace(0, gene_len, num=101, dtype=int)
            start = 0

            for i, b in enumerate(bins[1:]):
                final_vector[i] = coverage_arr[nl[start:b]].mean()
                start = b

            output_vector += final_vector

    if args.normalize:
        output_vector /= len(gdict)

    with open(args.outfile, 'wt') as ofs:
        for i in range(100):
            ofs.write('{}\t{}\n'.format(i, output_vector[i]))

    return 0


if __name__ == "__main__":
    parser = build_parser()
    args = parser.parse_args()

    if len(sys.argv) == 1:
        print(parser.format_help())
        sys.exit(0)

    sys.exit(main(args))
