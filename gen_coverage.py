#!/usr/bin/env python3

import argparse
import sys

import pysam

import gene_dict
import coverage_array
import coverage

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


def main(args):
    coverage_obj = coverage.Coverage()

    try:
        gtf_fs = open(args.gtffile, 'rt')
    except Exception as e:
        print('Failed to open GTF file: {}'.format(str(e)))
        return 1

    gtf = gtfparse.GtfParse(gtf_fs)

    if args.verbose:
        print('Reading GTF file...', end='', flush=True)

    coverage_obj.gene_dict_obj = gene_dict.GtfDict(gtf)

    if args.verbose:
        print('done.')

    try:
        bam_object = pysam.AlignmentFile(args.bamfile, 'rb')
    except Exception as e:
        print('Failed to open BAM file: {}'.format(str(e)))
        gtf_fs.close()
        return 1

    coverage_obj.bam_object = bam_object

    try:
        out_stream = open(args.outfile, 'wt')
    except Exception as e:
        print('Failed to open output stream: {}'.format(str(e)))
        gtf_fs.close()
        bam_object.close()
        return 1

    coverage_obj.out_stream = out_stream

    coverage_obj(args)
    
if __name__ == '__main__':
    parser = build_parser()
    args = parser.parse_args()

    if len(sys.argv) == 1:
        print(parser.format_help())
        sys.exit(0)

    sys.exit(main(args))
