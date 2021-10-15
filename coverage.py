import pysam

import gtfparse
import gene_dict
import coverage_array

import numpy as np


class CoverageException(Exception):
    pass


class Coverage:
    def __init__(self):
        self.bam_name_lengths = None
        self.bam_object = None
        self.max_gene_length = None
        self.min_gene_length = 100
        self.output_vector = np.zeros(100)
        self.gene_dict_obj = None
        self.out_stream = None

    def __call__(self, args):
        if self.bam_object is None:
            msg = "BAM object not instantiated"
            raise CoverageException(msg)

        if self.gene_dict_obj is None:
            msg = "GtfDict object not instantiated"
            raise CoverageException(msg)

        if self.out_stream is None:
            msg = "Output stream not instantiated"
            raise CoverageException(msg)

        nlzip = zip(self.bam_object.references, self.bam_object.lengths)
        self.bam_name_lengths = dict(nlzip)

        if args.min_gene_len:
            if args.min_gene_len >= 100:
                self.min_gene_length = args.min_gene_len
            else:
                msg = "Minimum gene length must be >= 100"
                raise CoverageException(msg)

        if args.max_gene_len:
            self.max_gene_length = args.max_gene_len

        number_of_genes = 0

        for chrom, gdict in self.gene_dict_obj:
            length = self.bam_name_lengths[chrom]

            if args.verbose:
                print("Processing {}: size {}".format(chrom, length))

            coverage_arr = coverage_array.CoverageArray(length)
            coverage_arr.populate(self.bam_object, chrom, args.unique)

            number_of_transcripts = 0

            for gid, transcripts in gdict.items():

                if args.single and len(transcripts) > 1:
                    continue

                nl = set()

                for tid, exons in transcripts.items():
                    for start, end, strand in exons:
                        for i in range(start, end):
                            nl.add(i)

                    number_of_transcripts += 1

                nl = sorted(nl)

                if strand == "-" or strand == "R":
                    nl.reverse()

                gene_length = len(nl)

                if gene_length < args.min_gene_len:
                    continue

                if self.max_gene_length is not None:
                    if gene_length > self.max_gene_length:
                        continue

                final_vector = np.zeros(100)

                bins = np.linspace(0, gene_length, num=101, dtype=int)

                start = 0

                for i, b in enumerate(bins[1:]):
                    final_vector[i] = coverage_arr[nl[start:b]].mean()
                    start = b

                self.output_vector += final_vector
                number_of_genes += 1

            if args.verbose:
                print("  {} transcripts considered".format(number_of_transcripts))

        if args.normalize:
            self.output_vector /= numer_of_genes

        for i in range(100):
            self.out_stream.write("{}\t{:0.4f}\n".format(i, self.output_vector[i]))
