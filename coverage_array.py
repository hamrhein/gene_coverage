import numpy as np


class CoverageArrayException(Exception):
    pass


class CoverageArray:
    def __init__(self, length):
        self.arr = np.zeros(length)

    def __getitem__(self, index):
        return self.arr[index]

    def populate(self, bam_object, reference, unique=False):
        align_iter = bam_object.fetch(reference)

        for algn in align_iter:
            if algn.is_unmapped:
                continue

            tags = dict(algn.get_tags())

            if unique and "NH" not in tags:
                msg = 'NH tags missing from BAM alignments'
                raise CoverageArrayException(msg)

            if unique and tags['NH'] > 1:
                continue

            for b in algn.blocks:
                self.arr[b[0]:b[1]] += 1

