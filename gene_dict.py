class GtfDictException:
    pass


class GtfDict:
    def __init__(self, parser_obj, source=None):
        self.gdict = dict()
        self.total_genes = 0
        self.total_transcripts = 0

        for gi in parser_obj:
            if gi['feature'] != 'exon':
                continue

            if source:
                if gi['source'] != source:
                    continue

            chrom = gi['seqname']
            
            exon = (int(gi['start']), int(gi['end']), gi['strand'])

            if chrom not in self.gdict:
                self.gdict[chrom] = dict()

            if gi.gene_id not in self.gdict[chrom]:
                self.gdict[chrom][gi.gene_id] = dict()
                self.total_genes += 1

            if gi.transcript_id not in self.gdict[chrom][gi.gene_id]:
                self.gdict[chrom][gi.gene_id][gi.transcript_id] = list()
                self.total_transcripts += 1

            self.gdict[chrom][gi.gene_id][gi.transcript_id].append(exon)

    def __iter__(self):
        return iter(self.gdict.items())
