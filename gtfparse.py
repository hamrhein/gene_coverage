
class GtfParseException(Exception):
    pass


class GtfItem:
    def __init__(self, fields):
        self.fields = fields

    def __str__(self):
        return str(self.fields)

    def __repr__(self):
        return "<GtfItem object at 0x{0:x}>".format(id(self))

    def __getitem__(self, attr):
        if attr in self.fields:
            return self.fields[attr]
        elif attr in self.fields["attributes"]:
            return self.fields["attributes"][attr]
        else:
            raise KeyError

    @property
    def gene_id(self):
        return self.fields["attributes"]["gene_id"]

    @property
    def transcript_id(self):
        return self.fields["attributes"]["transcript_id"]

    @property
    def attributes(self):
        return self.fields["attributes"]


class GtfParse:
    attrkeys = ["seqname", "source", "feature", "start",
            "end", "score", "strand", "frame"]

    def __init__(self, iterable):
        self.iterable = iterable
        self.encoding = "utf-8"

    def __iter__(self):
        return self

    def __next__(self):
        while True:
            try:
                line = next(self.iterable)
            except:
                raise StopIteration

            if isinstance(line, bytes):
                line = line.decode(self.encoding)

            ctok = line.split("#", 1)

            if len(ctok[0]) > 0:
                break

        tok = ctok[0].split("\t")
        atok = tok[8].replace('"', '').lstrip().split(";")[:-1]
        fields = dict(zip(self.attrkeys, tok[:-1]))
        fields["attributes"] = dict(a.lstrip().split(" ", 1) for a in atok)

        if len(ctok) > 1:
            fields["comments"] = ctok[1]
        else:
            fields["comments"] = ''

        return GtfItem(fields)

