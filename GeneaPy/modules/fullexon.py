from pyensembl.exon import Exon


class FullExon(Exon):
    """ Enhances the Exon object by adding exon/intron numbers"""

    def __init__(
        self,
        exon_id,
        contig,
        start,
        end,
        strand,
        gene_name,
        gene_id,
        position,
        number,
        exon,
    ):
        Exon.__init__(self, exon_id, contig, start, end, strand, gene_name, gene_id)
        self.position = position
        self.number = number
        self.exon = exon

    def __eq__(self, other):
        return self.__dict__ == other.__dict__

    def __str__(self):
        name = "Exon" if self.exon else "Intron"
        return (
            name
            + "(exon_id=%s, gene_name=%s, contig=%s, start=%d, end=%s, position=%s, number=%s)"
            % (
                self.id,
                self.gene_name,
                self.contig,
                self.start,
                self.end,
                self.position,
                self.number,
            )
        )

    @classmethod
    def from_pyexon(cls, Exon, position, number, exon):
        """ Parse a pyensembl Exon object instead"""
        return cls(
            Exon.id,
            Exon.contig,
            Exon.start,
            Exon.end,
            Exon.strand,
            Exon.gene_name,
            Exon.gene_id,
            position,
            number,
            exon,
        )
