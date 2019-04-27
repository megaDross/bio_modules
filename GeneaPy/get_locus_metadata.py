import argparse

import modules.custom_exceptions as ex
from modules.metadata import LocusMetaData


def output_all_metadata(infile, flank, outfile, hg, genome=None, gene_list=[]):
    """ Parse the metadata for all genomic positions
        detailed within infile and write to outfile
    """
    with open(outfile, "w") as out:
        write_header(out)
        with open(infile, "r") as f:
            for line in f:
                try:
                    position = line.rstrip("\n")
                    data = LocusMetaData.from_position(
                        genomic_position=position,
                        hg_version=hg,
                        flank=flank,
                        genome=genome,
                        gene_list=gene_list,
                    )
                    data_tuple = restructure_metadata(data)
                    out.write("\t".join(data_tuple) + "\n")
                except ex.NoGene as e:
                    print("ERROR: {}".format(e))
                    continue


def write_header(out):
    header = (
        "Query",
        "Genome Version",
        "Gene",
        "Gene ID",
        "Gene Location",
        "Type",
        "Transcript",
        "Exon ID",
        "Exon",
        "Intron",
        "Sequence Range",
        "Sequence",
    )
    out.write("\t".join(header) + "\n")


def restructure_metadata(data):
    query = "chr{}:{}".format(data.contig, data.position)
    gene_location = "{}:{}-{}".format(data.gene.contig, data.gene.start, data.gene.end)
    sequence = "".join(data.sequence.split("\n"))
    data_tuple = (
        query,
        data.hg_version,
        data.gene.name,
        data.gene.id,
        gene_location,
        data.gene.biotype,
        data.transcript.id,
        data.exon.id,
        data.exon.exon_no,
        data.exon.intron_no,
        data.seq_range,
        sequence,
    )
    return data_tuple


def get_parser():
    parser = argparse.ArgumentParser(
        description="Scrape a genomic positions meta-data from Ensembl"
    )
    parser.add_argument("-i", "--input", type=str, help="file with genomic positions")
    parser.add_argument("-p", "--position", type=str, help="genomic position")
    parser.add_argument(
        "-hg",
        "--genome_version",
        type=str,
        help="human genome verison (default=hg38)",
        default="hg38",
    )
    parser.add_argument(
        "-l",
        "--gene_list",
        nargs="+",
        help="list of genes to prefer if multiple genes found at the given position",
        default=[],
    )
    parser.add_argument(
        "-f",
        "--flank",
        type=int,
        help="desired number of bp to scrape from each side of the given genomic position (default=50)",
        default=50,
    )
    parser.add_argument(
        "-g", "--genome", type=str, help="path to genome FASTA file", default=None
    )
    parser.add_argument("-o", "--output", type=str, help="name of output file")
    return parser


def cli():
    parser = get_parser()
    args = vars(parser.parse_args())
    if args["input"]:
        output_all_metadata(
            args["input"],
            args["flank"],
            args["output"],
            args["genome_version"],
            args["genome"],
            args["gene_list"],
        )
    else:
        data = LocusMetaData.from_position(
            genomic_position=args["position"],
            hg_version=args["genome_version"],
            gene_list=args["gene_list"],
            flank=args["flank"],
            genome=args["genome"],
        )
        print(data)


if __name__ == "__main__":
    cli()
    # 'chr15:48733918'
