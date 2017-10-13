import GeneaPy.modules.custom_exceptions as ex
from GeneaPy.modules.common import correct_hg_version
import sys
import time
import requests
import argparse
import bs4
import re
import textwrap
if not sys.platform == 'cygwin':
    import pysam

# TODO: logging

def get_seq(location, hg_version='hg19', genome=None, upstream=None, 
            downstream=None, header=True):
    ''' Return a DNA sequence from a given genomic position or range. 

    Args:
        location: genomic position ('1:10000') or genomic range ('1:10000-10100')
        hg_version: human genome version
        genome: path to genome FASTA file
        upstream: bases upstream from location
        downstream: bases downstream from location
        header: give the sequence a FASTA header
    
    Returns:
        DNA sequence representative of the given genomic
        location or range

    Notes:
        if the location variable is a position then the
        upstream and downstream variables are required to 
        return a sequence. 
    '''
    hg_version = correct_hg_version(hg_version)
    seq_range = create_region(location, upstream, downstream)
    if genome:
        seq = get_sequence_locally(seq_range, genome)
    else:
        seq = get_sequence(seq_range, hg_version)
    # capatilise location base if it is a position
    if '-' not in location:
        seq = upper_pos(seq, upstream, downstream)
    seq = textwrap.fill(seq, width=50)

    if header:
        seq = ">{} {}\n{}".format(seq_range, hg_version, seq)
    return seq

def create_region(location, upstream, downstream):
    ''' Create a genomic range which begins and ends from 
        number of bases upstream and downstream from the
        given genomic location.
    '''
    # PySam doesnt like chr preceding position/range
    location = location.replace('chr', '')
    # check location variable isn't holding a seq range
    if all(x in location for x in [':', '-']):
        return location.replace('-',',')
    chrom, pos = location.split(':')
    start_pos = int(pos) - upstream
    end_pos = int(pos) + downstream
    seq_range = chrom+":"+str(start_pos)+","+str(end_pos)
    return seq_range

def get_sequence(seq_range, hg_version):
    ''' From a genomic range and human genome version, use UCSC DAS server
        to retrieve the sequence found in the given genomic range.

        http://www.biodas.org/documents/spec-1.53.html
    '''
    req = requests.get("http://genome.ucsc.edu/cgi-bin/das/"+hg_version+
                       "/dna?segment="+seq_range.replace('-', ','))
    req.raise_for_status()
    url = bs4.BeautifulSoup(req.text, features="xml").prettify()
    search = re.findall(r"[tacg{5}].*", url)
    seqs = [s for s in search if not s.strip("tacg")] 
    seq = "".join(seqs)
    if not seq:
        raise ex.NoSequence(seq_range)
    return seq

def get_sequence_locally(seq_range, genome_path):
    ''' Get the DNA sequence of the given genomic range from 
        a locally stored genome FASTA file.
    '''
    chrom, start, end = tuple(re.split(r"[:,]", seq_range))
    chrom = "".join(("chr",chrom))
    genome = pysam.FastaFile(genome_path)
    # -1 is required otherwise the first base is missing, no idea why
    seq = genome.fetch(chrom, int(start)-1, int(end))
    return seq

def upper_pos(seq, upstream, downstream):
    ''' Capatilise the position of interest in the sequence.'''
    before = seq[:upstream]
    var = seq[upstream]
    after = seq[upstream+1:len(seq)]
    altered_seq = "".join((before.lower(), var.upper(), after.lower()))
    return altered_seq 

def get_parser():
    parser = argparse.ArgumentParser(description='scrape DNA sequence covering a given genomic position/range')
    parser.add_argument('query', type=str, help='genomic position/range')
    parser.add_argument('-hg', '--genome_version', type=str, help='human genome version (default=hg19)', default='hg19')
    parser.add_argument('-g', '--genome', type=str, help='path to FASTA genome file')
    parser.add_argument('-u', '--upstream', type=int, help='number of base upstream from genomic position')
    parser.add_argument('-d', '--downstream', type=int, help='number of base downstream from genomic position')
    parser.add_argument('-r', '--header', action='store_true', help='fasta like header')
    return parser

def cli():
    parser = get_parser()
    args = vars(parser.parse_args())
    seq = get_seq(args['query'], args['genome_version'], args['genome'], 
                  args['upstream'], args['downstream'], args['header'])
    print(seq)


if __name__ == '__main__':
    cli()
