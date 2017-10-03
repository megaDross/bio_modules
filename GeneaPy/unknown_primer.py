import argparse
import requests 
import logging
import re 
import bs4
from modules import metadata
import modules.custom_exceptions as ex
from modules.common import correct_hg_version

logging.basicConfig(filename='unknown_primer.error.log', 
                    level=logging.ERROR,
                    format="%(asctime)s:%(levelname)s:%(message)s")

# TODO: logging issues
# TODO: unit testing

def unknown_primer(f_primer, r_primer, hg_version, primer_name,
                   max_size, min_perfect, min_good):
    ''' Use primer pairs to scrape in-silico PCR amplicon sequences
        from UCSC and gene/exon data from Ensembl.

    Args:
        primer_name: primer pair name
        f_primer: forward primer sequence
        r_primer: reverse primer sequence
        hg_version: human genome version
        max_size: maximum resulting amplicon size
        min_perfect: no. of bases that match exactly on 3' end of primers
        min_good: no. of bases on 3' end of primers where at least 2 out of 3 bases match 

    Returns:
        The in-silico generated amplicons metadata.
    ''' 
    hg_version = correct_hg_version(hg_version)
    check_input_errors(primer_name, f_primer, r_primer, hg_version)
    data = scrape_seq(primer_name, f_primer, r_primer, hg_version, 
                      max_size, min_perfect, min_good)
    header, seq = seperate_data(data)
    locus_metadata = get_metadata(header, seq, hg_version)
    all_data = (primer_name, f_primer, r_primer, hg_version) + locus_metadata
    return all_data

def check_input_errors(primer_name, f_primer, r_primer, hg_version):
    ''' Catch common input errors and alert the user'''
    if re.search(r'[^ATCG]',f_primer):
        raise ex.AmbigousBase(primer_name, f_primer)
    if re.search(r'[^ATCG]',r_primer):
        raise ex.AmbigousBase(primer_name, r_primer)
    valid = ['hg38', 'hg19', 'hg18', 'hg17', 'hg16']
    if hg_version not in valid: 
        raise ex.WrongHG(hg_version, valid)   
 
def scrape_seq(primer_name, f_primer, r_primer, hg_version, 
               max_size, min_perfect, min_good):
    ''' Use primer pairs to scrape in-silico PCR amplicon sequences
        from UCSC.
    '''
    url = 'https://genome.ucsc.edu/cgi-bin/hgPcr?hgsid=483' \
          '751629_vuLjoO4UVF9h4vF4TEp9U8OQiFd7&' \
          'org=Human&db={}&wp_target=genome&wp_f={}' \
          '&wp_r={}&Submit=submit&wp_size={}&wp_' \
          'perfect={}&wp_good={}&boolshad.wp_flipReverse=0'
    link = url.format(hg_version, f_primer, r_primer, max_size, min_perfect, min_good)
    req = requests.get(link)              
    req.raise_for_status()              

    entire_url = bs4.BeautifulSoup(req.text,"html.parser")
    pre_elements = entire_url.select('pre') 
    if not pre_elements:
        raise ex.NoAmplicon(primer_name)

    html_to_text = pre_elements[0].getText()         

    product_number = html_to_text.count(">")
    if product_number > 1:
        raise ex.MultipleAmplicons(primer_name, product_number)
  
    return html_to_text 
               
def seperate_data(text):
    ''' From scrapped in-silico PCR amplicon information, 
        seperate the amplicon sequence and its metadata.
    '''
    split_text = list(filter(None, text.split("\n")))
    seq = ''.join(split_text[1:])
    header = split_text[0]
    return (header, seq)
   
def get_metadata(header, seq, hg_version):
    ''' Gather metadata from the isPCR results and MetaData.
    '''
    pos_range = header.split(" ")[0]
    pos_range = pos_range[1:].replace("+", "-")
    gc = (seq.upper().count('C') + seq.upper().count('G')) / len(seq)
    size = "{}bp".format(len(seq))
    # get gene metadata from the middle of the amplicon
    chrom, start, end = re.split(':|-', pos_range)
    pos = int(start) - int(size.replace('bp', ''))
    data = metadata.LocusMetaData(chrom, pos, hg_version)
    if data.exon.exon:
        exon = data.exon.number
        intron = '-'
    else:
        exon = '-'
        intron = data.exon.number
    gene_metadata = (data.gene.name, exon, intron, size, 
                     pos_range, round(gc, 3)*100)
    return gene_metadata

def parse2output(args, header):
    ''' Write the results of parsing the input file
        contents through unknown_primer to an output 
        file.
    ''' 
    with open(args['output'], 'w') as out:
        out.write(header+"\n")
        with open(args['input'], 'r') as in_file:
            for line in in_file:
                try:
                    primer_name, f_primer, r_primer, hg_version = line.rstrip("\n").split("\t")
                    logging.info('Processing primer {}....'.format(primer_name))
                    metadata = unknown_primer(f_primer, r_primer, hg_version, 
                                              primer_name, args['max_size'], 
                                              args['min_perfect'], args['min_good'])
                    format_metadata = '\t'.join([str(x) for x in metadata])
                    out.write(format_metadata+"\n")

                except (ex.MultipleAmplicons, ex.NoAmplicon, ex.WrongHG, 
                        ex.AmbigousBase) as e:
                    logging.error(e.msg)
                    continue

                except ex.NoProteinCodingTranscript as e:
                    logging.error('No protein coding transcript was found within the in-silico amplicon generated by {}'.format(primer_name))

def print_metadata(args, header):
    ''' Print the results of parsing a primer pair
        through unknonw_primer.
    '''
    metadata = unknown_primer(args['primers'][0], args['primers'][1], 
                              args['genome_version'], 'query', args['max_size'],
                              args['min_perfect'], args['min_good'])
    print_metadata = '\t'.join([str(x) for x in metadata])
    print(header+"\n"+print_metadata)

def get_parser():
    parser = argparse.ArgumentParser(description='retrieve a primer pairs metadata from UCSC and Ensembl')
    parser.add_argument('-p', '--primers', nargs=2, type=str, help='primer pair')
    parser.add_argument('-i', '--input', type=str, help='a tab delimited file containing primer info')
    parser.add_argument('-hg', '--genome_version', type=str, help='human genome version (default=hg19)', default='hg19')
    parser.add_argument('-m', '--max_size', type=int, help='maximum in-silico amplicon size (default=4000)', default=4000)
    parser.add_argument('-mp', '--min_perfect', type=int, help='no. of bases that match exactly on 3end of primers (default=15)', default=15)
    parser.add_argument('-mg', '--min_good', type=int, help='no. of bases on 3end of primers where at least 2/3 bases match (default=15)', default=15)
    parser.add_argument('-o', '--output', type=str, help='output data to file name', default='output_primer_data.txt')
    return parser

def cli():
    parser = get_parser()
    args = vars(parser.parse_args())
    header = '\t'.join(('Primer', 'F_Primer', 'R_Primer', 'Genome', 'Gene', 'Exon',
                        'Intron', 'Product_Size', 'Primer_Range', 'GC%'))
    if args['input']:
        parse2output(args, header)
    else:
        print_metadata(args, header)  


if __name__ == '__main__':
    cli()
