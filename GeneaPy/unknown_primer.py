from __future__ import division
import requests,re, os, bs4, click
import useful
import get_gene_exon_info
from output import write_to_output

class AmbiguousBaseError(Exception):
    pass
    
class NoAmplicon(Exception):
    pass

class MultipleAmplicons(Exception):
    pass

class IncorrectVariant(Exception):
    pass

    
@click.command('unknown_primer')
@click.option('--primers',nargs=2,help="accepts 2 arguments only; forward primer reverse primer")
@click.option('--input_file',nargs=1, help="file used as input")
@click.option('--output_file',default=None,help="only available if the --input_file option is used")
@click.option('--hg_version',default="hg19",help="human genome version. default: hg19")

def unknown_primer(primers=None, input_file=None,output_file=None,
                   hg_version="hg19"):
                   
    ''' Takes primer sequences as input and uses the UCSC isPCR tool
        to generate amplicon sequence information.
        \b\n
        A file or string can be used as input.
        STRING: --primers option with 2 arguments (primer seqs) given.
        FILE: --input_file option deliminated file with the primer name, 
        F-primer seq and R-primer seq.
        \b\n                    
        Example:\b\n
           unknown_primer --input_file in.txt --output_file out.txt --hg_version hg38\n
           unknown_primer --primers TAACAGATTGATGATGCATG CCCATGAGTGGCTCCTAAA 
    '''
    header = "\t".join(("Primer", "F_Primer","R_Primer", "Gene", "Product_Size",
                        "Primer_Range","GC%","Number_Amplicons", "\n"))
    print(header)

    # determine whether the input is a file or string and process accordingly 
    if input_file:
        all_matched_primers = []
        for line in [line.rstrip("\n").split("\t") for line in open(input_file)]:
            primer_name = line[0]
            f_primer = line[1].upper()
            r_primer = line[2].upper()
            full_amplicon_info = get_all_primer_info(primer_name, hg_version, 
                                                     f_primer, r_primer)
            if not full_amplicon_info:
                full_amplicon_info = "\t".join((primer_name, "-", "-", "-", 
                                                "-", "-", "-", "-"))
            all_matched_primers.append(full_amplicon_info)
            print(full_amplicon_info)
    else:
        primer_name = "query"
        f_primer = primers[0]
        r_primer = primers[1]
        full_amplicon_info = get_all_primer_info(primer_name, hg_version, 
                                                 f_primer, r_primer)
        print(full_amplicon_info)

    if output_file:
        write_to_output(all_matched_primers, output_file, header)




def get_all_primer_info(primer_name, hg_version, f_primer, r_primer):
    ''' Scrape primer information from UCSC and gene information from
        Ensembl

        Returns a string
    '''
    amplicon_info = get_unknown_primer_info(primer_name, hg_version, 
                                            f_primer,r_primer)
    if amplicon_info:
        primer_range = amplicon_info[3]
        gene_name = get_gene_name(primer_range, hg_version)
        reorder_info = (primer_name,) + amplicon_info[:2] + (gene_name,) + amplicon_info[2:]
        full_amplicon_info = "\t".join(reorder_info)  
        return full_amplicon_info
    
    else:
        pass



def get_unknown_primer_info(primer_name, hg_version,f_primer=None,r_primer=None):
    ''' Generate an amplicon sequence from inputted primer sequences, which
        is further manipulated to gain the proposed amplicons metadata.
    
        Returns a tuple
    '''     
    try:
        # check for ambigous bases
        if re.search(r'[^ATCG]',f_primer) or re.search(r'[^ATCG]',r_primer):
            raise AmbiguousBaseError

        # generate amplicon sequence using isPCR tool
        req = requests.get('https://genome.ucsc.edu/cgi-bin/hgPcr?hgsid=483'
                           '751629_vuLjoO4UVF9h4vF4TEp9U8OQiFd7&org=Human&db='
                           hg_version+'&wp_target=genome&wp_f="+f_primer+"&wp_r='+r_primer+
                           '&Submit=submit&wp_size=4000&wp_perfect=15&wp_good=15&boolshad.wp_flipReverse=0')
        req.raise_for_status()                                # get the request error code if failed
        entire_url = bs4.BeautifulSoup(req.text,"html.parser")
        pre_elements = entire_url.select('pre')               # get all <pre> elements on webpage
        # if nothing between pre-elements, raise error
        if not pre_elements:
            raise NoAmplicon

        isPCR = pre_elements[0].getText()                     # get text for first <pre> element
        amplicon = "\n".join(isPCR.split("\n")[1:])           # split newlines, take 2nd to last and join back together

        # use amplicon sequence and header to get additional information
        amplicon_header = "\n".join(isPCR.split("\n")[:1])
        split_header = amplicon_header[1:].split(" ")
        region = split_header[0].replace("+","-").replace("chr","")
        amplicon_size = split_header[1]
        gc_percent = useful.gc_content(amplicon)
        product_number = len(pre_elements)
        if product_number > 1:
            raise MultipleAmplicons        


        # return scraped information
        output = (f_primer,r_primer,str(amplicon_size),
                  region,str(gc_percent)+"%",str(product_number))
        return output
   
    except AmbiguousBaseError:
        print("Skipping, non-ATGC base found in: "+primer_name)
    except NoAmplicon:
        print("No amplicon generated from isPCR for primer: "+primer_name)
    except MultipleAmplicons:
        print("The following primers generate more than one amplicon:"+primer_name)



def get_gene_name(primer_range, hg_version):
    ''' Retrieve gene, transcript and exon information from the genomic
        region in which the proposed amplicon derive from

        Returns a tuple of tuples
    '''
    all_info = get_gene_exon_info.gene_transcript_exon(primer_range, hg_version)
    gene_name, gene_id, gene_type, gene_range = all_info[0]
    transcript = all_info[1]
    exon_id, intron, exon = all_info[2]

    return gene_name



if __name__ == '__main__':
    unknown_primer()
