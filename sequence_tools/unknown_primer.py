from __future__ import division
import requests,re, os, bs4, click
from useful_tools import useful


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
@click.option('--delimiters',default="\t",help="--input_file option only. default: tab")

def unknown_primer(primers=None, input_file=None,output_file=None,
                   hg_version="hg19",delimiters=","):
                   
    ''' Takes primer sequences as input and uses the UCSC isPCR tool
        to generate amplicon sequence information.
        \b\n
        A file or string can be used as input.
        STRING: --primers option with 2 arguments (primer seqs) given.
        FILE: --input_file option deliminated file with the primer name, 
        F-primer seq anad R-primer seq.
        \b\n                    
        Example:\b\n
           unknown_primer --input_file in.txt --output_file out.txt --hg_version hg38\n
           unknown_primer --primers TAACAGATTGATGATGCATG CCCATGAGTGGCTCCTAAA 
    '''

    # if input is two arguments   
    if primers:
        input_file = ["query"+delimiters+primers[0]+delimiters+primers[1]]
        return process_primer_info(input_file,output_file,hg_version,delimiters)
    
    # if input is a file
    if os.path.isfile(input_file) is True:
        input_file = open(input_file,"r+")
        primer_info = process_primer_info(input_file,output_file,hg_version,delimiters)
       
        # if the --output_file is used, write the primer_info to a file
        if output_file is not None:
            output = open(output_file,"w")
            header = "\t".join(("primer_name","F-primer","R-primer","amplicon_size",
                               "genomic_range","gc%","number_amplicons","\n"))
            output.write(header)
            
            for primers in primer_info:
                output.write(primers+"\n")
            output.close()
            return output
        
        return primer_info
                

def process_primer_info(input_file,output_file,hg_version,delimiters):
    '''
    iterate through each line of the input and parse them into get_unknown_primer_info
    '''
    # adds all scrapped data to a list, which is written to an output, only --ouput_file is selected
    primer_info = []
    
    for primer in input_file:
        try:
            # splits up data into primer name, forward primer and reverse primer
            primer = primer.rstrip("\n\r").split(delimiters)
            primer_name = primer[0]
            f_primer = primer[1].upper()
            r_primer = primer[2].upper()
            
            # ensure input is bases only
            if re.search(r'[^ATCG]',f_primer) or re.search(r'[^ATCG]',r_primer):
                raise AmbiguousBaseError
            
                    
            amplicon_info = get_unknown_primer_info(hg_version,primer_name,f_primer,r_primer)
    
            print(amplicon_info)
            primer_info.append(amplicon_info)
        
        except IndexError:
            print(str(primer)+" is improperly deliminated")
        except AmbiguousBaseError:
            print("Skipping, non-ATGC base found in: "+primer_name)
        except NoAmplicon:
            print("No amplicon generated from isPCR for primer: "+primer_name)
        except MultipleAmplicons:
            print("The following primers generate more than one amplicon:"+primer_name)
        
    return primer_info
    
def get_unknown_primer_info(hg_version, primer_name="query",f_primer=None,r_primer=None):
    ''' Generate an amplicon sequence from inputted primer sequences, which
        is further manipulated t derive inofrmation from the sequence.
    '''         
        # generate amplicon sequence using isPCR tool
    req = requests.get("https://genome.ucsc.edu/cgi-bin/hgPcr?hgsid=483751629_vuLjoO4UVF9h4vF4TEp9U8OQiFd7&org=Human&db="+hg_version+"&wp_target=genome&wp_f="+f_primer+"&wp_r="+r_primer+"&Submit=submit&wp_size=4000&wp_perfect=15&wp_good=15&boolshad.wp_flipReverse=0")
    req.raise_for_status()                                # get the request error code if failed
    entire_url = bs4.BeautifulSoup(req.text,"html.parser")
    pre_elements = entire_url.select('pre')               # get all <pre> elements on webpage
        
    # if nothing between pre-elements, raise error
    if not pre_elements:
        raise NoAmplicon("No amplicon generated")
    isPCR = pre_elements[0].getText()                     # get text for first <pre> element
    amplicon = "\n".join(isPCR.split("\n")[1:])           # split newlines, take 2nd to last and join back together


    # use amplicon sequence and header to get additional information
    amplicon_header = "\n".join(isPCR.split("\n")[:1])
    split_header =  amplicon_header[1:].split(" ")
    region = split_header[0].replace("+","-").replace("chr","")
    amplicon_size = split_header[1]
    gc_percent = useful.gc_content(amplicon)
    product_number = len(pre_elements)
    if product_number > 1:
        raise MultipleAmplicons        


    # return scraped information
    output = (primer_name,f_primer,r_primer,str(amplicon_size),
              region,str(gc_percent)+"%",str(product_number))
    return "\t".join(output)


if __name__ == '__main__':
    unknown_primer()
