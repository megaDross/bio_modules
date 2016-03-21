from __future__ import division
import requests,re, os, bs4, click


# HG VERSION ALWAYS HG19 DESPITE WHAT YOU USE AS THE ARGUMENT FOR --HG_VERSION

class AmbiguousBaseError(Exception):
    pass
    
class NoAmplicon(Exception):
    pass

class MultipleAmplicons(Exception):
    pass

class IncorrectVariant(Exception):
    pass

@click.command()
@click.option('--primers',nargs=2)
@click.option('--input_file',nargs=1) # -1 means unlimited number of arguments
@click.option('--output_file',default=None,help="default: unknwon_primer_output.txt")
@click.option('--hg_version',default="hg19",help="human genome version. default: hg19")
@click.option('--delimiters',default="\t",help="character used to split F- & R-Primer, STRING only")

def unknown_primer(primers=None, input_file=None,output_file=None,
                   hg_version="hg19",delimiters=","):
                   
    ''' Takes primer sequences as input and uses the UCSC isPCR tool
        to generate amplicon sequence information.
        \b\n
        A file or string can be used as input.
        STRING: --primers option with 2 arguments (primer seqs) given.
        FILE: --input_file option deliminated file with the primer name, F-primer seq anad R-primer seq.
        \b\n                    
        Example:\b\n
        python unknown_primer.py --input_file in.txt --output_file out.txt --hg_version hg38
        python unknown_primer.py --primers CGATCGTTGC GCTAGTCGT 
    '''
    
    # if input is two arguments   
    if primers:
        input_file = ["query"+delimiters+primers[0]+delimiters+primers[1]]
        return process_primer_info(input_file,output_file,hg_version,delimiters)
    
    # if input is a file
    if os.path.isfile(input_file) is True:
        input_file = open(input_file,"r+")
        primer_info = process_primer_info(input_file,output_file,hg_version,delimiters)
        
        if output_file is not None:
            output = open(output_file,"w")
            for primers in primer_info:
                output.write(primers)
            output.close()
            return output
        
        return primer_info
        

def process_primer_info(input_file,output_file,hg_version,delimiters):
    primer_info = []
    for primer in input_file:
        try:
            primer = primer.rstrip("\n").split(delimiters)
            primer_name = primer[0]
            f_primer = primer[1].upper()
            r_primer = primer[2].upper()
            if re.search(r'[^ATCG]',f_primer) or re.search(r'[^ATCG]',r_primer):
                raise AmbiguousBaseError
            amplicon_info = get_unknown_primer_info(primer_name,f_primer,r_primer,hg_version)
    
            print amplicon_info
            primer_info.append(amplicon_info)

        except AmbiguousBaseError:
            print "Skipping, non-ATGC base found in: "+primer_name
        except NoAmplicon:
            print "No amplicon generated from isPCR for primer: "+primer_name
        except MultipleAmplicons:
            print "The following primers generate more than one amplicon:"+primer_name
        
    
def get_unknown_primer_info(primer_name="query",f_primer=None,r_primer=None,hg_version=None):
	''' Generate an amplicon sequence from inputted primer sequences, which
	    is further manipulated t derive inofrmation from the sequence.
	'''	
	
        # ensure input is only bases
        if re.search(r'[^ATCG]',f_primer)or re.search(r'[^ATCG]',r_primer):
            raise AmbiguousBaseError("Primers must contain ATGC bases only")
        
            
        # generate amplicon sequence
        req = requests.get("http://genome.ucsc.edu/cgi-bin/hgPcr?hg_version="+hg_version+\
                           "&wp_target=genome&wp_f="+f_primer+"&wp_r="+r_primer+\
                           "&wp_size=4000&wp_perfect=15&wp_good=15&boolshad"+\
                           ".wp_flipReverse=0")
        req.raise_for_status()                                # get the request error code if failed
        entire_url = bs4.BeautifulSoup(req.text,"html.parser")
        pre_elements = entire_url.select('pre')               # get all <pre> elements on webpage
        if not pre_elements:
            raise NoAmplicon("No amplicon generated")
        isPCR = pre_elements[0].getText()                     # get text for first <pre> element
        amplicon = "\n".join(isPCR.split("\n")[1:])           # split newlines, take 2nd to last and join back together
        
        
        # use amplicon sequence and header to get additional information
        amplicon_header = "\n".join(isPCR.split("\n")[:1])
        split_header =  amplicon_header[1:].split(" ")
        region = split_header[0].replace("+","-").replace("chr","")
        amplicon_size = split_header[1]
        gc_percent = round(((amplicon.count("G")+amplicon.count("C")+amplicon.count("c")+\
                    amplicon.count("g"))/ len(amplicon)) * 100,2)
        product_number = len(pre_elements)
        if product_number > 1:
            raise MultipleAmplicons        
        
        
        # return scraped information
        output = (primer_name,f_primer,r_primer,str(amplicon_size),
                  region,str(gc_percent)+"%",str(product_number))
        return "\t".join(output)


if __name__ == '__main__':
    unknown_primer()
