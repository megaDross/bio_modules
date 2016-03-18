from __future__ import division
import requests,re, os, bs4, click

class AmbiguousBaseError(Exception):
    pass
    
class NoAmplicon(Exception):
    pass

class MultipleAmplicons(Exception):
    pass

class IncorrectVariant(Exception):
    pass

@click.command()

@click.argument('input_file',nargs=1)
@click.option('--output_file',default="unknown_primer_output.txt",help="default: unknwon_primer_output.txt")
@click.option('--hg_version',default="hg19",help="human geome version. default: hg19")
@click.option('--delimiters',default=",",help="character used to split F- & R-Primer, STRING only")

def unknown_primer(input_file=None,output_file="unknown_primer_output.txt",
                   hg_version="hg19",delimiters=","):
                   
    ''' Takes primer sequences as input and uses the UCSC isPCR tool
        to generate amplicon sequence information.
        \b\n
        A file or string can be used as input.
        STRING: a forward and reverse primer deliminated by a comma.
        FILE: deliminated file withe the primer name, F-primer seq anad R-primer seq
        \b\n                    
        Example:\b\n
        python unknown_primer.py input.txt --output_file out.txt --hg_version hg38
        python unknown_primer.py CGATCGTTGC-GCTAGTCGT --delimiters -
    '''
    
               
    # if input is a file
    if os.path.isfile(input_file) is True:
        with open(output_file,"w") as output:
            header = ("Primer","F-Primer","R-Primer","Product_Size",
                        "Genomic_Region","GC%","Number_PCR_Products"+"\n")
            output.write("\t".join(header))
            
            for primer in open(input_file,"r"):
                try:
                    primer = primer.rstrip("\n").split("\t")
                    primer_name = primer[0]
                    f_primer = primer[1].upper()
                    r_primer = primer[2].upper()
                    amplicon_info = get_unknown_primer_info(primer_name,f_primer,r_primer,hg_version)
            
                    print amplicon_info
                    output.write(amplicon_info+"\n")
                except IndexError as e:
                    print e.args[0]+": Comma is defaulted to be used as a deliminator for the forward and reverse  primers. To override this use the --delimiters option"
                except AmbiguousBaseError:
                    print "Skipping, non-ATGC base found in: "+primer_name+"\n\ncomma is defaulted as a deliminate F&R primer. To override use the --delimiters option."
                except NoAmplicon:
                    print "No amplicon generated from isPCR for primer: "+primer_name
                except MultipleAmplicons:
                    print "The following primers generate more than one amplicon:"+primer_name
            
    # if input is a string    
    if os.path.isfile(input_file) is False:
        try:
            primer_name = "query"
            f_primer = input_file.split(delimiters)[0].upper()
            r_primer = input_file.split(delimiters)[1].upper()
            if re.search(r'[^ATCG,]',f_primer) or re.search(r'[^ATCG ]',r_primer):
                raise AmbiguousBaseError
            amplicon_info = get_unknown_primer_info(primer_name,f_primer,r_primer,hg_version)
            print amplicon_info
        
        except IndexError as e:
            print e.args[0]+": Comma is defaulted to be used as a deliminator for the forward and reverse  primers. To override this use the --delimiters option"
        except AmbiguousBaseError:
            print "Skipping, non-ATGC base found in: "+primer_name+"\n\ncomma is defaulted as a deliminate F&R primer. To override use the --delimiters option."
        except NoAmplicon:
            print "No amplicon generated from isPCR for primer: "+primer_name
        except MultipleAmplicons:
            print "The following primers generate more than one amplicon:"+primer_name


            
    
def get_unknown_primer_info(primer_name="query",f_primer=None,r_primer=None,hg_version="hg19"):
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
