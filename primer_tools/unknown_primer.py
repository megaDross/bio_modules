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

@click.option('--input_file')
@click.option('--output_file',default="unknown_primer_output.txt")
@click.option('--hg_version',default="hg19")

def unknown_primer(input_file=None,
                   output_file="unknown_primer_output.txt",
                   hg_version="hg19",):
                   
    ''' Takes primer sequences as input and uses the UCSC isPCR tool
        to generate amplicon sequence information.
        
        Arguments
        ------------------------------------------------------------------------
        
        input_file -- either a string where the forward and reverse primer sequence
                      deliminated by a comma or a tab deliminated file with the 
                      following information:
                         primer name   forward primer sequence   reverse primer sequence
        
        output_file -- created only if the inut given is a file opposed to a string.
                       If no name is specified for output_file, then it is automatically
                       named unknown_primer_output.txt
        
        hg_version -- refers to the human genome version i.e. hg19, hg38. 
                                defaulted to hg19 
        
        
        Returns
        ------------------------------------------------------------------------
        
        amplicon_info -- amplicon information given in a tab deliminated format:
                           primer_name, f_primer_seq, r_primer_seq,amplicon_size, 
                           genomic_region, amplicon_GC%, number_products
                           
        Examples
        ------------------------------------------------------------------------
        from primer_tools import unknown_primer
        
        unknown_primer.unknown_primer("TAACAGATTGATGATGCATGAAATGGG,CCCATGAGTGGCTCCTAAAGCAGCTGC")
        
        unkown_primer.unknown_primer("primer_file.txt","primer_information.txt","hg19")

        OR
        
        using click:
            
        python unknown_primer --input_file primer_file.txt --hg_version hg38
                          
    '''
    try:
                  
        # if input is a file
        if os.path.isfile(input_file) is True:
            with open(output_file,"w") as output:
                header = ("Primer","F-Primer","R-Primer","Product_Size",
                          "Genomic_Region","GC%","Number_PCR_Products"+"\n")
                output.write("\t".join(header))
                
                for primer in open(input_file,"r"):
                    primer = primer.rstrip("\n").split("\t")
                    primer_name = primer[0]
                    f_primer = primer[1].upper()
                    r_primer = primer[2].upper()
                    amplicon_info = get_unknown_primer_info(primer_name,f_primer,r_primer,hg_version)
            
                    print amplicon_info
                    output.write(amplicon_info+"\n")
               
        # if input is a string    
        if os.path.isfile(input_file) is False:
            primer_name = "query"
            f_primer = input_file.split(",")[0].upper()
            r_primer = input_file.split(",")[1].upper()
            if re.search(r'[^ATCG,]',f_primer) or re.search(r'[^ATCG ]',r_primer):
                raise AmbiguousBaseError
            amplicon_info = get_unknown_primer_info(primer_name,f_primer,r_primer,hg_version)
            print amplicon_info
            
    except IndexError as e:
        print e.args[0]+": Only a comma can be used as a deliminator for the forward and reverse  primers"
    except AmbiguousBaseError:
        print "Skipping, non-ATGC base found in: "+primer_name+"\n\nonly a comma can be used to deliminate F&R primer"
    except NoAmplicon:
        print "No amplicon generated from isPCR for primer: "+primer_name
    except MultipleAmplicons:
        print "The following primers generate more than one amplicon:"+primer_name
    
    
    
def get_unknown_primer_info(primer_name="query",f_primer=None,r_primer=None,hg_version="hg19"):
	''' Generate an amplicon sequence from inputted primer sequences, which
	    is further manipulated t derive inofrmation from the sequence
	    
	    Arguments
	    --------------------------------------------------------------------
	    primer_name -- defaulted to query
	    
	    f_primer -- forward primer sequence
	    
	    r-Primer -- reverse primer sequence
	    
            hg_version -- refers to the human genome version i.e. hg19, hg38. 
                                defaulted to hg19 
            
            Returns
            ------------------------------------------------------------------
            output -- amplicon information derived from the web-scraping
            
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
