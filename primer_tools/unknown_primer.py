from __future__ import division
import requests,re, os, csv, bs4, urllib2, click

class AmbiguousBaseError(Exception):
    pass
    
class NoAmplicon(Exception):
    pass

class MultipleAmplicons(Exception):
    pass

class IncorrectVariant(Exception):
    pass



def unknown_primer(input_file=None,
                   output_file="unknown_primer_output.txt",
                   DB="hg19",):
                   
    ''' For each primer pair output the: genomic region, amplicon size, number
        of amplicons generated and GC% of the amplicon along with the primer
        pair name to an output file
        
        
        DB         -     refers to the human genome version i.e. hg19, hg38
        input_file -     contains primer information for each pair on a newline,
                         which should conatin primer name, forward primer sequence,
                         reverse primer sequence (tab deliminated). 
                         
                         
    '''
    
    if os.path.isfile(input_file) is True:
        
        output = open(output_file,"w")
        header = ("Primer","F-Primer","R-Primer","Product_Size",
              "Genomic_Region","GC%","Number_PCR_Products"+"\n")
        output.write("\t".join(header))
        
        for primer in open(input_file,"r"):
            
            primer = primer.rstrip("\n").split("\t")
            primer_name = primer[0]
            f_primer = primer[1].upper()
            r_primer = primer[2].upper()
    
            try:
                answer = get_unknown_primer_info(primer_name,f_primer,r_primer,DB)
                print answer
                output.write(answer+"\n")
            except AmbiguousBaseError:
                print "Skipping invalid base in primer:"+primer[0]
            except NoAmplicon:
                print "No amplicon generated from isPCR for primer: "+primer[0]
            except MultipleAmplicons:
                print "The following primers generate more than one amplicon:"+primer[0]
        output.close()
        
        
    if os.path.isfile(input_file) is False:
        f_primer = input_file.split(",")[0]
        r_primer = input_file.split(",")[1]
        answer = get_unknown_primer_info("query",f_primer,r_primer,DB)
        print answer
    
def get_unknown_primer_info(primer_name="query",f_primer=None,r_primer=None,DB="hg19"):
	''' The DB, forward primer and reverse primer sequences are used as part 
	    of the UCSC isPCR link for webscraping, which are further filtered for 
	    pre elements containing the in silico amplicon sequence.
	    
	    This sequence is then further manipulated to attain the genomic position, 
	    amplicon size, number of amplicons generated and GC% of the amplicon.
	'''	
	
        
        if re.search(r'[^ATCG]',f_primer)or re.search(r'[^ATCG]',r_primer):
            raise AmbiguousBaseError("Primers must contain ATGC bases only")
            
        
        req = requests.get("http://genome.ucsc.edu/cgi-bin/hgPcr?db="+DB+\
        "&wp_target=genome&wp_f="+f_primer+"&wp_r="+r_primer+\
        "&wp_size=4000&wp_perfect=15&wp_good=15&boolshad.wp_flipReverse=0")
        
        req.raise_for_status()   # get the request error code if failed
        
        entire_url = bs4.BeautifulSoup(req.text,"html.parser")
        
        # if no pre_elemnts then no amplicon generated from the seqence
        pre_elements = entire_url.select('pre') # get all <pre> elements on webpage
        
        if not pre_elements:
            raise NoAmplicon("No amplicon generated")
        
        isPCR = pre_elements[0].getText()  # get text for first <pre> element
        
        
        amplicon = "\n".join(isPCR.split("\n")[1:]) # split newlines, take 2nd to last and join back together
        
        amplicon_header = "\n".join(isPCR.split("\n")[:1])
        split_header =  amplicon_header[1:].split(" ")
        
        region = split_header[0].replace("+","-").replace("chr","")
        amplicon_size = split_header[1]
        product_number = len(pre_elements)
        
        if product_number > 1:
            raise MultipleAmplicons
        
        # should be seperate module
        gc_percent = round(((amplicon.count("G")+amplicon.count("C")+amplicon.count("c")+\
                    amplicon.count("g"))/ len(amplicon)) * 100,2)
        
        output = (primer_name,f_primer,r_primer,str(amplicon_size),
                  region,str(gc_percent)+"%",str(product_number))
        return "\t".join(output)

unknown_primer("TAACAGATTGATGATGCATGAAATGGG,CCCATGAGTGGCTCCTAAAGCAGCTGC")
#unknown_primer("T.txt")

