from __future__ import division
import requests,re, csv, bs4

class AmbiguousBaseError(Exception):
    pass
    
class NoAmplicon(Exception):
    pass

class MultipleAmplicons(Exception):
    pass

# unknown_primer wont write to the file, file empty after processing




def unknown_primer(DB,input_file,output_file):
    ''' For each primer pair output the: genomic region, amplicon size, number
        of amplicons generated and GC% of the amplicon along with the primer
        pair name to an output file
        
        
        DB         -     refers to the human genome version i.e. hg19, hg38
        input_file -     contains primer information for each pair on a newline,
                         which should conatin primer name, forward primer sequence,
                         reverse primer sequence (tab deliminated). 
                         
                         
    '''
    output = open(output_file,"w")
    header = ("Primer","Genomic_Position","Product_Size","Number_PCR_Products","GC%"+"\n")
    output.write("\t".join(header))
    
    for primer in open(input_file,"r"):
        primer = primer.rstrip("\n").split("\t")

        try:
            print get_unknown_primer_info(DB,primer)
            output.write(get_unknown_primer_info(DB,primer)+"\n")
        except AmbiguousBaseError:
            print "Skipping invalid base in primer:"+primer[0]
        except NoAmplicon:
            print "No amplicon generated from isPCR for primer: "+primer[0]
        except MultipleAmplicons:
            print "The following primers generate more than one amplicon:"+primer[0]
    output.close()
    
    
def get_unknown_primer_info(DB, input_file):
	''' The DB, forward primer and reverse primer sequences are used as part 
	    of the UCSC isPCR link for webscraping, which are further filtered for 
	    pre elements containing the in silico amplicon sequence.
	    
	    This sequence is then further manipulated to attain the genomic position, 
	    amplicon size, number of amplicons generated and GC% of the amplicon.
	'''	
	     
	f_primer = input_file[1].upper()
        r_primer = input_file[2].upper()
        
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
        
        gc_percent = round(((amplicon.count("G")+amplicon.count("C")+amplicon.count("c")+\
                    amplicon.count("g"))/ len(amplicon)) * 100,2)
        
        output = (input_file[0],region,str(amplicon_size),str(product_number),str(gc_percent)+"%")
        return "\t".join(output)

print unknown_primer("hg19","TAAD_Primer_Validation_Database.txt","moomoo.txt")


def region_extractor(input_file, output_file, delimiters=None, number_upstream=None, number_downstream=None):
        ''' Produces a sequence 20bp upstream and downstream from
        the defined variant position in the input_file

        The input_file should be in the following format:
                sample_name\tchromosome_number:nucleotide_number
        i.e.    sample17\t2:189851842

        The delimiters argument is defaulted to tab
        number_upstream & number_downstream defaulted to 20
        '''
        if delimiters is None:
                delimiters = "\t"
                
        if number_upstream is None:
                number_upstream = 20
        
        if number_downstream is None:
                number_downstream = 20
        
        alt_T=[]
        # change dir manually in python
        T = open(input_file,"r+")
        output=open(output_file,"w")


        reader = csv.reader(T, delimiter=delimiters)
        for changes in reader:
            nospace = changes[1].replace(" ","")
            chrom = nospace.split(":")[0]
            pos = nospace.split(":")[1]
            start_pos = int(pos) - number_downstream
            end_pos = int(pos) + number_upstream
            alt_T.append(changes[0]+"\t"+chrom+":"+str(start_pos)+","+str(end_pos))

        for alted_T in alt_T:
                region = alted_T.split("\t")[1]
                test = urllib2.urlopen("http://genome.ucsc.edu/cgi-bin/das/hg19/dna?segment="+region)
                search = test.read()
                seq = re.findall(r"[tacg{40}].*",search)[4]
                downstream = seq[0:number_downstream]
                var = seq[number_downstream]
                upstream = seq[number_upstream+1:len(seq)]
                answer = downstream+"-"+var+"-"+upstream
                print alted_T.split("\t")[0]+"\t"+region+"\t"+answer+"\n"
                output.write(alted_T.split("\t")[0]+"\t"+region+"\t"+answer+"\n")
                        
        output.close()


		

region_extractor("all_vars_in.txt","GIVE_ME_ANSWER.txt")