"""
utilise various UCSC tools via the command line
"""
from __future__ import division
import urllib2, re, csv

# unknown_primer wont write to the file, file empty after processing

def unknown_primer(DB, input_file, output_file, delimiters="\t"):
	'''Generates coordinates for which a primer pair bind to within the genome, 
	and gives the PCR products size and number of potential PCR products
		
	Input file should be in the following format:
		primer_name\tF-primer\tR_primer

	The delimiters argument is defaulted to tab, it is assumed a header is present.
	'''
	output_file = open(output_file,"w")
	reader = csv.reader(open(input_file,"r+"), delimiter=delimiters)
	skip_header = reader.next()
	
	for primer in reader:
		ispcr = urllib2.urlopen("http://genome.ucsc.edu/cgi-bin/hgPcr?db="+DB+\
		"&wp_target=genome&wp_f="+primer[1]+"&wp_r="+primer[2]+\
		"&wp_size=4000&wp_perfect=15&wp_good=15&boolshad.wp_flipReverse=0")
		search = ispcr.read()

        	if re.search(r'No matches.*',search):
                	output = primer[0]+"\t"+primer[1]+"\t"+primer[2]+"\t"+"none"+"\t"+"no match"+"\n"
                	return output
                	output_file.write(output)

        	elif re.search(r'No matches.*',search) is None:
                	
                	# generate genomic region co-ordinates and PCR product size
                	shortened = re.search(r'chr.*[0-9]&',search).group()
                	region = shortened[:len(shortened)-1]
                	size = str((int(re.split(r'[-]',region)[1]) - int(re.split(r'[:-]',region)[1]))+1)
                	
                	# number of PCR products
			get_product_number = re.findall(r'>chr.*[0-9].*',search)
			product_number = len(get_product_number)
			
			# prep for getting amplicon info
			product_seq = re.findall(r"[ACTGatcg]{25,54}.*",search, re.DOTALL)
                        no_special_chars = re.sub(r'([\[\]\'\, ,\\n])',"",str(product_seq))    
                        
                        # sequence of amplicon/PCR-product
                        amplicon = re.search(r"[ACTG]{10,40}[actg]{25,1000}[ACTG]{10,40}",\
                        no_special_chars, re.DOTALL).group()
                        
                        # GC% of the amplicon
                        GC_percent = ((amplicon.count("G")+amplicon.count("C")+amplicon.count("c")+\
                        amplicon.count("g"))/ len(amplicon)) * 100
			
			output = primer[0]+"\t"+primer[1]+"\t"+primer[2]+"\t"+size+"\t"+region\
			+"\t"+str(product_number)+"\t"+str(GC_percent)+"\n"
                	return output
                	output_file.write(output)
        output_file.close()


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


def isPCR(DB, input_file, output_file, delimiters=None):
	''' Use the UCSC isPCR tool to get primer information
	it is basically EXACTLY the same as unknown_primer
	consider deleting 
	'''
	if delimiters is None:
                delimiters = "\t"

        test = open(output_file,"w")
        reader = csv.reader(open(input_file,"r+"), delimiter=delimiters)

        for primer in reader:
                ispcr = urllib2.urlopen("http://genome.ucsc.edu/cgi-bin/hgPcr?db="+DB+\
                "&wp_target=genome&wp_f="+primer[1]+"&wp_r="+primer[2]+\
                "&wp_size=4000&wp_perfect=15&wp_good=15&boolshad.wp_flipReverse=0")
                search = ispcr.read()
		
		all = re.findall(r'>chr.*[0-9].*',search)
		ans = str(all).replace("</A>","").replace(">","").replace("[","").replace("]","").replace("'","")
		answer = str(ans).split(",")
		print primer[0] + "\t" + ans + "\t"+str(len(all))
		
