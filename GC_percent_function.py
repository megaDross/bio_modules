from __future__ import division

import urllib2
import re
import csv



def GC_content(DB, input_file, delimiters=None):
	''' Give GC% of an amplicon produced from inputting primer pairs into UCSC 
	isPCR program 
	'''
	if delimiters is None:
		delimiters = "\t"

	reader = csv.reader(open(input_file),"r+"), delimiter="\t")

	# skip first line containing column information
	header = reader.next()

	for primer in reader:
    		ispcr = urllib2.urlopen("http://genome.ucsc.edu/cgi-bin/hgPcr?db="+DB+\
    		"&wp_target=genome&wp_f="+primer[1].upper()+"&wp_r="+primer[2].upper()+\
    		"&wp_size=4000&wp_perfect=15&wp_good=15&boolshad.wp_flipReverse=0")
    		search = ispcr.read()
    
    		if re.search(r'No matches.*',search):
        		output = primer[0]+"\t"+primer[1]+"\t"+primer[2]+"\t"+"none"+"\t"+"no match"+"\n"
        		print output
                	
		elif re.search(r'No matches.*',search) is None:
        		# find the begining of the amplicon sequence and everything that occurs after
		        # it, hence re.DOTALL
			product_seq = re.findall(r"[ACTGatcg]{25,54}.*",search, re.DOTALL)
        		# remove the special characters []'. \n from product_seq
        		no_special_chars = re.sub(r'([\[\]\'\, ,\\n])',"",str(product_seq))    
        		# the amplicon seq begins and ends with upper case bases and the bases
        		# between are lower case, therefore the following code will capture the
        		# amplicon sequence only
        		amplicon = re.search(r"[ACTG]{10,40}[actg]{25,1000}[ACTG]{10,40}",\
        		no_special_chars, re.DOTALL).group()
        
       	 		# simple math to get the GC %        
       	 		GC_percent = ((amplicon.count("G")+amplicon.count("C")+amplicon.count("c")+\
        		amplicon.count("g"))/ len(amplicon)) * 100

        		print GC_percent
	        	print len(amplicon)
