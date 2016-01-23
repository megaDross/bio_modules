"""
utilise various UCSC tools via the command line
"""
from __future__ import division
import urllib2, re, csv

# unknown_primer wont write to the file, file empty after processing

def unknown_primer(DB, input_file, output_file, delimiters="\t"):
	
	test = open(output_file,"w")
	read_input = open(input_file,"r+")
	skip_header = read_input.next()
	ispcr_results = []
	
	for f in read_input:
	        primer = f.split(delimiters)
		ispcr = urllib2.urlopen("http://genome.ucsc.edu/cgi-bin/hgPcr?db="+DB+\
		"&wp_target=genome&wp_f="+primer[1]+"&wp_r="+primer[2]+\
		"&wp_size=4000&wp_perfect=15&wp_good=15&boolshad.wp_flipReverse=0")
		search = ispcr.read()
		ispcr_results.append(search)
                print ispcr_results
                
print unknown_primer("hg19","input_UCSC.txt","output_txt")
		