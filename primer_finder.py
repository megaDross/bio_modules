#!/usr/bin/python
import re
import csv
import itertools


def primer_finder(input_file, primer_database, output_file, delimiters=None):
	'''Takes an input file consisting of variant positions (input_file) and finds
	primers within a primer_database which cover said variant. Allowing one to find
	possible primers and outputs them into a new_file.

	Input file should be in he following format:
		variant_name	Variant_position
	primer_database format:
		Primer_name	F-primer	R-primer	Product_size	Region
	
	The delimiters argument is defaulted to tab
	'''
	####################### PREPERATION ############################
	
	# change dir to where the txt/csv files reside and open the files that contain primer
	# information and one which conatins new variants and variant position. use csv module
	# to split columns by tab and allow one to manipulate each column. open two empty
	# lists, one empty dictionaries and writable empty file.
	
	if delimiters is None:
	   delimiters = "\t"

	primer_file = open(primer_database,"r")
	variant_file = open(input_file,"r")
	
	product_file =open(output_file,"w")
	
	primers_pos = csv.reader(primer_file, delimiter="\t")
	variant_pos = csv.reader(variant_file, delimiter="\t")
	
	var_area = {}
	primer_area = {}
	primerz = []
	product = []
	
	
	#################### PRIMER POSITION LIST #######################
	
	# select fifth column of primer info file (region) and only select if begins with "chr".
	# split into three variables chromosome, primer start and end site. create a new loop
	# which produces a int for every number between start and stop site for each primer pair
	# with the the chromosome number prefixed to each number and append them to the empty
	# primerz list, to create a list element for every genome position in which a primer pair
	# covers. 
	for region in primers_pos:
    	   #prim = region[0]
    	   #reg = region[4]
    	   #primer_area[prim] = reg
    	   if re.search(r'chr', region[4]):
	       chromosome = re.split(r'[:-]',region[4])[0]
	       start = int(re.split(r'[:-]',region[4])[1])
	       stop = int(re.split(r'[:-]',region[4])[2])
	       for numbers in range(start,stop):
                   	primer_area = chromosome + ":" + str(numbers)
                   	primerz.append(region[0] + " " + primer_area)     
	
            	
	#################### VARIANT DICTIONARY ##############################
	
	# prodce a dictionary that contains all varinat names as keys and the equivalent 
	# ivar poisition as items and store them within the empty dictionary.
	for var in variant_pos:
           varz = var[1]
    	   name = var[0]
    	   var_area[name] = varz
    	
    	
	
    	
	
	#################### PRIMER FINDER ####################################
	
	# loop both the newly created primerz list and the newly created dictionary, with the keys
	# placed in the variable var_name and the items in the variable var_pos. If the var_pos is found
	# in the list of primers then append the variants name along with the primer pair name that
	# can be used to produce an amplicon for said variant name.
	# 
	# itertools allows one to produce loops with multiple variables. itertools.product is like zip
	# expect the two loops don't have to be the same length and unlike itertools.izip it doesn't
	# produce NoneTypes for extra entries in shortest lengthed loop. Instead it loops the shortest length
	# loop to the same size as the longest length loop (print var_name in the loop to understand)
	#
	# produce a loop that prints product list minus duplicate entries (set) and also writes them to the
	# empty file.
	for (primer_pos),(var_name,var_pos) in itertools.product(primerz,var_area.items()):
    	   if re.search(var_pos,primer_pos):
	       matched_primer_pos = primer_pos.split()[0]
	       product.append("variant: "+var_name +", "+ " primers: " + str(matched_primer_pos))
	
	
	
	for match in list(set(product)):
    	   print match
    	#matche_primer = re.split(":",match)[2][1:]
    	   product_file.write(match+"\n")
        product_file.close()
	
	
	#for that in product:
 	#   print re.split(":",that)[2][1:]
	
	
	
	# Below shows how one can find matching entries in two lists using set().intersection()
	a = ["red","yellow","blue"]
	b = ["yellow","green","purple","red"]
	
	intersect = list(set(a).intersection(b))
	print "\n" + str(intersect)
	
	
