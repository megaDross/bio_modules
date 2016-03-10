import re, os.path


#exception handeling and use string as input

def primer_finder(input_file,primer_database,output_file=None,delimiters="\t"):
    ''' Takes variant postion(s) and matches it with a primer in the primer 
        database. 
       
        input_file should be in the following format (tab deliminated):
            variant_name    variant_position
        
        primer_database should be in the following format (tab deliminated):
            primer name    F-primer    R-primer    product size    genomic region 
        
        output_file is optional and delimiters is tab by default
    '''
    primer_file = open(primer_database)
    all_primer_pos = get_all_primer_pos(primer_file)
    
    if os.path.isfile(input_file) is False:
        var_pos = re.sub(r'[^0-9:]','',input_file)
        matched_primers = match(var_pos,all_primer_pos)
        return matched_primers
            
    if output_file is not None:
        output_file = open(output_file,"w")
        output_file.write("Variant_Name"+"\t"+"Primer_Name")
        
        var_file = open(input_file)
        for var in var_file:
            var_name = var.split("\t")[0]
            var_pos = re.sub(r'[^0-9:]','',var.split("\t")[1])
            matched_primers = match(var_pos,all_primer_pos,var_name)
            
            print matched_primers
            output_file.write("".join((matched_primers,"\n")))
        output_file.close()
    
    else:    
        var_file = open(input_file)
        for var in var_file:
            var_name = var.split("\t")[0]
            var_pos = re.sub(r'[^0-9:]','',var.split("\t")[1])
            matched_primers = match(var_pos,all_primer_pos,var_name)
            print matched_primers
            
    
        
        
def get_all_primer_pos(primer_file):
    ''' Generate evey genomic position within each primer and concatenate it with
        the primer name and pack within a list; each list contains every position 
        within a primer pair. Append every list to an empty list and unpack the list of
        lists with a nested list comprehension.
    '''
    header = primer_file.next()
    all_primer_pos = []
    
    for i in primer_file:
        
        i = i.split("\t")
        primer_name = i[0]
        primer_range = re.split(r'[:-]',i[4])   # split into three components
        chrom = re.sub(r'[^0-9]','',primer_range[0])  # remove any non-integers
        start = int(primer_range[1])
        stop = int(primer_range[2])
        
        all_pos_in_primer = [ " ".join((primer_name,chrom+":"+str(i))) 
                            for i in range(start,stop)]
        all_primer_pos.append(all_pos_in_primer)
        
    all_primers_pos_unpacked = [x for i in all_primer_pos for x in i]   # unpack list of lists
    return all_primers_pos_unpacked
        
        
        
def match(var_pos,primer_info,var_name=None):
    ''' Match a given variant position against every genomic position covered
        in all the primers in the primer database.
    '''
    
    for i in primer_info:
        
        if var_name is None:   # used if a string given as input instead of a file
            var_name = var_pos
        
        primer_name = i.split(" ")[0]
        primer_pos = i.split(" ")[1]
        
        
        ## if no match then print cant find matching primer
        ## proving to be challenging, not working as expected
        ## if not Non, if != "" if not "" etc. not working, just returns None, regardless
        
        if var_pos == primer_pos:
            answer = ",".join((var_name,primer_name))
            return answer
                 
        
    
    
    

print primer_finder("2:189859300","TAAD_Primer_Validation_Database.txt")
#print primer_finder("in.txt","TAAD_Primer_Validation_Database.txt")
