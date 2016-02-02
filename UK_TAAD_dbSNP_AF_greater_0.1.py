import re, collections
# tabix <dbSNP.vcf> <positions>
# the tabixed file is the UK_vars vcf which was filtered through with dbSNP vcf

def extract_relevant_info(file_name):
    '''extracts the relevant information from the tabixed output file
    '''
    dbSNP_UK_Matched_Common = []
    for row in open(file_name, 'r'):
        
        info = row.split("\t")[7]
        pos = row.split("\t")[0]+":"+row.split("\t")[1]
        rs = row.split("\t")[2]
        alt = row.split("\t")[4] 
        
        search_CAF = re.findall(r"CAF=.*;",info)
        CAF = ','.join(search_CAF).split(",")   # comma seperated allele frequency (ref,alt)
        
        if len(CAF) > 1:
                AF = re.sub(";","",CAF[1])   # alternative allele frequency 
                if AF != "." and float(AF) > 0.29: 
                    dbSNP_UK_Matched_Common.append(pos+"\t"+rs+"\t"+str(AF)+"\t"+alt.split(",")[0])
                    
                    ## THERE MUST BE A BETTER WAY TO DO THE MULTIPLE IF STATEMENTS BELOW
                    
                    if len(alt.split(",")) > 1:   # if the alt variable has more than one value, use the second
                        dbSNP_UK_Matched_Common.append(pos+"\t"+rs+"\t"+str(AF)+"\t"+alt.split(",")[1])
                    if len(alt.split(",")) > 2:   #if the alt variable has more than one value, use the third
                        dbSNP_UK_Matched_Common.append(pos+"\t"+rs+"\t"+str(AF)+"\t"+alt.split(",")[2])
                    if len(alt.split(",")) > 3:   #if the alt variable has more than one value, use the fourth
                        dbSNP_UK_Matched_Common.append(pos+"\t"+rs+"\t"+str(AF)+"\t"+alt.split(",")[3])
    
    return dbSNP_UK_Matched_Common
    
    
    
def get_pos_alt_pairs(file_name):
    """Returns a set of the unique pos/alt pairs from all rows in the given file"""
    
    pairs = set()     # eliminate replicates 
    
    for row in file_name:
            pos_alt = row.split('\t')[0:4:3]   #extract 0 and 3 indices
            pos_alt_pair = tuple(map(str.strip,pos_alt))  # remove newlines from the list
            pairs.add(pos_alt_pair)   #append these pairs into empty set
    
    return pairs
    

def get_matching_rows(file_name, valid_pos_alt_pairs):
    """Returns a collection of matching rows.

    Compares the pos/alt from each row in the given file against a master pos/alt list.
    If a match is found, the row's fields are further validated to check if it
    is a row of interest.  If so, it gets added to the returned collection.
    """

    matching_rows = []
   
    with open(file_name, 'r') as input_file:
        for row in input_file:
            if not row.startswith('#'):    # removes headers
                
                cells = row.split('\t')
                pos = '{0}:{1}'.format(cells[0], cells[1])
                alt = cells[4]
                pos_alt_pair = (pos, alt)
                
                if pos_alt_pair in valid_pos_alt_pairs:
                    cell_values = [cell.split(':') for cell in cells[16:]]    # this creates a list of lists containing sample info fields
                    remove_empty_cells = filter(lambda a: "./." not in a, cell_values)
                    
                    for cleaned_cell_values in remove_empty_cells:
                        if len(cleaned_cell_values) == 5:
                            
                            DP = int(cleaned_cell_values[2])
                            GT = cleaned_cell_values[0]
                            AD_ref = int(cleaned_cell_values[1].split(",")[0])
                            AD_alt = int(cleaned_cell_values[1].split(",")[1])
                            
                            GQ = int(cleaned_cell_values[3])
                            
                            if (GT == "1/1" or GT == "2/2" or GT == "3/3" or GT == "4/4") and DP > 49 and GQ >29:    
                                matching_rows.append(pos_alt_pair)  # append the (pos,alt) associated with each individually filtered sample info field
                                print cleaned_cell_values
                                
                            if (GT == "0/1" or GT == "1/0") and AD_alt > 49 and GQ > 29:
                                matching_rows.append(pos_alt_pair)
                                print cleaned_cell_values
                                
                            if (GT == "2/0" or GT =="0/2"):
                                AD_second_alt = int(cleaned_cell_values[1].split(",")[2])
                                if AD_second_alt > 49 and GQ > 29:
                                    matching_rows.append(pos_alt_pair)
                            
                            if (GT == "3/0" or GT == "0/3"):
                                AD_third_alt = int(cleaned_cell_values[1].split(",")[3])
                                if AD_third_alt > 49 and GQ > 29:
                                    matching_rows.append(pos_alt_pair)
                            
                            if (GT == "4/0" or GT == "0/4"):
                                AD_fourth_alt = int(cleaned_cell_values[1].split(",")[4])
                                if AD_fourth_alt > 49 and GQ > 29:
                                    matching_rows.append(pos_alt_pair)
                            
                                
                            
                                
                            

    matching_rows_count = collections.Counter(matching_rows)   # count no. of duplicates in list; no. of variants common in our cohort
    return matching_rows_count#[('20:45354291', 'A')]



          
def process_files(file_name_1, file_name_2):
    """Returns rows from the second file that match both pos/alt and a valid cell value"""
    output = open("hopeful_output.txt","w")
    output.write("variant_pos"+"\t"+"ALT"+"\t"+"sample_number"+"\n")
    
    
    dbSNP_UK_Matched_Common = extract_relevant_info(file_name_1)
    valid_pos_alt_pairs = get_pos_alt_pairs(dbSNP_UK_Matched_Common)
    matching_rows_count = get_matching_rows(file_name_2, valid_pos_alt_pairs)
    
    most_common_variants = []    
    
    for pos,alt in matching_rows_count.most_common():
        most_common_variants.append(str(pos[0]) +","+str(pos[1])+"," + str(alt))
        output.write(str(pos[0]) +"\t"+str(pos[1])+"\t" + str(alt)+"\n")
    
    output.close()
    return most_common_variants
    
    
    
    
    
print process_files("UK_TAAD_Variants_in_dbSNP.vcf","var.both.TAAD_UK_Nov2015.filters.vcf")
