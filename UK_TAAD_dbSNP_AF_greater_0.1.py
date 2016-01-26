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
                if AF != "." and float(AF) > 0.1: 
                    dbSNP_UK_Matched_Common.append(pos+"\t"+rs+"\t"+str(AF)+"\t"+alt)
    
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
                        if len(cleaned_cell_values) == 5 and "1/1" in cleaned_cell_values[0] and int(cleaned_cell_values[2]) > 49:    # must have a 1/1 genotype and alternative AD > 49
                            matching_rows.append(pos_alt_pair)  # append the (pos,alt) associated with each individually filtered sample info field

    matching_rows_count = collections.Counter(matching_rows)   # count no. of duplicates in list; no. of variants common in our cohort
    return matching_rows_count#[('20:45354291', 'A')]



          
def process_files(file_name_1, file_name_2):
    """Returns rows from the second file that match both pos/alt and a valid cell value"""
    
    dbSNP_UK_Matched_Common = extract_relevant_info(file_name_1)
    valid_pos_alt_pairs = get_pos_alt_pairs(dbSNP_UK_Matched_Common)
    matching_rows_count = get_matching_rows(file_name_2, valid_pos_alt_pairs)
    
    most_common_variants = []    
    
    for pos,alt in matching_rows_count.most_common():
        most_common_variants.append(str(pos[0]) +","+str(pos[1])+"," + str(alt))

    return most_common_variants
    
    
    
    
    
print process_files("UK_TAAD_Variants_in_dbSNP.vcf","test.txt")
