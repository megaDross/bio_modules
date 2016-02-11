import requests, json
from get_variant_information import *

# Changes to be made

# 1-get rid of global variables
# 2-find a way to get intron names/nums
# 3-exception handeling
# 4-structure
# 5-function comment clarity



def do_it(input_file):
    '''From a list of variant aliases, get its row number to allow us to retrieve
       its transcript id and variant position from the variant database. Parse this
       information into generate
    '''
    if input_file.endswith(".txt"):
            var_alias_list = [ var.rstrip() for var in open(input_file)]
    else:
            var_alias_list = [input_file]
            
    variants = open_variant_validation_spreadsheets()[0]
    
    for variant_alias in var_alias_list:
        
        try:
            variant_row = get_row(variant_alias,variants)
            variant_info = get_variant_info(variant_row,variants,variant_header)
            pos = variant_header.get("Variant_Position")
            transcript_hgvs = variant_header.get("HGVSc")
            
            print generate_exon_numbering(variant_alias, transcript_hgvs, pos)
            
        except:                                 # BAD!
            print variant_alias +"\t"+"Nope"
            continue



def generate_exon_numbering(variant_alias,transcript_hgvs,var_position):
    ''' 
    '''
    global transcript
    global pos 
    
    transcript = transcript_hgvs.split(":")[0][:-2]
    pos = var_position
    
    
    generate_all_exon_data = request_ensembl(transcript)
    exon_region = get_exon_information()
    exon_id = get_exon_id(exon_region)
    exon_num = exon_number(exon_id)
    last_exon = total_exons()
    
    return variant_alias+"\t"+pos+"\t"+transcript+"\t"+str(exon_num)+"/"+str(last_exon)



def request_ensembl(transcript):
    ''' Retrieve all expanded exon information associated with the transcipt ID
        from the Ensembl REST API.
        
        The resulting list of dicts will be stored as a global variable
    '''
    global all_exon_data
    
    url = "http://grch37.rest.ensembl.org/overlap/id/"
    ext = "?feature=exon;content-type=application/json;expand=1"
    req = requests.get(url+transcript+ext)
    req.raise_for_status()
    all_exon_data = json.loads(req.text)
    
    
    
def get_exon_information():
    '''Get every Exon_ID and its start and stop position which reside within 
       the inputted transcript ID 
       
       A conditional is used to extract only the dictionarys containing the 
       transcript ID. The start, stop positions and exon_id for each matching 
       dictionary is then stored in a tuple named exon_region and returned. 
    '''
    exon_region = set()
    for exon_dicts in all_exon_data:
        if exon_dicts.get("Parent") == transcript:
            start = exon_dicts.get("start")
            end = exon_dicts.get("end")
            exon = exon_dicts.get("exon_id")
            exon_start_end = (exon,start,end)
            exon_region.add(exon_start_end)
            
    return exon_region
    
    

def get_exon_id(exon_region):  
    ''' Get the exon id for which the variant position is within.
    
        Use the exon_region tuple to iterate between the start and stop 
        positions of each exon within the transcript until one number 
        matches the variant position number, then return the exon_id 
        associated with that position.
    '''  
    exon_id = []      
    for exons in exon_region:
        for x in range(exons[1],exons[2]):
            if x == int(pos.split(":")[1]):
                exon_id.append(exons[0])
    
    return exon_id
    


def exon_number(exon_id):
    ''' Returns the exon number in which the variant is within.
    
        use the exon_id to retrieve the matching dictionary and extract
        the rank (exon number) from that dictionary.
    '''
    for y in all_exon_data:
        if y.get("exon_id") == exon_id[0] and y.get("Parent") == transcript:
            exon_number = y.get("rank")
   
    return exon_number

        
    
def total_exons():
    ''' Get the last exon number in the transcript
    
        List comprehension which extracts all exon numbers found in all the 
        dictionarys with the transcript id. Extract and return the last element
        as the total_exons in the generated list.
    '''
    all_exon_ranks = [i.get("rank") for i in all_exon_data if i.get("Parent") == transcript]
    total_exons = all_exon_ranks[-1]
    return total_exons



    

print do_it("all_vars_alias.txt")


