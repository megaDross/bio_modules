import requests, json
from get_variant_information import *

# Changes to be made

# 2-find a way to get intron names/nums
# 3-exception handeling
# 4-structure
# 5-function comment clarity



def get_exon_number(input_file):
    '''From a list of variant aliases, get its row number to allow us to retrieve
       its transcript id and variant position from the variant database. Parse this
       information into generate_exon_numbering to return the exon number in which
       the variant resides within.
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
            print variant_alias +"\t"+"Something Happened >:("
            continue



def generate_exon_numbering(variant_alias,transcript_hgvs,var_position):
    ''' Return the exon number from which the variant is within and the total
        number of exons in the given transcript
    '''
    
    transcript = transcript_hgvs.split(":")[0][:-2]
    pos = var_position
    
    
    exon_dics = request_ensembl(transcript)                 # request from REST API
    exon_region = get_exon_information(exon_dics,transcript)# get all exon_id, start and stop exon pos
    exon_id = get_exon_id(exon_region,pos)                  # filter for exon id in which variant is within
    exon_num = exon_number(exon_dics,exon_id,transcript)    # use th exon_id to get exon number
    last_exon = total_exons(exon_dics,transcript)           # get total exons of transcript
    
    return variant_alias+"\t"+pos+"\t"+transcript+"\t"+str(exon_num)+"/"+str(last_exon)



def request_ensembl(transcript):
    ''' Retrieve all expanded exon information associated with the transcipt ID
        from the Ensembl REST API.
        
        The resulting list of dicts will be stored as a global variable
    '''
    
    url = "http://grch37.rest.ensembl.org/overlap/id/"
    ext = "?feature=exon;content-type=application/json;expand=1"
    req = requests.get(url+transcript+ext)
    req.raise_for_status()
    exon_dics = json.loads(req.text)
    return exon_dics
    
    
def get_exon_information(exon_dics,transcript):
    '''Get every Exon_ID and its start and stop position which reside within 
       the inputted transcript ID 
       
       A conditional is used to extract only the dictionarys containing the 
       transcript ID. The start, stop positions and exon_id for each matching 
       dictionary is then stored in a tuple named exon_region and returned. 
    '''
    exon_region = set()
    for exon_dicts in exon_dics:
        if exon_dicts.get("Parent") == transcript:
            start = exon_dicts.get("start")
            end = exon_dicts.get("end")
            exon = exon_dicts.get("exon_id")
            exon_start_end = (exon,start,end)
            exon_region.add(exon_start_end)
            
    return exon_region
    
    

def get_exon_id(exon_region,pos):  
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
    


def exon_number(exon_dics,exon_id,transcript):
    ''' Returns the exon number in which the variant is within.
    
        use the exon_id to retrieve the matching dictionary and extract
        the rank (exon number) from that dictionary.
    '''
    for y in exon_dics:
        if y.get("exon_id") == exon_id[0] and y.get("Parent") == transcript:
            exon_number = y.get("rank")
   
    return exon_number

        
    
def total_exons(exon_dics,transcript):
    ''' Get the last exon number in the transcript
    
        List comprehension which extracts all exon numbers found in all the 
        dictionarys with the transcript id. Extract and return the last element
        as the total_exons in the generated list.
    '''
    all_exon_ranks = [i.get("rank") for i in exon_dics if i.get("Parent") == transcript]
    total_exons = all_exon_ranks[-1]
    return total_exons


#l = ["WYN11-XX-62","WYN10-XX-00","WYN8-XX-03"]

print get_exon_number("all_vars_alias.txt")


