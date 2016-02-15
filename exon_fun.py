import requests, json
from get_variant_information import *

# Changes to be made

# 3-exception handeling
# 4-structure
# 5-function comment clarity



def get_exon_number(input_file):
    '''From a list of variant aliases, get its row number to allow us to retrieve
       its transcript id and variant position from the variant database. Parse this
       information into generate_exon_numbering to return the exon number in which
       the variant resides within.
       
       The input file should consist of variant_aliases seperated by a newline
       
       Use the get_variant_information module to retireve the aformentioned 
       row information. Below is a tree to show the codes structure:
       
       get_exon_number
                |
            generate_exon_numbering
                    |
                request_ensembl--->all_exon_regions--->get_exon_id--------
                   ----->exon_num & last_exon ------------------------> OUT   
                    |___ intron_region ---> intron_num & last_exon ---> OUT
    '''
    
    exon_info = open("exon_info.txt","w")
    exon_info.write("Variant_alias"+"\t"+ 
                    "pos"+"\t"+"transcript_id"+"\t"+
                    "intron_num"+"\t"+"exon_num"+"\n")
    
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
            
            exon_info.write(generate_exon_numbering(variant_alias, transcript_hgvs, pos)+"\n")
            print generate_exon_numbering(variant_alias, transcript_hgvs, pos)
            
        except Exception,e:                                 # BAD!
            exon_info.write(variant_alias +"\t"+"Something Happened"+"\n")
            print variant_alias +"\t"+"Something Happened"
            
            continue

    exon_info.close()
    
    

def generate_exon_numbering(variant_alias,transcript_hgvs,var_position):
    ''' Return the exon number from which the variant is within and the total
        number of exons in the given transcript. If no exon_id is found, then 
        retrieve an intron number for the variant.
        
        This essentially links all the below functions to generate an exon or intron number
    '''
    transcript = transcript_hgvs.split(":")[0][:-2]
    pos = var_position
    
    exon_dics = request_ensembl(transcript)                 # request from REST API
    exon_region = all_exon_regions(exon_dics,transcript)# get all exon_id, start and stop exon pos
    exon_id = get_exon_id(exon_region,pos) # filter for exon id in which variant is within
    if not exon_id:
        sorted_exon_regions = sorted(exon_region)
        intron_region = all_intron_regions(sorted_exon_regions)
        intron_num = intron_number(intron_region,pos)
        
        if intron_num is None:
            return variant_alias+"\t"+pos+"\t"+transcript+"\t"+"NO INTRON/EXON MATCHED"
        else:
            last_exon = total_exons(exon_dics,transcript)
            last_intron = int(last_exon)-1
            return variant_alias+"\t"+pos+"\t"+transcript+"\t"+str(intron_num)+"/"+str(last_intron)+"\t"+"-"
    else:
        exon_num = exon_number(exon_dics,exon_id,transcript)    # use th exon_id to get exon number
        last_exon = total_exons(exon_dics,transcript)           # get total exons of transcript
        
        
        return variant_alias+"\t"+pos+"\t"+transcript+"\t"+"-"+"\t"+str(exon_num)+"/"+str(last_exon)



def request_ensembl(transcript):
    ''' Retrieve all expanded exon information associated with the transcipt ID
        from the Ensembl REST API.
        
        The resulting list of dicts will be stored as a global variable
    '''
    try:
        url = "http://grch37.rest.ensembl.org/overlap/id/"
        ext = "?feature=exon;content-type=application/json;expand=1"
        req = requests.get(url+transcript+ext)
        req.raise_for_status()
        exon_dics = json.loads(req.text)
        return exon_dics
    except requests.exceptions.RequestException as e:
        return e
        
        
def all_exon_regions(exon_dics,transcript):
    '''Get every Exon_ID and its start and stop position which reside within 
       the inputted transcript ID 
       
       A conditional is used to extract only the dictionarys containing the 
       transcript ID. The start, stop positions and exon_id for each matching 
       dictionary is then stored in a tuple named exon_region and returned. 
    '''
    exon_region = set()
    for exon_dicts in exon_dics:
        if exon_dicts.get("Parent") == transcript:
            rank = exon_dicts.get("rank")
            start = exon_dicts.get("start")
            end = exon_dicts.get("end")
            exon = exon_dicts.get("exon_id")
            exon_start_end = (rank,exon,start,end)
            exon_region.add(exon_start_end)
    
    return exon_region
    
def all_intron_regions(sorted_exon_regions):
    ''' Get all intron numbers, start and stop positions
    ''' 
    intron_region = set()
    
    for i in range(0,len(sorted_exon_regions)):
        for x in range(0,len(sorted_exon_regions)):
            if x == i +1:
                intron_number = str(sorted_exon_regions[i][0])
                end_pos_previous_exon = str(sorted_exon_regions[i][3])
                start_pos_next_exon = str(sorted_exon_regions[x][2])
                intron_info = (int(intron_number),end_pos_previous_exon,\
                start_pos_next_exon,str(int(end_pos_previous_exon)-int(start_pos_next_exon)))
                intron_region.add(intron_info)
    
    return intron_region

def intron_number(intron_region,pos):
    ''' Returns the intron number in which the variant is within
    '''
    for intron in intron_region:
        if int(intron[3]) > 0:
            for x in range(int(intron[2]),int(intron[1])):   # no idea why the start and end pos are swapped, hence the need for if/else
                if x == int(pos.split(":")[1]):              # good candidate for recursion
                    intron_num= intron[0]
                    return intron_num
        else:
            for x in range(int(intron[1]),int(intron[2])):
                if x == int(pos.split(":")[1]):
                    intron_num= intron[0]
                    return intron_num  

            
                
def get_exon_id(exon_region,pos):  
    ''' Get the exon id for which the variant position is within.
    
        Use the exon_region tuple to iterate between the start and stop 
        positions of each exon within the transcript until one number 
        matches the variant position number, then return the exon_id 
        associated with that position.
    '''  
    exon_id = []      
    for exons in exon_region:
        for x in range(exons[2],exons[3]):
            if x == int(pos.split(":")[1]):
                exon_id.append(exons[1])
    
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
#print get_exon_number("HUK10-BO")

