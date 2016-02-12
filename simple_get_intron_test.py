import requests, json
from get_variant_information import *

transcript = "ENST00000552516"
pos = "chr9:101907012"


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
            start = exon_dicts.get("start")
            end = exon_dicts.get("end")
            exon = exon_dicts.get("exon_id")
            rank = exon_dicts.get("rank")
            exon_start_end = (rank,exon,start,end)
            exon_region.add(exon_start_end)
            
    return exon_region

def all_intron_regions(exon_region):
    ''' Get all intron numbers, start and stop positions
    ''' 
    intron_region = set()
       
    for i in range(0,len(sorted_tuple)):
        for x in range(0,len(sorted_tuple)):
            if x == i +1:
                intron_number = str(sorted_tuple[i][0])
                end_pos_previous_exon = str(sorted_tuple[i][3])
                start_pos_next_exon = str(sorted_tuple[x][2])
                intron_region.add(intron_number +"\t"+end_pos_previous_exon +"\t"+start_pos_next_exon)
    
    return intron_region

def intron_number(intron_region,pos):
    ''' Returns the intron number in which the variant is within
    '''
    for introns in intron_region:
        intron = introns.split("\t")
        for x in range(int(intron[1]),int(intron[2])):
            if x == int(pos.split(":")[1]):
                intron_number= intron[0]
    return intron_number
                
    

recieve_tuple = all_exon_regions(request_ensembl(transcript),transcript)
sorted_tuple = sorted(recieve_tuple)
all_intron_region = all_intron_regions(sorted_tuple)
number =intron_number(all_intron_region,pos)
print number