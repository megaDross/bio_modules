''' Gets Exon/Intron number (pyensembl cannot do this). It is very rough, old and inefficient. Consider rewriting.'''
import requests
import sys

# NOTE: Response 400 in request_ensembl() occurs when a transcript has NO exons.
#       Try wrapping RequestException in a custom exception and using that as a 
#       trigger to return (-, -, -).

def main(transcript, hg_version, pos):
    ''' Get exon information from a genomic position.

    Args:
        transcript: ensembl transcript ID
        hg_version: human genome version (hg19 or hg39)
        pos: genomic coordinate

    Returns:
        (Ensembl Exon ID, Exon Number, Intron Number)
    '''
    exon_dics = request_ensembl(hg_version, transcript)
    exon_region = all_exon_regions(transcript, exon_dics)
    exon_id = get_exon_id(pos, exon_region) # filter for exon id in which variant is within
    if not exon_id:
        sorted_exon_regions = sorted(exon_region)
        intron_region = all_intron_regions(sorted_exon_regions)
        intron_num = intron_number(pos, intron_region)
        if intron_num is None:
            return ('-', '-',"-")
        else:
            last_exon = total_exons(transcript, exon_dics)
            last_intron = int(last_exon)-1
            return ("-",str(intron_num)+"/"+str(last_intron),"-")
    else:
        exon_num = exon_number(transcript, exon_dics,exon_id)    # use th exon_id to get exon number
        last_exon = total_exons(transcript, exon_dics)           # get total exons of transcript
    return (exon_id[0], "-" ,str(exon_num)+"/"+str(last_exon))

def request_ensembl(hg_version, transcript):
    ''' Retrieve all expanded exon information associated with a given transcript ID
        in JSON format
    '''
    try:
        if hg_version == "hg19":
            grch = "grch37."
        elif hg_version == "hg38":
            grch = ""
        else:
            # NOTE: raise a custom exception instead of printing
            print("Incompatible human genome version")
            sys.exit(1) 
        url = "".join(("http://",grch,"rest.ensembl.org/overlap/id/", transcript,
                       "?feature=exon;content-type=application/json;expand=1"))
        req = requests.get(url)
        req.raise_for_status()
        exon_dics = req.json()
        return exon_dics
    
    except requests.exceptions.RequestException as e:
        if str(e.response) == '<Response [400]>':
            return e.response

def all_exon_regions(transcript, exon_dics):
    '''Get every exon ID and its start and stop position which resides within the
       given transcript
    '''
    exon_region = set()
    for exon_dicts in exon_dics:
        if exon_dicts.get("Parent") == transcript:
            rank = exon_dicts.get("rank")
            start = exon_dicts.get("start")
            end = exon_dicts.get("end")
            exon = exon_dicts.get("exon_id")
            exon_start_end = (rank, exon, start, end)
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
                intron_info = (int(intron_number),
                               end_pos_previous_exon,
                               start_pos_next_exon,
                               str(int(end_pos_previous_exon)-int(start_pos_next_exon)))
                intron_region.add(intron_info)
    return sorted(intron_region)

def intron_number(pos, intron_region):
    ''' Returns the intron number in which the variant is within
    '''
    for intron in intron_region:
        if int(intron[3]) > 0:
            # no idea why the start and end pos are swapped, hence the need for if/else
            for x in range(int(intron[2]),int(intron[1])):          
                if x == int(pos.split(":")[1]):
                    intron_num= intron[0]
                    return intron_num
        else:
            for x in range(int(intron[1]),int(intron[2])):
                if x == int(pos.split(":")[1]):
                    intron_num= intron[0]
                    return intron_num  
                        
def get_exon_id(pos, exon_region):  
    ''' Get the exon id for which the variant position is within.
    
        Use the exon_region tuple to iterate between the start and stop 
        positions of each exon within the transcript until one number 
        matches the variant position number, then return the exon_id 
        associated with that position.
    '''  
    exon_id = []      
    for exons in sorted(exon_region):
        for x in range(exons[2],exons[3]):
            if x == int(pos.split(":")[1]):
                exon_id.append(exons[1])
    return exon_id
    
def exon_number(transcript,exon_dics,exon_id):
    ''' Returns the exon number in which the variant is within.
    
        use the exon_id to retrieve the matching dictionary and extract
        the rank (exon number) from that dictionary.
    '''
    for y in exon_dics:
        if y.get("exon_id") == exon_id[0] and y.get("Parent") == transcript:
            exon_number = y.get("rank")
    return exon_number

def total_exons(transcript,exon_dics):
    ''' Get the last exon number in the transcript
    
        List comprehension which extracts all exon numbers found in all the 
        dictionarys with the transcript id. Extract and return the last element
        as the total_exons in the generated list.
    '''
    all_exon_ranks = [i.get("rank") for i in exon_dics 
                      if i.get("Parent") == transcript]
    total_exons = all_exon_ranks[-1]
    return total_exons
