from Bio import SeqIO
import os
import  subprocess


def handle_seq_file(seq_file, seq_dir):
    ''' Decide how to process seq_files
    
        -id is used instead of -if as the -if only works on a file containing paths to the ab1 files requiring processing; cannot directly select a file in ttuner.
    '''
    ttuner = "/home/david/bin/tracetuner_3.0.6beta/rel/Linux_64/ttuner"

    if seq_file.endswith(".seq"):
        seq_file = "".join([x.strip("\n") for x in open(seq_file) 
                                          if not x.startswith(">")])

    elif seq_file.endswith("ab1") and os.path.isfile(seq_file+".seq") \
        and os.path.isfile(seq_file+".tab"):
        seq_file = "".join([x.rstrip("\n") for x in open(seq_file+".seq")
                                        if not x.startswith(">")]) 

    elif seq_file.endswith("ab1") and not os.path.isfile(seq_file+".tab"):
        subprocess.call([ttuner, "-tabd", seq_dir, "-id", 
                         seq_dir, "-mix"], stdout=open(os.devnull, 'wb'),
                         stderr=open(os.devnull, 'wb'))
        subprocess.call([ttuner, "-sd", seq_dir, "-id", seq_dir], 
                        stdout=open(os.devnull, 'wb'), stderr=open(os.devnull, 'wb'))

    else:
        seq_file = seq_file

    return seq_file




def get_matching_seq_file(query, directory):
    ''' Find a file name that best matches given query and return the
        sequence string in matched file
    '''
    length = len(query)
    cut = query.replace("_","-").split("-")
    reversed_query = "_".join((cut[1],cut[0])) if len(cut) > 1 else query
    
    # stop recursion
    if length < 6:
        return []

    # appending is required in case there is more than one match, returns first match
    store_matches = []
    for f in [x for x in os.listdir(directory) if x.endswith(".ab1")]:
        if query.replace("-","_") in f or query.replace("_","-") in f or \
           reversed_query in f or reversed_query.replace("_","-") in f:
            file_match = directory+f
            store_matches.append(file_match)
    
    # trying to find matches that have both _ and -
    sorted_matches = sorted(store_matches)

    if not sorted_matches:
        return get_matching_seq_file(query[:-1], directory)

    return sorted_matches





def convert_ab1_to_seq(ab1):
    ''' Extract the sequence string from .ab1 file and return
        
        INCOMPATIBLE WITH TTUNER!!!!!!!
        WARNING: This method needs ALOT of testing
    '''
    handle = open(ab1, "rb")
    new_file_name = ab1.replace("ab1", "seq")
    new_file = open(new_file_name, "w")

    # Use SeqIO and re to get the sequence from the ab1 file
    for record in SeqIO.parse(handle, "abi"):
            raw_data = record.annotations.get("abif_raw")
            sequence = raw_data.get("PBAS2")
            new_file.write(sequence)
            new_file.close()
            return(new_file_name)




def compare_nucleotides(base_1, base_2):
        ''' Compare two nucleotides
        '''
        if base_1 != base_2:
            return("the nucleotides given are DIFFERENT",2)

        elif base_1== base_2:
            return("the nucleotides given are the SAME",1)
        
        else:
            return "not found"

