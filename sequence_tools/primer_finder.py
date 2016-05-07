import re, os.path, sys, click
from useful_tools import useful
from useful_tools.process_file import ProcessIO

# gives the dir/PATH of this program, will be used to find the defaulted primer_database
# which are hg19 primers only.
file_path = useful.cwd_file_path(__file__) 

@click.command('primer_finder') 
@click.argument('input_file',nargs=1, required=False)
@click.option('--output_file',default=None,help="output; defaulted as matching_primers_output.txt")
@click.option('--primer_database',default=file_path+"useful_tools/TAAD_Primer_Validation_Database.txt",help="defaulted to TAAD primer DB")  
@click.option('--delimiters', default="\t",help="defaulted to tab")

def matching_primer(input_file,
                    output_file = None,
                    primer_database = file_path+"TAAD_Primer_Validation_Database.txt",
                    delimiters = "\t"):
                    
    ''' Takes variant postion(s) as input and matches it with an appropriate primer
        pair in a given primer database. 
        \b\n
        input can be either: FILE; tab deliminated file with variant name and variant
        position, STRING;  variant position
        \b\n
        Examples:\b\n
            primer_finder 15:48729400 --primer_database primer_database.txt
            primer_finder --input_file in.txt --output_file out.txt 
        
    '''
    try:
        # get all genomic locations within primer pairs, from all primers in the database
        all_primer_pos = get_all_primer_pos(primer_database)
    except IOError as e:
        print(primer_database+": "+e.strerror)
        sys.exit(0)
        
    try:

        process_io = ProcessIO(input_file,output_file)
        input_file = process_io.process_input()
        for info in input_file:
            var_name = info.split("\t")[0]
            var_pos = info.split("\t")[1].rstrip("\n")
            matched_primers = match(var_pos,all_primer_pos,var_name)
            print(matched_primers)
   
    except IndexError as e:
        return e.args[0]+"\n\n"+input_file+\
        " is incorrectly deliminated or an incorrect delimiter was specified."+"\n"    
            
def get_all_primer_pos(primer_database):
    ''' Generate a list of every genomic position witin each primer pair given 
        in the primer database.
    '''
    try:
        primer_file = open(primer_database,"r")
        header = next(primer_file)    # skip database header
        
        # list of lists, where every single list is every genomic position within a primer pair
        # followed by bp distance of genomic position from the forward and reverse primer position
        all_primer_pos = []            
        
        # generates all_primer_pos
        for i in primer_file:
            
            # split into chromosome, start position and stop position in which primer occupies
            i = i.split("\t")
            primer_name = i[0]
            primer_range = re.split(r'[:-]',i[4])   # split into three components
            chrom = re.sub(r'[^0-9]','',primer_range[0])  # remove any non-integers
            start = int(primer_range[1])
            stop = int(primer_range[2])
            
            # use start and stop to iterate through every position in primer pair and
            # output primer_name, primer positions, distance of primer position from the
            # forward and reverse primer respectively and append to an empty list.
            all_pos_in_primer = [ " ".join((primer_name,chrom+":"+str(i),
                                            str(i-start),str(stop-i))) 
                                for i in range(start,stop)]
            all_primer_pos.append(all_pos_in_primer)
        
        # unpacks list of lists into a single list
        all_primers_pos_unpacked = [x for i in all_primer_pos for x in i] 
        return all_primers_pos_unpacked
        
    except IndexError:
        raise 
        
        
def match(var_pos,primer_info,var_name=None):
    ''' Match a given variant position against every genomic position covered
        in all the primer pairs in the primer database.
    '''
    # used if a string given as input instead of a file
    if var_name is None:
        var_name = var_pos
    
    # stored matched primer information
    answer = []
    
    # searches and generates matched primer pair for the variant position given
    for i in primer_info:   
        primer_name = i.split(" ")[0]
        primer_pos = i.split(" ")[1]
        variant_distance_f = i.split(" ")[2]
        variant_distance_r = i.split(" ")[3]    
        if var_pos == primer_pos:
            match = "\t".join((var_name,primer_name, variant_distance_f,
                               variant_distance_r,"\n"))
            answer.append(match)  
    
    # returns no match error if no primer pair is found, else return answer as string
    if not answer:
        return  "no match found for: "+var_name
    else:
        return "".join(answer)
                 
        
        

# required for click to be useable, dont know why
if __name__ == '__main__':
    matching_primer()
    
