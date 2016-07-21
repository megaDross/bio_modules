import re, os.path, click
import useful
from output import write_to_output

# gives the dir/PATH of this program, will be used to find the defaulted primer_database
# which are hg19 primers only.
file_path = useful.cwd_file_path(__file__) 


@click.command('primer_finder') 
@click.argument('input_file',nargs=1, required=False)
@click.option('--output_file',default=None,help="output; defaulted as matching_primers_output.txt")
@click.option('--primer_database',default=file_path+"TAAD_Primer_Validation_Database.txt",help="defaulted to TAAD primer DB")  
@click.option('--delimiters', default="\t",help="defaulted to tab")
def main(input_file, output_file = None,
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
        # allows one to pipe in an argument at the cmd, requires required=False in 
        # @click.argument()
        if not input_file:
            input_file = input()
       
       # get all genomic locations within primer pairs, from all primers in the database
        all_primer_pos = get_all_primer_pos(primer_database)        
        header = "\t".join(("Variant","Primer","Dist_from_F","Dist_from_R","\n"))

        # determine input type and process accordingly
        if os.path.isfile(input_file) is True:
            all_matched_primers = []
            for line in [line.rstrip("\n").split("\t") for line in open(input_file)]:
                var_name = line[0]
                var_pos = line[1].replace("chr","")
                matched_primers = match(var_pos,all_primer_pos,var_name)
                all_matched_primers.append(matched_primers)
                print(matched_primers)

        else:
            position = input_file.replace("chr","")
            matched_primers = match(position,all_primer_pos,"query")
            print(header[:-1])
            print(matched_primers)

        
        if output_file:
            write_to_output(all_matched_primers, output_file, header)

    except IndexError as e:
        return e.args[0]+"\n\n"+input_file+\
        " is incorrectly deliminated or an incorrect delimiter was specified."+"\n"    





def get_all_primer_pos(primer_database):
    ''' Generate a list of every genomic position witin each primer pair given 
        in the primer database.
    '''
    try:
        primer_file = open(primer_database)
        header = next(primer_file)    # skip database header
        
        # list of lists, where every single list is every genomic position within a primer pair
        # followed by bp distance of genomic position from the forward and reverse primer position
        all_primer_pos = []     

        # generates all_primer_pos
        for i in [i.rstrip("\n").split("\t") for i in primer_file]:
            primer_name = i[0]
            primer_range = re.split(r'[:-]',i[4])   # split into three components
            chrom = re.sub(r'[^0-9]','',primer_range[0])  # remove any non-integers
            start = int(primer_range[1])
            stop = int(primer_range[2])
            
            # use start and stop to iterate through every position in primer pair and
            # output primer_name, primer positions, distance of primer position from the
            # forward and reverse primer respectively and append to an empty list.
            all_pos_in_primer = tuple(" ".join((primer_name,chrom+":"+str(i),
                                                str(i-start),str(stop-i))) 
                                                for i in range(start,stop))
            
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
                               variant_distance_r))
            answer.append(match)  
    
    # returns no match error if no primer pair is found, else return answer as string
    if not answer:
        return  "no match found for: "+var_name
    else:
        return "".join(answer)
                 
    
        

# excute only from outside file
if __name__ == '__main__':
    main()
    
