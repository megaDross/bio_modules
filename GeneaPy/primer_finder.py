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
@click.option('--distance', type=int, help="number of bp a primer must be from the position")
@click.option('--size', type=int, help="maximum desired amplicon size")
@click.option('--gc', type=int, help="maximum desired GC content of amplicon")
def main(input_file, distance, size, gc, output_file = None,
         primer_database = file_path+"TAAD_Primer_Validation_Database.txt"):
                    
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
        header = "\t".join(("Variant","Primer", "Position", "Gene_Name", "Amplicon_Size", "GC%", 
                            "Amplicon_Number", "Dist_from_F","Dist_from_R","\n"))
        print(header[:-1])

        # determine input type and process accordingly
        if os.path.isfile(input_file) is True:
            all_matched_primers = []
            for line in [line.rstrip("\n").split("\t") for line in open(input_file)]:
                var_name = line[0]
                var_pos = line[1].replace("chr","")
                filtered_primer_pos = filter_positions(distance, size, gc,
                                                       all_primer_pos)
                matched_primers = match(var_pos,filtered_primer_pos,var_name)
                all_matched_primers.append(matched_primers)
                print(matched_primers)

        else:
            output_file = None
            position = input_file.replace("chr","")
            filtered_primer_pos = filter_positions(distance, size, gc, all_primer_pos)
            matched_primers = match(position, filtered_primer_pos,"query")
            print(matched_primers)

        
        if output_file:
            write_to_output(all_matched_primers, output_file, header)

    except IndexError as e:
        return e.args[0]+"\n\n"+input_file+\
        " is incorrectly deliminated or an incorrect delimiter was specified."+"\n"    





def get_all_primer_pos(primer_database):
    ''' Generate a list of every genomic position witin each primer pair given 
        in the primer database.

        for every line in the primer database, split the genomic range up into
        chrom, start and stop. Iterate through the start and stop positions to 
        generate a list of every position in the primer and append to an empty
        list. A list of lists of every genomic position in every primer pair
        in a primer database is generated which is subsequently unpacked into
        a single list of tuples. Each tuple contains: 
            ('name', 'position in primer', 'distance from F', 'distance from R')
    '''
    try:
        primer_file = open(primer_database)
        header = next(primer_file)    # skip database header
        
        # list of lists containing every position in every primer pair
        all_primer_pos = []

        # generates all_primer_pos
        for i in [i.rstrip("\n").split("\t") for i in primer_file]:
            primer_name = i[0]
            gene_name = i[3]
            product_size = i[4]
            GC = i[6]
            num_amplicons = i[7]
            primer_range = re.split(r'[:-]',i[5])   # split into three components
            chrom = re.sub(r'[^0-9]','',primer_range[0])  # remove any non-integers
            start = int(primer_range[1])
            stop = int(primer_range[2])      

            # generate every position for current primer pair and append to empty list
            all_pos_in_primer = ["\t".join((primer_name, chrom+":"+str(i), gene_name, 
                                            product_size, GC, num_amplicons,
                                            str(i-start), str(stop-i))) 
                                            for i in range(start,stop)]
            
            all_primer_pos.append(all_pos_in_primer)
            
        # unpacks list of lists into a single list of tuples
        all_primers_pos_unpacked = [tuple(x.split("\t")) for i in all_primer_pos 
                                    for x in i]
        return all_primers_pos_unpacked
        
    except IndexError:
        raise 

       
def filter_positions(distance, size, gc, primer_info):
    ''' filter out primer positions that do not meet the limitations
        set by the user

        Each tuple in the primer_info list should look like:
            (primer name, primer position, product size, gene_name
             GC% of amplicon, number amplicons generated,
             distance of position from F primer, distance from R)
    '''
    # filter out primers that do not generate unique amplicons
    primer_info = [x for x in primer_info
                   if int(x[5]) == 1]

    if distance:
        primer_info = [x for x in primer_info 
                      if int(x[6]) > distance 
                      and int(x[7]) > distance]

    if size:
        primer_info = [x for x in primer_info if 
                       int(x[3].split("bp")[0]) < size]

    if gc:
        primer_info = [x for x in primer_info if
                       int(x[4].split(".")[0]) < gc]
    
    return(primer_info)




def match(var_pos, primer_info, var_name=None):
    ''' Match a given variant position against every genomic position covered
        in all the primer pairs in the primer database
    '''
    # used if a string given as input instead of a file
    if var_name is None:
        var_name = var_pos
    
    # stored matched primer information
    answer = []

    # searches and generates matched primer pair for the variant position given
    for line in primer_info:  
        primer_name, primer_pos, product_size, gene_name, GC, \
            num_amplicons, variant_distance_f, variant_distance_r \
            = line
        
        # append all matching primers to the list
        if var_pos == primer_pos:
            match = "\t".join(((var_name,) + line))
            answer.append(match) 

    # returns no match error if no primer pair is found, else return answer as string
    if not answer:
        return  "\t".join((var_name, "-", "-", "-", "-", "-", "-", "-"))
    else:
        return "\n".join(answer)
                 
    
        

# excute only from outside file
if __name__ == '__main__':
    main()
    
