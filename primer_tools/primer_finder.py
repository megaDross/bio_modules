import re, os.path, sys
import useful 
import click

# gives the dir/PATH of this program, will be used to find the defaulted primer_database
file_path = useful.cwd_file_path(__file__) 


# allows one to use this as a command line tool
# renders the program unuseable in a python terminal
@click.command() #decorator converts the function into a command

@click.option('--input_file', 
              help="input; should be a variant position string or a file")
              
@click.option('--output_file',
              default="matching_primers_output.txt",
              help="output; defaulted as matching_primers_output.txt")
              
@click.option('--primer_database',
              default=file_path+"TAAD_Primer_Validation_Database.txt",
              help="defaulted to TAAD primer DB")
              
@click.option('--delimiters', default="\t",help="defaulted to tab")


def matching_primer(input_file,
                    output_file = "matching_primers_output.txt",
                    primer_database = file_path+"TAAD_Primer_Validation_Database.txt",
                    delimiters = "\t"):
                    
    ''' Takes variant postion(s) as input and matches it with an appropriate primer
        pair in a given primer database. 
        
        Arguments
        -------------------------------------------------------------------------
        input_file -- should be a variant position string or a file in the following 
                      format (tab deliminated):
                            variant_name    variant_position
        
        primer_database -- should be in the following format (tab deliminated):
                              primer name    F-primer    R-primer    product size    genomic region 
        
        delimiters -- defaulted to tab
        
        output_file -- created only if the input given is a file opposed to a string. 
                       If no name is specified for output_file, then it is automatically
                       named matching_primers_output.txt
        
        Returns
        -------------------------------------------------------------------------
        matched_primers -- primer pair(s) that flank the given variant position(s).
                           The distance of the variant from the forward and reverse 
                           primer is given in bases/nucleotides. the tab deliminated 
                           format outputted:
                             variant_name   matched_primer_name distance_f  distance_r
        
        Examples
        --------------------------------------------------------------------------
        from primer_tools import primer_finder 
        
        primer_finder.matching_primer("15:48729400")
        
        primer_finder.matching_primer("variant_positions.txt,"matched_primers.txt")
        
        OR
        
        using click:
        
        python primer_finder --input_file 15:48729400 --primer_database primer_database.txt
        
    '''
    try:
        # get all genomic locations within primer pairs, from all primers in the database
        all_primer_pos = get_all_primer_pos(primer_database)
        
        # if input is not an existing file and contains numbers and a colon 
        if os.path.isfile(input_file) is False and re.search(r"[0-9]:",input_file):             
            var_pos = re.sub(r'[^0-9:]','',input_file)      # remove any characters that isn't a number or colon
            matched_primers = match(var_pos,all_primer_pos)
            click.echo(matched_primers)
            return matched_primers
            
    
        # if input is a file 
        if os.path.isfile(input_file) is True:
            with open(output_file,"w") as output_file:
                output_file.write("\t".join(("Variant_Name","Primer_Name",
                                            "distance_F(bp)","distance_R(bp)",
                                            "\n"))
                                )
                var_file = open(input_file)
                
                for var in var_file:
                    var_name = var.split("\t")[0]
                    var_pos = re.sub(r'[^0-9:]','',var.split("\t")[1])
                    matched_primers = match(var_pos,all_primer_pos,var_name)
                    
                    print matched_primers
                    output_file.write("".join((matched_primers,"\n")))
        
        # give an error message
        else:
            return input_file + " is an invalid input"+"\n\n\t"+\
                                "String or existing file can only be used as input."+"\n\t"+\
                                "String must contain numbers and a colon."
    except IOError as e:
        return e.strerror+"\n"+primer_database+" does not exist."
        
    except IndexError as e:
        return e.args[0]+"\n\n"+ primer_database+" or "+input_file+" is incorrectly deliminated or an incorrect delimiter was specified."+"\n"    
             
                        
            
def get_all_primer_pos(primer_database):
    ''' Generate a list of every genomic position witin each primer pair given 
        in the primer database.
    
        Arguments
        -----------------------------------------------------------------------
        primer_database -- should be in the following format (tab deliminated):
                              primer name    F-primer    R-primer    product size    genomic region 
        
        Returns
        -----------------------------------------------------------------------
        all_primers_pos_unpacked -- the aformentioned list in the following format:
                                        primer_name genomic_position distance_postion_from_F_primer  distance_position_from_R_primer
    '''
    try:
        primer_file = open(primer_database,"r")
        header = primer_file.next()    # skip database header
        
        # list of lists, where every single list is every genomic position within a primer pair
        # followed by bp distance of genomic position from the forward and reverse primer position
        all_primer_pos = []            
        
        # generates all_primer_pos
        for i in primer_file:
            
            i = i.split("\t")
            primer_name = i[0]
            primer_range = re.split(r'[:-]',i[4])   # split into three components
            chrom = re.sub(r'[^0-9]','',primer_range[0])  # remove any non-integers
            start = int(primer_range[1])
            stop = int(primer_range[2])
            
            
            all_pos_in_primer = [ " ".join((primer_name,chrom+":"+str(i),str(i-start),str(stop-i))) 
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
        
        Arguments
        -----------------------------------------------------------------------
        var_pos -- variant position.
        
        primer_info -- the returned list from get_all_primer_pos() function.
        
        var_name -- the variant alias for the variant position. Defaulted to None.
        
        Returns
        ------------------------------------------------------------------------
        answer -- matched primers with distance of F/R primer from variant position.
                  if no match is found then "No match found" is returned. 
        
        
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
            match = "\t".join((var_name,primer_name, variant_distance_f,variant_distance_r))
            answer.append(match)  
    
    # returns no match error if no primer pair is found, else return answer as string
    if not answer:
        return  "no match found for: "+var_name
    else:
        return "".join(answer)
                 
        
        

# required for click to be useable, dont know why
if __name__ == '__main__':
    matching_primer()
    