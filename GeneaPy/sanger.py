from __future__ import division, print_function
import os, click, sys, traceback, logging, datetime
import get_gene_exon_info
import useful
from Ensembl import ScrapeEnsembl
from check_sanger import CompareSeqs 
from output import write_to_output
from UCSC import ScrapeSeq
import subprocess
import config
import get_AB1_file
import vcf2input

# get the absolute path to this file
file_path = useful.cwd_file_path(__file__)
home = os.path.expanduser("~")

# read the config file
if os.path.isfile(home+"/.config/GeneaPy/GeneaPy.config") is True:
    ttuner_path, default_hg_version, genome_path = config.read_config("hg19")
else:
    ttuner_path, default_hg_version, genome_path = (None, None, None)

# create a basic logging file
logging.basicConfig(filename=os.path.expanduser('~')+"/.config/GeneaPy/GeneaPy.log", level=logging.DEBUG)

logging.info("Itialisation at {}".format(datetime.datetime.now().strftime("%d/%m/%Y -  %I:%M%p")))

@click.command('main')
@click.argument('input_file',nargs=1, required=False)
@click.option('--output_file',default=None, help='give an output file, requires input to be a file')
@click.option('--upstream', default=20, help="number of bases to get upstream, default: 20") 
@click.option('--downstream',default=20, help="number of bases to get downstream, default: 20")
@click.option('--hg_version',default=default_hg_version, help="human genome version. default: hg19")
@click.option('--seq_dir', default=file_path[:-8]+"test/test_files/", help="path to directory containing AB1/seq files")
@click.option('--genome', default=genome_path, help="use locally stored genome instead of UCSC")
@click.option('--download/--donot', default='n', help="download ttuner and/or human genome fasta files")
@click.option('--configure/--noconfig', default='n', help="configure get_seq")
@click.option('--vcf', default=None, help='convert vcf to input file')

def main(input_file, output_file=None, upstream=20, downstream=20,
         hg_version=default_hg_version,  download='n', configure='n',seq_dir=file_path[:-8]+"test/test_files/", genome=genome_path, vcf=None):
    '''
    From a genomic postion, genomic range or tab-deliminated file produce a
    reference sequence that can be compared with a sanger trace along with 
    said position/ranges gene/transcript/exon information. 
    ''' 
    # whether to download ttuner and or human genomes
    if download:
        subprocess.call([file_path+"download.sh"])
        config.config() 
        sys.exit()

    # write a config file
    if configure:
        config.config()
        sys.exit()
   
    # convert vcf to input file
    if vcf:
        vcf2input.vcf2input(vcf, "temp.tsv")
        sys.exit(0)

    # ensures the wrong human genome isnt used
    if hg_version not in genome:
        genome = None


    # intiate the classes in UCSC 
    reference = ScrapeSeq(input_file, upstream, downstream, hg_version)
    
    # if the arg given is a file, parse it in line by line otherwise assume its a string
    all_scrapped_info = parse_file(input_file, output_file, upstream, downstream,
                                   hg_version, seq_dir, reference, genome)
    

    # write the header and each element per line to the file
    if output_file:
        header = "\t".join(("Name", "Position", "Seq_Range", "AB1", "NGS", "Ref", "Seq",
                            "Result", "\n"))
        write_to_output(all_scrapped_info, output_file, header)



def parse_file(*args):    
    ''' Parse a file line by line into the get_seq function
    '''
    input_file, output_file, upstream, downstream, hg_version, \
            seq_dir, reference, genome = args
    
    # append all returned data to this list
    all_scrapped_info = []

    for line in [line.rstrip("\n").split("\t") for line in open(input_file)]:
        seq_name = line[0]
        var_pos = line[1]
        
        logging.info("\n\n\t We are looking for {} in sample {}\n\n".format( var_pos, seq_name))

        # Determine the type of mutation
        mutation = line[2]
        ref_base = line[2].split("/")[0]
        alt_answer, mut_type = determine_mutation(mutation)

        logging.debug(" We are looking for a {}\n".format(mutation))

        # check each individual line of the file for CUSTOM ERRORS
        error_check = reference.handle_argument_exception(var_pos)

        # find a list of files with seq_name in its title.
        seq_file = get_AB1_file.get_matching_seq_file(seq_name, seq_dir)

        # produce seq & tab files for each matched file, otherwise return blank line
        if seq_file:
            convert = [get_AB1_file.handle_seq_file(x, seq_dir) 
                       for x in seq_file]
        else:
            all_scrapped_info.append("\t".join((seq_name, var_pos, "-", "-", "-", "-", "-", "0")))    
        
        # intialise the check_sanger CompareSeq class for every found seq_file
        sanger = [CompareSeqs(upstream, downstream, ref_base, alt_answer, 
                              mut_type, x, seq_dir) 
                  for x in seq_file]

        # parse it all into get_seq()
        sequence_info = [get_seq(seq_name, var_pos, mutation, reference,  
                                hg_version, genome, x) for x in sanger]
        
        logging.debug("Sequence Info: ".format(sequence_info))

        # filter out sequences where no seq file was found
        filtered_answers = [x for x in sequence_info if "-" != x[1].split("\t")[4]]
        logging.debug("Filtered Answer: {}".format(filtered_answers))
       
        # check if any of the filtered_answer list contains the het call of interest
        found_answer = find_variant(filtered_answers, ref_base, alt_answer)

        # if a het call was found in some of the matching seq files, then print and return its values. Otherwise, return no seq_file matched values
        if found_answer:
            print_out, answer = found_answer
            print(print_out)
            logging.debug("printing: {}".format(print_out))
            all_scrapped_info.append(answer)

        elif filtered_answers:
            print_out, answer = filtered_answers[0]
            print(print_out)
            logging.debug("printing: {}".format(print_out)) 
            # append each lines returned data to an emty list
            all_scrapped_info.append(answer)

        elif sequence_info:
            print_out, answer = sequence_info[0]
            print(print_out)
            logging.debug("printing: {}".format(print_out)) 
            all_scrapped_info.append(answer)

        # reset variable for next line/iteration
        found_answer = None

    # return all scrapped data
    return all_scrapped_info



def determine_mutation(mutation):
    ''' Determine which type of mutation is in each line of the
        input file. Assumes the mutation is a het.

        This will eventually be utilised in the check_sangers het calls to find the correct mutation. 
    '''
    ref, alt = tuple(mutation.split("/"))

    # determine whether mut is a SNP, deletion or insertion
    if len(ref) == len(alt):
        logging.info("We are looking for a SNP")#.format(name))
        return (alt, "snp")
    elif len(ref) > len(alt):
        logging.info("We are looking for a DELETION")#.format(name))
        return (alt, "d")
    elif len(ref) < len(alt):
        logging.info("we are looking for an INSERTION")#.format(name))
        return (alt, "i")



def find_variant(filtered_answers, ref_base, alt_answer):
    ''' From a list of get_seq() outputs, check that a given
        het variant called vai NGS is within one of the lists 
        elements and, if found, return said element.
    '''
    # required in case filtered_answers is an empty list
    found_answer = None

    index = 0
    for answer in filtered_answers:
        filtered_call = answer[1].split("\t")[6]
        logging.debug("Filtered Call: {}, NGS Call: {}".format(filtered_call, 
                                                          ref_base+"/"+alt_answer))
        if ref_base+"/"+alt_answer == filtered_call:
            found_answer = filtered_answers[index]  
            break
        else:
            found_answer = None
        index += 1

    return found_answer



def get_seq(seq_name, var_pos, mutation, reference,  hg_version,  genome, sanger=None):
    ''' Get the reference sequence from a given position, possibly compare
        to a sanger squence and get gene/exon information. 

        sanger and gene/exon information is dependent upon whether the sanger 
        and pyensembl is or is not None
    '''
    # default variables to - incase the below conditions are not met
    ref_base = sanger_base = seq_file_used = "-"
    compare_result = "0"
    
    # get the seq file used for the analysis
    seq_file_used = sanger.seq_filename.split("/")[-1]
    
    # check if var_pos is a GENOMIC REGION, else construct one from var_pos
    seq_range = reference.create_region(var_pos)
    
    # use loacally stored genome or UCSC to get the genomic ranges DNA sequence
    if genome:
        sequence = reference.get_region_info_locally(seq_range, genome)
    else:
        sequence = reference.get_region_info(seq_range)
    
    # if CompareSeqs class has been intiated try and find a matching .seq file
    if sanger:
        
        sanger_sequence = sanger.match_with_seq_file(sequence)
        logging.debug("{}\n".format(sanger_sequence))

        # if .seq file found compare the ref var_pos base and the sanger var_pos base
        if sanger_sequence:
            upstream_seq, downstream_seq, ref_base, sanger_base, \
                poly_het = sanger_sequence

            # decipher het call
            het_call = determine_het_call(sanger_sequence) 
               
            # if a het was found use it in the full sequence else use the ref base
            if het_call:
                sanger_base = het_call 
                full_seq = "".join((upstream_seq,het_call,downstream_seq))
            else:
                full_seq = "".join((upstream_seq,sanger_base,downstream_seq))

            statement, compare_result = get_AB1_file.compare_nucleotides(ref_base,sanger_base) 


    # construct a header
    header = reference.header_option(seq_name,var_pos,seq_range,sequence)

    # print reference sequence
    if not sanger:
        print_out ="\n".join((header, sequence,"\n"))

    # print reference and sanger sequence
    elif sanger_sequence :
        print_out = "\n".join((header,"Reference Sequence:\t"+sequence,
                               "Sanger Sequence:\t"+full_seq,
                               statement,"\n"))

    # if no matching seq file 
    elif not sanger_sequence:
        print_out = "\n".join((header,"Reference Sequence:\t"+sequence,
                               "Sanger Sequence:\tNo Match Found", "\n"))

    # print results and return everything for outputing to a file
    all_data = (seq_name, var_pos, seq_range, seq_file_used, mutation, ref_base,
                sanger_base, str(compare_result))
    answer = "\t".join(all_data)
    
    return (print_out, answer)
   


def determine_het_call(sanger_sequence):
    ''' Determine het call from the output of
        the CompareSeqs.match_with_seq_file()
    '''
    upstream_seq, downstream_seq, ref_base, sanger_base, \
                    poly_het = sanger_sequence

    # if poly_het is a tuple then it is an indel else process the position
    if isinstance(poly_het, tuple):
        alternate_bases = [sanger.seq_file[x] 
                           for x in range(poly_het[0], poly_het[1]+1)]
        indel = "".join(alternate_bases)
        if poly_het[-1] == "i":
            het_call = "/".join((ref_base,indel))
        if poly_het[-1] == "d":
            het_call = "/".join((indel, ref_base))
    
    # string means it is a deletion which was found in the het call dict
    elif isinstance(poly_het, str):
        het_call = poly_het

    else:
        het_call = None
              
    return het_call
               

if __name__ == '__main__':
    main()         


