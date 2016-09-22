from __future__ import division, print_function
import os, click, sys, traceback
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
if os.path.isfile(home+"/.config/GeneaPy.config") is True:
    ttuner_path, default_hg_version, genome_path = config.read_config("hg19")
else:
    ttuner_path, default_hg_version, genome_path = (None, None, None)


@click.command('main')
@click.argument('input_file',nargs=1, required=False)
@click.option('--output_file',default=None, help='give an output file, requires input to be a file')
@click.option('--upstream', default=20, help="number of bases to get upstream, default: 20") 
@click.option('--downstream',default=20, help="number of bases to get downstream, default: 20")
@click.option('--hg_version',default=default_hg_version, help="human genome version. default: hg19")
@click.option('--header/--no_header',default='n',help="header gives metadata i.e. sequence name etc.")
@click.option('--seq_file', default=None, help="match sequence with AB1/seq file contents")
@click.option('--seq_dir', default=file_path[:-8]+"test/test_files/", help="path to directory containing AB1/seq files")
@click.option('--genome', default=genome_path, help="use locally stored genome instead of UCSC")
@click.option('--ensembl/--ne', default='n', help="scrape gene & exon information")
@click.option('--download/--donot', default='n', help="download ttuner and/or human genome fasta files")
@click.option('--configure/--noconfig', default='n', help="configure get_seq")
@click.option('--vcf', default=None, help='convert vcf to input file')

def main(input_file, output_file=None, upstream=20, downstream=20,
         hg_version=default_hg_version, header="n", seq_file=None, download='n', 
         configure='n',seq_dir=file_path[:-8]+"test/test_files/", ensembl="n", 
         genome=genome_path, vcf=None):
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
        vcf2input.vcf2input(vcf, "temp.tsv", 4)
        input_file = "temp.tsv"

    # ensures the wrong human genome isnt used
    if hg_version not in genome:
        genome = None

    # allows one to pipe in an argument at the cmd
    input_file = input() if not input_file else input_file

    # intiate the classes in UCSC 
    reference = ScrapeSeq(input_file, upstream, downstream, hg_version, ensembl, header)
    
    # if the arg given is a file, parse it in line by line otherwise assume its a string
    if os.path.isfile(input_file) is True:
        all_scrapped_info = parse_file(input_file, output_file, upstream, downstream,
                                       hg_version, header, seq_file, seq_dir, reference,  ensembl, genome)
    else:
        parse_string(input_file, output_file, upstream, downstream, hg_version, header,
                     seq_file, seq_dir, reference, ensembl, genome)
    
    # delete the converted tsv file after usage
    if vcf:
        os.remove("temp.tsv")

    # write the header and each element per line to the file
    if output_file:
        header = "\t".join(("Name", "Position", "Seq_Range", "Gene_Name", 
                            "Gene_ID", "Type", "Gene_Range", "Transcript","Exon_ID",
                            "Intron", "Exon", "AB1", "Ref", "Seq", "Result", "\n"))
        write_to_output(all_scrapped_info, output_file, header)




def parse_file(*args):    
    ''' Parse a file line by line into the get_seq function
    '''
    input_file, output_file, upstream, downstream, hg_version, header, \
           seq_file, seq_dir, reference,  \
        ensembl, genome = args
    
    # print a warning if --seq_file is used with input as a file
    print('WARNING: --seq_file argument ignored. Automatically searching seq_files'
          ' directory for a matching file\n' if seq_file else '')
    
    # append all returned data to this list
    all_scrapped_info = []

    for line in [line.rstrip("\n").split("\t") for line in open(input_file)]:
        seq_name = line[0]
        var_pos = line[1]
        
        # Determine the type of mutation
        if len(line) > 2:
            mutation = line[2]
            ref_base = line[2].split("/")[0]
            print(line[2].split("/"))
            alt_answer, mut_type = determine_mutation(mutation)
        else:
            mut_type = None
            alt_answer = None

        # check each individual line of the file for CUSTOM ERRORS
        error_check = reference.handle_argument_exception(var_pos)
        # intialise the class in Ensembl.py
        pyensembl = ScrapeEnsembl(var_pos, hg_version) if ensembl else None
        # find a list of files with seq_name in its title, if ab1 matched then convert it to a .seq and .tab file
        seq_file = get_AB1_file.get_matching_seq_file(seq_name, seq_dir)
        convert = [get_AB1_file.handle_seq_file(x, seq_dir) 
                   for x in seq_file if seq_file]
        # intialise the check_sanger class for every found seq_file
        sanger = [CompareSeqs(upstream, downstream, alt_answer, mut_type, x, seq_dir) 
                  for x in seq_file]
        # parse it all into get_seq()
        sequence_info = [get_seq(seq_name, var_pos, reference,  
                                hg_version, pyensembl, genome, x) for x in sanger]
        # filter out sequences where no seq file was found
        filtered_answer = [x for x in sequence_info if "-" != x[1].split("\t")[11]]
        
        # find the index in filtered_answer which contains the mutation detected by NGS and store in a var called found_answer
        index = 0
        for i in filtered_answer:
            filtered_call = i[1].split("\t")[13]

            if ref_base+"/"+alt_answer == filtered_call:
                found_answer = filtered_answer[index]  

            index += 1

        # if a het call was found in some of the matching seq files, then print and return its values. Otherwise, return no seq_file matched values
        if found_answer:
            print_out, answer = found_answer
            print(print_out)
            all_scrapped_info.append(answer)

        elif filtered_answer:
            print_out, answer = filtered_answer[0]
            print(print_out)
            # append each lines returned data to an emty list
            all_scrapped_info.append(answer)

        elif sequence_info:
            print_out, answer = sequence_info[0]
            print(print_out)
            all_scrapped_info.append(answer)

        # reset variable for next line
        found_answer = None

    # return all scrapped data
    return all_scrapped_info



def determine_mutation(mutation):
    ''' Determine which type of mutation is in each line of the
        input file

        This will eventually be utilised in the check_sangers het calls to find the correct mutation. 
    '''
    ref, alt = tuple(mutation.split("/"))

    # determine whether mut is a SNP, deletion or insertion
    if len(ref) == len(alt):
        print("We are looking for a SNP in {}")#.format(name))
        return (alt, "snp")
    elif len(ref) > len(alt):
        print("We are looking for a DELETION in {}")#.format(name))
        return (alt, "d")
    elif len(ref) < len(alt):
        print("we are looking for an INSERTION in {}")#.format(name))
        return (alt, "i")





def parse_string(*args):
    ''' Parse a string into the get_seq function
    '''
    input_file, output_file, upstream, downstream, hg_version, header, \
           seq_file, seq_dir, reference,  \
        ensembl, genome = args
 
    # check the input for CUSTOM ERRORS and intialise the Ensembl class
    error_check = reference.handle_argument_exception(input_file)
    pyensembl = ScrapeEnsembl(input_file, hg_version) if ensembl else None
    
    # intialiase the check_sanger class and parse it all into the get_seq()
    #if seq_file and seq_file.endswith("ab1"):
    #    seq_file = CompareSeqs.convert_ab1_to_seq(seq_file)
    sanger = CompareSeqs(upstream,downstream, seq_file, seq_dir) if seq_file else None
    sequence_info = get_seq("query", input_file, reference,  
                             hg_version, pyensembl, genome, sanger)

    # this condition fixes a subscripting error that occurs if ensembl_error() is triggered
    if isinstance(sequence_info, tuple):
        print(sequence_info[0])
    else:
        print(sequence_info)
    

    

def get_seq(seq_name, var_pos, reference,  hg_version, pyensembl, genome, sanger=None):
    ''' Get the reference sequence from a given position, possibly compare
        to a sanger squence and get gene/exon information. 

        sanger and gene/exon information is dependent upon whether the sanger 
        and pyensembl is or is not None
    '''
    try:
        # default variables to - incase the below conditions are not met
        gene_name = gene_id = gene_type = gene_range = transcript = exon_id = \
            exon_id = intron = exon = ref_base = sanger_base = seq_file_used = "-"
        compare_result = "0"

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
            # if .seq file found compare the ref var_pos base and the sanger var_pos base
            if sanger_sequence:
                seq_file_used = sanger.seq_filename.split("/")[-1]
                upstream_seq, downstream_seq, ref_base, sanger_base, \
                        var_index = sanger_sequence
                
                # if var_index is a tuple then it is an insertion else process the position
                if isinstance(var_index, tuple):
                    alternate_bases = [sanger.seq_file[x] for x in range(var_index[0], var_index[1]+1)]
                    insertion = "".join(alternate_bases)
                    if var_index[-1] == "i":
                        het_call = "/".join((ref_base,insertion))
                    if var_index[-1] == "d":
                        het_call = "/".join((insertion, ref_base))

                else:
                    alternate_bases = sanger.get_het_call(var_index)
                    het_call = sanger.base_caller(alternate_bases, ref_base)
                    

                # if a het was found use it in the full sequence else use the ref base
                if het_call:
                    sanger_base = het_call 
                    full_seq = "".join((upstream_seq,het_call,downstream_seq))
                else:
                    full_seq = "".join((upstream_seq,sanger_base,downstream_seq))

                compare = get_AB1_file.compare_nucleotides(ref_base,sanger_base) 

                statement = compare[0]
                compare_result = compare[1]

        # get gene information for the variant position
        if pyensembl:
            all_info = get_gene_exon_info.gene_transcript_exon(var_pos, hg_version)
            gene_name, gene_id, gene_type, gene_range = all_info[0]
            transcript = all_info[1]
            exon_id, intron, exon = all_info[2]
             
        # determine whether to give a HEADER
        header = reference.header_option(seq_name,var_pos,seq_range,sequence, gene_name)

        # print reference sequence (no options)
        if 'sanger_sequence' not in locals():
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
        all_data = (seq_name, var_pos, seq_range, gene_name,gene_id, gene_type,
                    gene_range, transcript, exon_id, intron, exon, seq_file_used, ref_base, sanger_base, str(compare_result))
        answer = "\t".join(all_data)
        
        return (print_out, answer)
   

    except ValueError as e:
        t, v, tb = sys.exc_info()
        sam = traceback.format_exception(t, v, tb)
        msg = " ".join(("ERROR: No %s information found for",
                        seq_name,"at",var_pos,"in",hg_version))
        ensembl_error(sam[1], msg)
        exon_id = intron = exon = "-"



def ensembl_error(error, msg):
    ''' Handle ValueErrors that occur within get_seq()
        and determine whether they affect gene, transcript 
        or exona information. Print to the screen and write 
        to an error file.
    '''   
    with open(file_path + "stderror_get_seq.txt", "a") as out:
        if "gene" in error:
            msg = msg % "gene"
        elif "transcript" in error:
            msg = msg % "transcript"
        elif "exon" in error:
            msg = msg % "exon"
        else:
            msg = "unknown cause of ValueError"
            
        print(msg)
        out.write(msg+"\n")

               
if __name__ == '__main__':
    main()         

# main("in.txt","b", 2, 20,"hg19","\t","y")

