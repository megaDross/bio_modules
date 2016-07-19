from __future__ import division, print_function
import os, sys,re, click, requests, bs4, json
import warnings
from Ensembl import ScrapeEnsembl, get_exon_number
from check_sanger import CompareSeqs 
from transcription_translation import ProteinRNA
from useful_tools.output import write_to_output
from useful_tools import useful
from UCSC import ScrapeSeq

# get the absolute path to this file
file_path = useful.cwd_file_path(__file__)
 
@click.command('main')
@click.argument('input_file',nargs=1, required=False)
@click.option('--output_file',default=None, help='give an output file, requires input to be a file')
@click.option('--upstream', default=20, help="number of bases to get upstream, default: 20") # default to an int makes option accept int only
@click.option('--downstream',default=20, help="number of bases to get downstream, default: 20")
@click.option('--hg_version',default="hg19", help="human genome version. default: hg19")
@click.option('--header/--no_header',default='n',help="header gives metadata i.e. sequence name etc.")
@click.option('--transcribe/--nr',default='n',help="transcribe into RNA sequence")
@click.option('--translate/--np',default='n',help="translate RNA seq into protein seq")
@click.option('--rc/--no_rc',default='n',help="reverse complement the DNA")
@click.option('--seq_file', default=None, help="match sequence with .seq file contents")
@click.option('--seq_dir', default=file_path+"seq_files/")
@click.option('--ensembl/--ne', default='n', help="scrape gene & exon information")

def main(input_file, output_file=None, upstream=20, downstream=20, hg_version="hg19",
         header="n", transcribe="n", translate="n", rc='n', seq_file=None,
         seq_dir=file_path+"seq_files/", ensembl="n"):
    '''
    From a genomic postion, genomic range or tab-deliminated file produce a
    reference sequence that can be compared with a sanger trace along with 
    said position/ranges gene/transcript/exon information. 
         \b\n
    A file or string can be used as input. STRING: either a variant position 
    or a genomic range deliminated by a comma. FILE: deliminated file with 
    the variant name and the variant position
         \b\n
    Example:\b\n
        get_seq chr1:169314424 --upstream 200 --downstream 200\n
        get_seq chr1:169314424,169314600 --hg_version hg38\n
        get_seq 1:169314424 --seq_file dir/test.seq 
        get_seq input.txt --output_file output.txt --header --seq_dir dir/\n
    ''' 
    # allows one to pipe in an argument at the cmd
    input_file = input() if not input_file else input_file
    
    # intiate the classes in UCSC and transcribe_translate respectively
    reference = ScrapeSeq(input_file,  upstream, downstream, hg_version, ensembl, header)
    trans = ProteinRNA(transcribe, translate, rc)
    
    # if the arg given is a file, parse it in line by line otherwise assume its a string
    if os.path.isfile(input_file) is True:
        all_scrapped_info = parse_file(input_file, output_file, upstream, downstream,
                                      hg_version, header, transcribe, translate, rc,
                                      seq_file, seq_dir, reference, trans, ensembl)
    else:
        parse_string(input_file, output_file, upstream, downstream, hg_version, header,
                     transcribe, translate, rc, seq_file, seq_dir, reference, trans,
                     ensembl)
    
    # write the header and each element per line to the file
    if output_file:
        header = "\t".join(("Name", "Position", "Seq_Range", "Gene_Name", 
                            "Gene_ID", "Type", "Gene_Range", "Transcript","Exon_ID",
                            "Intron", "Exon", "Ref", "Seq", "Result", "\n"))
        write_to_output(all_scrapped_info, output_file, header)




def parse_file(*args):    
    ''' Parse a file line by line into the get_seq function
    '''
    input_file, output_file, upstream, downstream, hg_version, header, \
    transcribe, translate, rc, seq_file, seq_dir, reference, trans, ensembl = args
    
    # print a warning if --seq_file is used with input as a file
    print('WARNING: --seq_file argument ignored. Automatically searching seq_files'
          ' directory for a matching file\n' if seq_file else '')
    
    # append all returned data to this list
    all_scrapped_info = []

    for line in [line.rstrip("\n").split("\t") for line in open(input_file)]:
        seq_name = line[0]
        var_pos = line[1]
        # check each individual line of the file for CUSTOM ERRORS
        error_check = reference.handle_argument_exception(var_pos)
        # intialise the class in Ensembl.py
        pyensembl = ScrapeEnsembl(var_pos, hg_version) if ensembl else None
        # find a .seq file with seq_name in its title
        seq_file = CompareSeqs.get_matching_seq_file(seq_name, seq_dir)
        # intialise the check_sanger class
        sanger = CompareSeqs(upstream, downstream, seq_file)
        # parse it all into get_seq()
        sequence_info = get_seq(seq_name, var_pos, reference, trans, 
                                hg_version, pyensembl, sanger)
        # append each lines returned data to an emty list
        all_scrapped_info.append(sequence_info)
    # return all scrapped data
    return all_scrapped_info




def parse_string(*args):
    ''' Parse a string into the get_seq function
    '''
    input_file, output_file, upstream, downstream, hg_version, header, \
    transcribe, translate, rc, seq_file, seq_dir, reference, trans, ensembl = args
 
    # check the input for CUSTOM ERRORS and intialise the Ensembl class
    error_check = reference.handle_argument_exception(input_file)
    pyensembl = ScrapeEnsembl(input_file, hg_version) if ensembl else None
    
    # intialiase the check_sanger class and parse it all into the get_seq()
    sanger = CompareSeqs(upstream,downstream, seq_file) if seq_file else None
    get_seq("query", input_file, reference, trans, hg_version, pyensembl, sanger)

    

    

def get_seq(seq_name, var_pos, reference, trans, hg_version, pyensembl, sanger=None):
    # default variables incase the below conditions are not met
    gene_name = gene_id = gene_type = gene_range = transcript = exon_id = \
    exon_id = intron = exon = ref_base = sanger_base = "-"
    compare_result = "0"

    # check if var_pos is a GENOMIC REGION, else construct one from var_pos
    seq_range = reference.create_region(var_pos)
    
    # use UCSC to get the genomic ranges DNA sequence
    sequence = reference.get_region_info(seq_range)
    
    # if CompareSeqs class has been intiated try and find a matching .seq file
    if sanger:
        sanger_sequence = sanger.match_with_seq_file(sequence)
        # if .seq file found compare the ref var_pos base and the sanger var_pos base
        if sanger_sequence:
            ref_base = sanger_sequence[1] 
            sanger_base = sanger_sequence[2] 
            compare = CompareSeqs.compare_nucleotides(ref_base,sanger_base) 
            statement = compare[0]
            compare_result = compare[1]


    # determine whether to transcribe or translate to rna or protein EXPERIMENTAL
    sequence = str(trans.get_rna_seq(sequence))
    sequence = str(trans.get_protein_seq(sequence))
     
    # get gene information for the variant position
    if pyensembl:
        gene_info = pyensembl.get_gene_info()
        if gene_info:
            gene_name, gene_id, gene_type, gene_range = gene_info
            transcript = pyensembl.get_canonical_transcript(gene_name)

            if transcript:
                exon_info = get_exon_number(transcript, hg_version, var_pos)
                exon_id, intron, exon = exon_info
            

    # determine whether to give a HEADER
    if not pyensembl:
        header = reference.header_option(seq_name,var_pos,seq_range,sequence)
    else:
        header = reference.header_option(seq_name,var_pos,seq_range,sequence, gene_name)




    # print reference sequence (no options)
    if 'sanger_sequence' not in locals():
        print("\n".join((header, "Reference Sequence:\t"+sequence,"\n")))

    # print above and sanger sequence
    elif sanger_sequence :
        print("\n".join((header,"Reference Sequence:\t"+sequence,
                         "Sanger Sequence:\t"+sanger_sequence[0],
                         statement,"\n")))
    
    # if no matching seq file 
    elif not sanger_sequence:
        print("\n".join((header,"Reference Sequence:\t"+sequence,
                         "Sanger Sequence:\tNo Match Found", "\n")))


    # return everything in a tab delimited format 
    return("\t".join((seq_name, var_pos, seq_range, gene_name,
                      gene_id, gene_type, gene_range, transcript, 
                      exon_id, intron, exon, ref_base,
                      sanger_base, str(compare_result))))

 

               
if __name__ == '__main__':
    main()         

# main("in.txt","b", 2, 20,"hg19","\t","y")

