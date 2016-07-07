from __future__ import division, print_function
import os, sys,re, click, requests, bs4, json
from Ensembl import ScrapeEnsembl, get_exon_number
from check_sanger import CompareSeqs 
from transcription_translation import ProteinRNA
from useful_tools.output import write_to_output
from useful_tools import useful
from UCSC import ScrapeSeq

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

# what use is there in just transcribing and translating for the sake of it? perhaps use an intiation codon finder which can be used o run through a DNA sequence pior to transcripion and transcribing from said site. Also something which checks the dna seq is a multiple of 3 would also be useful. Perhaps something where I could push the rna/protein seq to find which gene exon etc. is or maybe nBLAST, pBLAST etc which I am sure BioPython will have something for parsing to.

def main(input_file, output_file=None, upstream=20, downstream=20, hg_version="hg19",
          header="n", transcribe="n", translate="n",
            rc='n', seq_file=None, seq_dir=file_path+"seq_files/"):
        '''
    Produce a sequence using the UCSC DAS server from an inputted genomic 
	postion and defined number of bases upstream and downstream from said 
	position. A genomic range can be used in place of a genomic position
	and renders the upstream/downstream options irrelevant. An input file
	can have a mixture of genomic positions and genomic ranges. genomic
    ranges is not compatible with the --seq_file option.
         \b\n
    A file or string can be used as input. STRING: either a variant position 
    or a genomic range deliminated by a comma. FILE: deliminated file with 
    the variant name and the variant position
         \b\n
    Example:\b\n
        get_seq chr1:169314424 --dash --upstream 200 --downstream 200\n
        get_seq chr1:169314424,169314600 --hg_version hg38\n
        get_seq input.txt --output_file output.txt --header\n
        ''' 
        # allows one to pipe in an argument at the cmd
        if not input_file:
            input_file = input()
        
        # parse arguments into the ScrapeSeq and ProteinRNA Class 
        reference = ScrapeSeq(input_file,  upstream,
                              downstream, hg_version, header)
        trans = ProteinRNA(transcribe, translate, rc)
        
        # get the path to this file
        
        # if the arg given is a file, parse it in line by line
        if os.path.isfile(input_file) is True:
            all_scrapped_info = []
            if seq_file:
                print("\nWARNING:--seq_file argument ignored. Automatically" 
                      +" searching seq_files directory for a matching file\n")
            for line in [line.rstrip("\n").split("\t") for line in open(input_file)]:
                seq_name = line[0]
                var_pos = line[1]
                ensembl = ScrapeEnsembl(var_pos, hg_version)
                seq_file = CompareSeqs.get_matching_seq_file(seq_name, seq_dir)
                sanger = CompareSeqs(upstream, downstream, seq_file)
                sequence_info = get_seq(seq_name, var_pos, reference, trans, 
                                        sanger, hg_version, ensembl)
                all_scrapped_info.append(sequence_info)

        else:
            ensembl = ScrapeEnsembl(input_file, hg_version)
            sanger = CompareSeqs(upstream,downstream, seq_file)
            get_seq("query", input_file, reference, trans, sanger, hg_version, ensembl)

        if output_file:
            header = "\t".join(("Name", "Position", "Seq_Range", "Gene_Name", 
                                "Gene_ID", "Type", "Gene_Range", "Transcript",
                                "Exon_ID", "Intron", "Exon", "Ref", "Seq", "Result",
                                "\n"))
            write_to_output(all_scrapped_info, output_file, header)




        

def get_seq(seq_name, var_pos, reference, trans, sanger, hg_version, ensembl):
        # adds all scrapped data to a list, which is written to an output file if the 
        # option is selected
       # try:
            # check each individual line of the file for CUSTOM ERRORS
            error_check = reference.handle_argument_exception(var_pos)
            
            # check if var_pos is a GENOMIC REGION, else construct one from var_pos
            seq_range = reference.create_region(var_pos)
            
            # use UCSC to get the genomic ranges DNA sequence
            sequence = reference.get_region_info(seq_range)
            
            # assess whether the var pos base in the sanger trace (seq_file) is
            # different to the reference base 
            sanger_sequence = sanger.match_with_seq_file(sequence)
            if isinstance(sanger_sequence, tuple):# if sanger_sequence is tuple else
                ref_base = sanger_sequence[1] 
                sanger_base = sanger_sequence[2] 

                # compare the reqerence var_pos base and the sanger var_pos base
                compare = CompareSeqs.compare_nucleotides(ref_base,sanger_base) 
            
            # detrmine whether to transcribe or translate to RNA or PROTEIN
            sequence = str(trans.get_rna_seq(sequence))
            sequence = str(trans.get_protein_seq(sequence))

            # determine whether to give a HEADER
            header = reference.header_option(seq_name,var_pos,
                                           seq_range,sequence)
            
            # get gene information for the variant position
            gene_info = ensembl.get_gene_info()
            if gene_info:
                gene_name, gene_id, gene_type, gene_range = gene_info
                transcript = ensembl.get_canonical_transcript(gene_name)
                exon_info = get_exon_number(transcript, hg_version, var_pos)
                exon_id, intron, exon = exon_info

            

            if isinstance(sanger_sequence, tuple) and gene_info:
                print("\n".join((header,"Reference Sequence:\t"+sequence,
                                 "Sanger Sequence:\t"+sanger_sequence[0],
                                 compare[0],"\n")))
                return("\t".join((seq_name, var_pos, seq_range, gene_name,
                                  gene_id, gene_type, gene_range, transcript, 
                                  exon_id, intron, exon, ref_base,
                                  sanger_base, str(compare[1]))))

            elif not gene_info and isinstance(sanger_sequence, tuple):
                print("\n".join((header, "Reference Sequence:\t"+sequence,
                                 "Sanger Sequence:\t"+sanger_sequence[0],
                                 compare[0], "\n")))
            elif not gene_info:
                print("\n".join((header, "Reference Sequence:\t"+sequence,"\n")))

            else:
                print("\n".join((header,"Reference Sequence:\t"+sequence, "\n")))
                return("\t".join((seq_name, var_pos, seq_range, gene_name, gene_id,  
                                  gene_type, gene_range, "-", "-",str(0))))


        #except WrongHGversion:
        #    print("Human genome version "+hg_version+" not recognised")
        #    sys.exit(0)
       # except TypographyError:
        #    print("Only one colon and no more than one comma/dash is allowed for "
         #         +var_pos+" in "+seq_name+"\n")    
       # except ErrorUCSC:
        #    print(var_pos+" in "+seq_name+" is not recognised by UCSC"+"\n")


               
if __name__ == '__main__':
    main()         

# main("in.txt","b", 2, 20,"hg19","\t","y")

