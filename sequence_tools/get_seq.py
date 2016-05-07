from __future__ import division, print_function
import os, sys,re, click, requests, bs4
from useful_tools.transcription_translation import transcription, translation
from useful_tools.process_file import ProcessIO

class WrongHGversion(Exception):
    pass

class TypographyError(Exception):
    pass
    
class ErrorUCSC(Exception):
    pass

   
@click.command('get_seq')
@click.argument('input_file',nargs=1, required=False)
@click.option('--output_file',default=None, help='give an output file, requires input to be a file')
@click.option('--upstream', default=20, help="number of bases to get upstream, default: 20") # default to an int makes option accept int only
@click.option('--downstream',default=20, help="number of bases to get downstream, default: 20")
@click.option('--hg_version',default="hg19", help="human genome version. default: hg19")
@click.option('--dash/--no_dash',default='n', help="dashes flanking the variant position base. default: --no_dash") # the slash in the option makes it a boolean
@click.option('--header/--no_header',default='n',help="header gives metadata i.e. sequence name etc.")
@click.option('--transcribe/--nr',default='n',help="transcribe into RNA sequence")
@click.option('--translate/--np',default='n',help="translate RNA seq into protein seq")
@click.option('--rc/--no_rc',default='n',help="reverse complement the DNA")

# what use is there in just transcribing and translating for the sake of it? perhaps use an intiation codon finder which can be used o run through a DNA sequence pior to transcripion and transcribing from said site. Also something which checks the dna seq is a multiple of 3 would also be useful. Perhaps something where I could push the rna/protein seq to find which gene exon etc. is or maybe nBLAST, pBLAST etc which I am sure BioPython will have something for parsing to.

def get_seq(input_file, output_file=None, upstream=20, downstream=20, hg_version="hg19"
            , dash="n", header="n", transcribe="n", translate="n",
            rc='n'):
        '''
    Produce a sequence using the UCSC DAS server from an inputted genomic 
	postion and defined number of bases upstream and downstream from said 
	position. A genomic range can be used in place of a genomic position
	and renders the upstream/downstream options irrelevant. An input file
	can have a mixture of genomic positions and genomic ranges.
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
        # parse all arguments into the Processing Class 
        process = Processing(input_file,output_file,upstream,downstream,hg_version
                             ,dash,transcribe,translate,rc,header)
        
        # parse IO into ProcessIO Class to determine and process input type accordingly
        process_io = ProcessIO(input_file,output_file)

        # determine whether input_file is a string or file and alter accordingly
        input_file = process_io.process_input()
        
        # adds all scrapped data to a list, which is written to an output file if the 
        # option is selected
        sequence_data = []
        
        for changes in input_file:
            try:
                # split data up into name and position, remove ambigious characters 
                # from position
                seq_name = changes.split("\t")[0]
                var_pos = re.sub(r'[^0-9:,-]','',changes.split("\t")[1])
                
                # check each individual line of the file for CUSTOM ERRORS
                error_check = process.handle_argument_exception(var_pos)
                
                # check if var_pos is a GENOMIC REGION, else construct one from var_pos
                seq_range = process.create_region(var_pos)
                
                # use UCSC to get the genomic ranges DNA sequence
                answer = process.get_region_info(seq_range)
                
                # detrmine whether to transcribe or translate to RNA or PROTEIN
                sequence = process.get_rna_seq(answer)
                sequence = process.get_protein_seq(sequence)
                
                # determine whether to give a HEADER
                header_sequence = process.header_option(seq_name,var_pos,
                                                        seq_range,sequence)
                # append the output of each iteration into an empty list
                if output_file:
                    sequence_data.append(header_sequence)

            except WrongHGversion:
                print("Human genome version "+hg_version+" not recognised")
                sys.exit(0)
            except TypographyError:
                print("Only one colon and no more than one comma/dash is allowed for "
                      +var_pos+" in "+seq_name+"\n")    
            except ErrorUCSC:
                print(var_pos+" in "+seq_name+" is not recognised by UCSC"+"\n")
        
        # if OUTPUT_FILE option is selected then write all seq data to a file
        if output_file:
            process_io.write_to_output(sequence_data)
        
        # DO I WANT TO RETURN THIS?
        return sequence_data


@transcription
def transcribe_dna(dna):
    ''' Transcribe dna into rna 
    '''
    return dna

@translation
def translate_rna(rna):
    return rna




class Processing():

    def __init__(self,input_file,output_file,upstream, downstream, hg_version, 
                 dash, transcribe, translate,rc,header):
        
        self.input_file = input_file
        self.output_file = output_file
        self.upstream = upstream
        self.downstream = downstream
        self.hg_version = hg_version
        self.dash = dash
        self.transcribe = transcribe
        self.translate = translate
        self.rc = rc
        self.header = header

              
    def handle_argument_exception(self,var_pos):
        ''' Stores custom exceptions
        '''        
        
        if self.hg_version not in ["hg16","hg17","hg18","hg19","hg38"]:
            raise WrongHGversion("Human genome version "+self.hg_version+
                                 " not recognised")
            sys.exit(0)
            
        
        if var_pos.count(",") > 1 or var_pos.count("-") >1:
            raise TypographyError("too many commas in "+self.input_file)
            
            
        if var_pos.count(":") < 1 or var_pos.count(":") >1:
            raise TypographyError("A single colon is required to seperate"+\
                                  "the chromosome and position numbers in the"+\
                                  "variant position: "+self.input_file)
                                             
                
    def create_region(self,var_pos):
        ''' use the variant position given, add and subtract the 
            numbers given in upstream and downstream respectively
            from the given variant position to return a genomic range.
        '''
        # check if var_pos is a GENOMIC REGION, else construct one from var_pos
        if re.search(r"[,-]",var_pos):
            var_pos = var_pos.replace("-",",")
            return var_pos

        else:                    
            nospace = var_pos.replace(" ","")
            chrom = nospace.split(":")[0]
            pos = nospace.split(":")[1]
            start_pos = int(pos) - self.upstream
            end_pos = int(pos) + self.downstream
            seq_range = chrom+":"+str(start_pos)+","+str(end_pos)
            return seq_range
                
                
    def get_region_info(self, seq_range):
        ''' From a genomic range and human genome version, use UCSC DAS server
            to retrieve the sequence found in the given genomic range.
        '''
        # scrape for the sequence associated with the seq_range AKA genomic region 
        req = requests.get("http://genome.ucsc.edu/cgi-bin/das/"+self.hg_version+
                               "/dna?segment="+seq_range)
        req.raise_for_status()
        url = bs4.BeautifulSoup(req.text, features="xml").prettify()
        search = re.findall(r"[tacg{5}].*",url)
        
        # filters for elements which only contain nucleotides and concatenate
        seqs = [s for s in search if not s.strip("tacg")] 
        seq = "".join(seqs)
        if not seq:
            raise ErrorUCSC
        
        
        # flank the base associated with the variant position with dashes
        if self.dash:
            downstream = seq[:self.upstream]
            var = seq[self.upstream]
            upstream = seq[self.upstream+1:len(seq)]
            answer = "".join((downstream,"-",var,"-",upstream))
            return answer
    
        # return sequence without dashes
        if not self.dash:
            return seq


    def get_rna_seq(self,sequence):
        ''' determine whether to transcribe rna,
            based upon options selected
        '''
        # perorm transcription, reverse complement and/or translation depending
        # on which options have been selected
        if self.transcribe:
         
            if self.rc:
                rna = transcribe_dna(sequence.replace("-",""),"rc")
                return rna

            else:
                rna = transcribe_dna(sequence.replace("-",""))
                return rna
        else:
            # DNA
            return sequence


        
    def get_protein_seq(self,rna):
        ''' determine whether to translate protein,
            based upon options selected
        '''
        if self.translate:
            protein = translate_rna(rna)
            return protein
        else:
            return rna
    
    
    
    def header_option(self,seq_name,var_pos,seq_range,sequence):
        ''' determine whether to place a header
            above the returned sequence, based 
            upon options selected by user
        '''
        # concatenate the name and outputs from Class, determine whether to 
        # add a header
        if self.header:
            sequence = " ".join((">",seq_name,var_pos,seq_range,
                                 "\n",str(sequence),"\n"))
        else:
            sequence = "".join((str(sequence),"\n"))

        # output sequences to the screen and append to a list
        print(sequence)
        return(sequence)

         
                
if __name__ == '__main__':
    get_seq()         

# get_seq("in.txt","b", 2, 20,"hg19","\t","y")

