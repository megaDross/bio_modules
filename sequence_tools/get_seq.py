from __future__ import division, print_function
import os, sys,re, urllib2, click

class WrongHGversion(Exception):
    pass

class TypographyError(Exception):
    pass
    
class ErrorUCSC(Exception):
    pass
    
    
@click.command('get_seq')
@click.argument('input_file',nargs=1)
@click.option('--output_file',default=None, help='give an output file, requires input to be a file')
@click.option('--upstream', default=20, help="number of bases to get upstream, default: 20") # default to an int makes option accept int only
@click.option('--downstream',default=20, help="number of bases to get downstream, default: 20")
@click.option('--hg_version',default="hg19", help="human genome version. default: hg19")
@click.option('--delimiters',default="\t", help="file delimiter. default: tab")
@click.option('--dash/--no_dash',default='n', help="dashes flanking the variant position base. default: --no_dash") # the slash in the option makes it a boolean
@click.option('--header/--no_header',default='n',help="header gives metadata i.e. sequence name etc.")

def get_seq(input_file, output_file=None, upstream=20, downstream=20, hg_version="hg19", delimiters="\t",dash="n",header="n"):
        ''' 
        Produce a sequence using the UCSC DAS server from an inputted genomic 
	postion and defined number of bases upstream and downstream from said 
	position. A genomic range can be used in place of a genomic position
	and renders the upstream/downstream options irrelevant. An input file
	can have a mixture of genomic positions and genomic ranges.
         \b\n
        A file or string can be used as input.
	STRING: either a variant position or a genomic range deliminated by a comma
	FILE: deliminated file with the variant name and the variant position
         \b\n
        Example:\b\n
           get_seq chr1:169314424 --dash --upstream 200 --downstream 200\n
           get_seq chr1:169314424,169314600 --hg_version hg38\n
           get_seq input.txt --output_file output.txt --header --delimiters ,\n
        '''
        
        # if input is not a file, create a list
        if os.path.isfile(input_file) is False:
            input_file = ["query"+delimiters+input_file]
            return get_seq_data(input_file, output_file, upstream, downstream, hg_version, delimiters,dash,header)        
        
        # if input is a file, open file
        if os.path.isfile(input_file) is True:
            input_file = open(input_file,"r+")
            sequence_data = get_seq_data(input_file, output_file, upstream, downstream, hg_version, delimiters,dash,header)
            
            # if the --output_file option is used, write the sequence_data to a file
            if output_file is not None:
                output = open(output_file,"w")
                for seqs in sequence_data:
                    output.write(seqs)
                output.close()
                return output
                
            return sequence_data


def get_seq_data(input_file, output_file, upstream, downstream, hg_version, delimiters,dash,header):
        '''
        Read get_seq description, should be here (for the sake of click --help page it isn't). 
        ''' 
        # adds all scrapped data to a list, which is written to an output file if option is selected
        sequence_data = []
        
        for changes in input_file:
            try:
                # split data up into name and position, remove bad characters from position
                changes = changes.split(delimiters)
                seq_name = changes[0]
                var_pos = re.sub(r'[^0-9:,-]','',changes[1])
                
                # check each individual line of the FILE for custom errors
                process = Processing(var_pos,output_file,upstream,downstream,hg_version,delimiters,dash)
                error_check = process.handle_argument_exception()
                
                # check if var_pos is a genomic region, else construct one from var_pos
                if re.search(r"[,-]",var_pos):
                    var_pos = var_pos.replace("-",",")
                    seq_range = var_pos
                else:   
                    seq_range = process.create_region(var_pos)
                
                # use UCSC to get the genomic ranges sequence
                answer = process.get_region_info(seq_range)
                
                # concatenate the name and outputs from Class, determine whether to add a header
                if header:
                    sequence = " ".join((">",seq_name,var_pos,seq_range,"\n",answer,"\n"))
                if not header:
                    sequence = "".join((answer,"\n"))
                print(sequence)
                sequence_data.append(sequence)
                
            except WrongHGversion:
                print("Human genome version "+hg_version+" not recognised")
                sys.exit(0)
            except TypographyError:
                print("Only one colon and no more than one comma/dash is allowed for "
                        +var_pos+" in "+seq_name+"\n")    
            except ErrorUCSC:
                print(var_pos+" in "+seq_name+" is not recognised by UCSC"+"\n")
            
        return sequence_data
            
        
class Processing():

    def __init__(self,input_file,output_file,
                 upstream, downstream, hg_version, 
                 delimiters,dash):
        
            self.input_file = input_file
            self.output_file = output_file
            self.upstream = upstream
            self.downstream = downstream
            self.hg_version = hg_version
            self.delimiters = delimiters
            self.dash = dash
        
        
    def handle_argument_exception(self):
            ''' Stores custom exceptions
            '''        
            
            if self.hg_version not in ["hg16","hg17","hg18","hg19","hg38"]:
                raise WrongHGversion("Human genome version "+self.hg_version+" not recognised")
                sys.exit(0)
                
            
            if self.input_file.count(",") > 1 or self.input_file.count("-") >1:
                raise TypographyError("too many commas in "+self.input_file)
                
                
            if self.input_file.count(":") < 1 or self.input_file.count(":") >1:
                raise TypographyError("A single colon is required to seperate "+ 
                "the chromosome and position numbers in the variant position: "+
                self.input_file)
                                             
                
    def create_region(self,var_pos):
            ''' use the variant position given, add and subtract the 
                numbers given in upstream and downstream respectively
                from the given variant position to return a genomic range.
            '''
            
            # create a genomic range
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
        test = urllib2.urlopen("http://genome.ucsc.edu/cgi-bin/das/"+self.hg_version+"/dna?segment="+seq_range)
        search = test.read()
        refined_search = re.findall(r"[tacg{5}].*",search)
        confirm_error = re.findall(r"nnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnn",search)
        if confirm_error:
            raise ErrorUCSC
        
        # filters for elements which only contain nucleotides and concatenate
        seqs = [s for s in refined_search if not s.strip("tacg")] 
        seq = "".join(seqs)
        
        
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
            
                
if __name__ == '__main__':
    get_seq()         

#get_seq("in.txt","b", 2, 20,"hg19","\t","y")

