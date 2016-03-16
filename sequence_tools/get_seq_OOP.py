from __future__ import division
import os, sys,re, urllib2, click


class Processing():
    '''
    '''
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
            '''
            '''        
            if self.hg_version not in ["hg16","hg17","hg18","hg19","hg38"]:
                raise click.ClickException("Human genome version "+self.hg_version+" not recognised")
                sys.exit(0)
                
            if not int(self.upstream) or not int(self.downstream):
                raise click.ClickException("Ensure the string given for upstream "+
                "and/or downstream are integers:\n\tupstream: "+str(self.upstream)+
                "\t downstream: "+str(self.downstream))
                sys.exit(0)
                
            if os.path.isfile(self.input_file) is False:
                if self.input_file.count(",") > 1:
                    return("too many commas in "+self.input_file)
                
                if self.input_file.count(":") < 1 or self.input_file.count(":") >1:
                     return("A single colon is required to seperate"+ 
                     "the chromosome and position numbers in the variant position: "+
                     self.input_file)
                     
                    
                    
            if self.dash not in ["Y","N","y","n"]:
                raise click.ClickException("Only Y or N can be used with the --dash argument")
                sys.exit(0)
                
                
    def create_region(self,var_pos):
            ''' use the variant position given, add and subtract the 
                numbers given in upstream and downstream respectively
                to return a genomic range.
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
            
            Returns
            --------------------------------------------------------------------
            answer -- sequence associated with the given seq_range AKA genomic range
            
        '''
        # scrape for the sequence associated with the seq_range AKA genomic region 
        test = urllib2.urlopen("http://genome.ucsc.edu/cgi-bin/das/"+self.hg_version+"/dna?segment="+seq_range)
        search = test.read()
        refined_search = re.findall(r"[tacg{5}].*",search)
        confirm_error = re.findall(r"nnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnn",search)
        
        # BELOW IS REDUNDANT AS IT IS HANDLED BELOW AND IN THE TOP FUNCTION
        if confirm_error:
            return "Variant not recognised by UCSC"
            
        # filters for elements which only contain nucleotides and concatenate
        seqs = [s for s in refined_search if not s.strip("tacg")] 
        seq = "".join(seqs)
        
        
        # flank the base associated with the variant position with dashes
        try:
            if self.dash == "Y" or self.dash == "y":
                downstream = seq[:self.upstream]
                var = seq[self.upstream]
                upstream = seq[self.upstream+1:len(seq)]
                answer = "".join((downstream,"-",var,"-",upstream))
                return answer
        
        except IndexError:
            return "Variant not recognised by UCSC"
        
        # return sequence without dashes
        if self.dash == "N" or self.dash == "n":
            return seq
        
        # error message
        else:
            print "only Y or N can be given for --dash option"
             

#### exception handeling and commenting still needed
#### CREATE ANOTHER FUNCTION THAT HANDLES ALL ARGUMENT ERRORS


@click.command()
@click.option('--input_file')
@click.option('--output_file',default="get_seq_output.txt")
@click.option('--upstream', default=20)
@click.option('--downstream',default=20)
@click.option('--hg_version',default="hg19")
@click.option('--delimiters',default="\t")
@click.option('--dash',default="Y")

def get_seq(input_file, 
            output_file="get_seq_output.txt",
            upstream=20, 
            downstream=20, 
            hg_version="hg19", 
            delimiters="\t",
            dash="Y"
            ):
                     
        ''' Produce a sequence using the UCSC DAS server from a variant postion and 
            number of bases upstream and downstream from said position.
        
            Arguments
            --------------------------------------------------------------------
            input_file -- a variant name and variant position tab deliminated file
            
            output_file -- created only if the input is a file opposed to a string
                           if no name is sepcified for output_file, then it is 
                           automatically named get_seq_output.txt
                           
            upstream -- number of bases to get upstream from the given variant position
            
            downstream -- number of bases to get downstream from the given varant position
            
            hg_version -- the version of the human genome to scrape the sequence from
                          i.e. hg19, hg38
                          
            dash -- choose whether to include dashes flanking the variant positions base
                    or not i.e. N = no, Y = yes. defaulted to Y
                          
            Returns
            --------------------------------------------------------------------
            sequence -- sequence associated with the given seq_range AKA genomic range
                        with a dash flanking the base associated with the variant position
                        
            Examples
            --------------------------------------------------------------------
            
            
        '''
            
        # handle any errors in the arguments    
        process = Processing(input_file,output_file,upstream,downstream,hg_version,delimiters,dash)
        process.handle_argument_exception()
        
        # if input_file is a file
        if os.path.isfile(input_file) is True:
            output=open(output_file,"w")
            for changes in open(input_file,"r+"):
                
                changes = changes.split(delimiters)
                seq_name = changes[0]
                var_pos = changes[1]
                
                process = Processing(var_pos,output_file,upstream,downstream,hg_version,delimiters,dash)
                error_check = process.handle_argument_exception()
                if error_check is not None:
                    print error_check
                    continue
                
                seq_range = process.create_region(var_pos)
                if seq_range is None:
                    continue
                    
                answer = process.get_region_info(seq_range)
                if answer is None:
                    continue
                    
                sequence = "\t".join((seq_name,seq_range,answer))
                print sequence
                output.write(sequence+"\n")
                
            output.close()
        
        # if input_file is a genomic range
        if os.path.isfile(input_file) is False and re.search(r",",input_file):
            dash = "N"
            seq_range = re.sub(r'[^0-9:,]','',input_file) 
            sequence = process.get_region_info(seq_range)
            click.echo(sequence)
            return sequence
        
        # if input_file is a genomic position
        if os.path.isfile(input_file) is False:
            var_pos = re.sub(r'[^0-9:]','',input_file)
            seq_range = process.create_region(var_pos)
            sequence = process.get_region_info(seq_range)
            click.echo(sequence)
            return sequence


    
        

                
if __name__ == '__main__':
    get_seq()         

#get_seq("in.txt","b", 2, 20,"hg19","\t","y")

