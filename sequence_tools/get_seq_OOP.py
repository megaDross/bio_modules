from __future__ import division
import os, sys,re, urllib2, click

#### exception handeling and commenting still needed

@click.command('get_seq')
@click.argument('input_file',nargs=1)
@click.option('--output_file',default="get_seq_output.txt", help='default: get_seq_output.txt')
@click.option('--upstream', default=20, help="number of bases to get upstream, default: 20") # default to an int makes option accept int only
@click.option('--downstream',default=20, help="number of bases to get downstream, default: 20")
@click.option('--hg_version',default="hg19", help="human geome version to get the sequence. default: hg19")
@click.option('--delimiters',default="\t", help="file delimiter. default: tab")
@click.option('--dash/--no_dash',default='n', help="dashes flanking the variant position base. default: --no_dash") # the slash in the option makes it a boolean

def get_seq(input_file, output_file="get_seq_output.txt", upstream=20, downstream=20, hg_version="hg19", delimiters="\t",dash="n"):
                     
        ''' 
        Produce a sequence using the UCSC DAS server from a variant postion and 
        number of bases upstream and downstream from said position.
         \b\n
        A file or string can be used as input; STRING: either a variant position
        or a genomic range deliminated by a comma FILE: deliminated file with the
        varinat name and the variant position   
         \b\n
        Example:\b\n
           python get_seq_OOP.py chr1:169314424 --dash --hg_version hg38\n
           python get_seq_OOP.py chr1:169314424,169314600\n
           python get_seq_OOP.py input.txt --output_file output.txt --dash\n
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
                
                # check each individual line of the file for errors
                process = Processing(var_pos,output_file,upstream,downstream,hg_version,delimiters,dash)
                error_check = process.handle_argument_exception()
                
                # if error is returned
                if error_check is not None:
                    print error_check
                    continue
                
                # if create_region is unsuccessful then move to next element
                seq_range = process.create_region(var_pos)
                if seq_range is None:
                    continue
                
                # if get_region_info is unsuccessful then move to next element
                answer = process.get_region_info(seq_range)
                if answer is None:
                    continue
                
                # concatenate the name and outputs from Class
                sequence = "\t".join((seq_name,seq_range,answer))
                print sequence
                output.write(sequence+"\n")
                
            output.close()
        
        
        # if input_file is a genomic range
        if os.path.isfile(input_file) is False and re.search(r",",input_file):
            # remove unwanted characters from the range
            seq_range = re.sub(r'[^0-9:,]','',input_file) 
            sequence = process.get_region_info(seq_range)
            click.echo(sequence)
            return sequence
        
        
        # if input_file is a genomic position
        if os.path.isfile(input_file) is False:
            # remove unwanted characters from var position
            var_pos = re.sub(r'[^0-9:]','',input_file)
            seq_range = process.create_region(var_pos)
            sequence = process.get_region_info(seq_range)
            click.echo(sequence)
            return sequence



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
                raise click.ClickException("Human genome version "+self.hg_version+" not recognised")
                sys.exit(0)
                
            if os.path.isfile(self.input_file) is False:
                if self.input_file.count(",") > 1:
                    return("too many commas in "+self.input_file)
                
                if self.input_file.count(":") < 1 or self.input_file.count(":") >1:
                     return("A single colon is required to seperate"+ 
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
        
        # BELOW IS REDUNDANT AS IT IS HANDLED BELOW AND IN THE TOP FUNCTION
        #if confirm_error:
        #    return "Variant not recognised by UCSC"
            
        # filters for elements which only contain nucleotides and concatenate
        seqs = [s for s in refined_search if not s.strip("tacg")] 
        seq = "".join(seqs)
        
        
        # flank the base associated with the variant position with dashes
        try:
            if self.dash:
                downstream = seq[:self.upstream]
                var = seq[self.upstream]
                upstream = seq[self.upstream+1:len(seq)]
                answer = "".join((downstream,"-",var,"-",upstream))
                return answer
        
        except IndexError:
            return "Variant not recognised by UCSC"
        
        # return sequence without dashes
        if not self.dash:
            return seq
            
                
if __name__ == '__main__':
    get_seq()         

#get_seq("in.txt","b", 2, 20,"hg19","\t","y")

