from __future__ import division
import os,re, urllib2, click

class IncorrectVariant():
    pass

#### exception handeling and commenting still needed

@click.command()
@click.option('--input_file')
@click.option('--output_file',default="region_extractor_output.txt")
@click.option('--upstream', default=20)
@click.option('--downstream',default=20)
@click.option('--hg_version',default="hg19")
@click.option('--delimiters',default="\t")
@click.option('--dash',default="Y")

def region_extractor(input_file, 
                     output_file="region_extractor_output.txt",
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
                           automatically named region_extractor_output.txt
                           
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
        # if input_file is a file
        if os.path.isfile(input_file) is True:
            output=open(output_file,"w")
            for changes in open(input_file,"r+"):
                
                changes = changes.split(delimiters)
                seq_name = changes[0]
                var_pos = changes[1]
                
                seq_range = create_region(var_pos,upstream,downstream)
                answer = get_region_info(seq_range,upstream,downstream,hg_version,dash)
                sequence = "\t".join((seq_name,seq_range,answer))
                print sequence
                output.write(sequence+"\n")
                
            output.close()
        
        # if input_file is a genomic range
        if os.path.isfile(input_file) is False and re.search(r",",input_file):
            dash = "N"
            var_pos = re.sub(r'[^0-9:,]','',input_file) 
            sequence = get_region_info(var_pos,upstream,downstream,hg_version,dash)
            click.echo(sequence)
            return sequence
        
        # if input_file is a genomic position
        if os.path.isfile(input_file) is False:
            var_pos = re.sub(r'[^0-9:]','',input_file) 
            seq_range = create_region(var_pos,upstream,downstream)
            sequence = get_region_info(seq_range,upstream,downstream,hg_version,dash)
            click.echo(sequence)
            return sequence
        
        


def create_region(var_pos,upstream=20, downstream=20):
            ''' use the variant position given, add and subtract the 
                numbers given in upstream and downstream
                to return a genomic range.
                
                var_pos -- variant position
                
                upstream -- number of bases to get upstream from the given variant position
            
                downstream -- number of bases to get downstream from the give n variant position
            
            '''
            
            # create a genomic range 
            nospace = var_pos.replace(" ","")
            chrom = nospace.split(":")[0]
            pos = nospace.split(":")[1]
            start_pos = int(pos) - upstream
            end_pos = int(pos) + downstream
            seq_range = chrom+":"+str(start_pos)+","+str(end_pos)
            return seq_range
            
    
        
def get_region_info(seq_range, upstream,downstream,hg_version,dash):
        ''' From a genomic range and human genome version, use UCSC DAS server
            to retrieve the sequence found in the given genomic range.
            
            Arguments
            -------------------------------------------------------------------
            seq_range -- a genomic range i.e. chr1:2400000,2408002
            
            upstream -- number of bases to get upstream from the given variant position
            
            downstream -- number of bases to get downstream from the give n variant position
            
            hg_version -- version of the human genome i.e. hg19, hg38
            
            dash -- choose whether to include dashes flanking the variant positions base
                    or not i.e. N, Y. defaulted to Y
            
            Returns
            --------------------------------------------------------------------
            answer -- sequence associated with the given seq_range AKA genomic range
            
        '''
        # scrape for the sequence associated with the seq_range AKA genomic region 
        test = urllib2.urlopen("http://genome.ucsc.edu/cgi-bin/das/"+hg_version+"/dna?segment="+seq_range)
        search = test.read()
        refined_search = re.findall(r"[tacg{5}].*",search)
        confirm_error = re.findall(r"nnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnn",search)
        if confirm_error:
            raise IncorrectVariant("Variant not recognised by UCSC")
            
        # filters for elements which only contain nucleotides
        seqs = [s for s in refined_search if not s.strip("tacg")] 
        seq = "".join(seqs)
        
        
        # flank the base associated with the variant position with dashes
        if dash == "Y" or dash == "y":
            downstream = seq[:upstream]
            var = seq[upstream]
            upstream = seq[upstream+1:len(seq)]
            answer = downstream+"-"+var+"-"+upstream
            return answer
        
        # retrun sequence without dashes
        if dash == "N" or dash == "n":
            return seq
        
        # error message
        else:
            print "only Y or N can be given for --dash option"
        
                
if __name__ == '__main__':
    region_extractor()         


