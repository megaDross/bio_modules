import sys, requests, bs4, re 

# pysam doesnt work with cygwin
if sys.platform == "cygwin":
    pass
else:
    import pysam


class WrongHGversion(Exception):
    pass

class TypographyError(Exception):
    pass
    
class ErrorUCSC(Exception):
    pass




class ScrapeSeq():
    ''' A collection of methods that are required for scrapping the reference
        sequence from UCSC
    '''
    def __init__(self,input_file,upstream, downstream, hg_version, header="Y"):
        self.input_file = input_file
        self.upstream = upstream
        self.downstream = downstream
        self.hg_version = hg_version
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
        ''' use the variant position given, add and subtract the numbers given 
            in upstream and downstream respectively from the given variant position 
            to return a genomic range.
        '''
        # PySam doesnt like chr preceding position/range
        var_pos = var_pos.replace('chr', '')

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

            http://www.biodas.org/documents/spec-1.53.html
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
        
        # split up scrapped sequence and make var upper
        downstream = seq[:self.upstream]
        var = seq[self.upstream]
        upstream = seq[self.upstream+1:len(seq)]
        answer = "".join((downstream.lower(),var.upper(),upstream.lower()))
        
        return answer

   
    def get_region_info_locally(self, seq_range, genome_path):
        ''' Use a locally downloaded genome instead of scrapping from UCSC
        '''
        # get the sequence associated with a given genomic range from a given fasta file
        chrom, start, end = tuple(re.split(r"[:,]", seq_range))
        chrom = "".join(("chr",chrom))
        genome = pysam.FastaFile(genome_path)

        # -1 is required otherwise the first base is missing, no idea why
        seq = genome.fetch(chrom, int(start)-1, int(end))

        # split up sequence and make var upper
        downstream = seq[:self.upstream]
        var = seq[self.upstream]
        upstream = seq[self.upstream+1:len(seq)]
        answer = "".join((downstream.lower(),var.upper(),upstream.lower()))
        
        return answer
       
    
    def header_option(self,seq_name,var_pos,seq_range,sequence, gene_name=None):
        ''' determine whether to place a header above the returned sequence, based 
            upon options selected by user
        '''
        # concatenate the name and outputs from Class, determine whether to 
        # add a header
        if self.header and gene_name:
            header = " ".join((">",seq_name,var_pos,seq_range,gene_name))
        elif self.header:
            header = " ".join((">",seq_name,var_pos,seq_range))
        else:
            header = ""

        # output sequences to the screen and append to a list
        return(header)


