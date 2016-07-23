import re, requests, sys
from pyensembl import EnsemblRelease



class ScrapeEnsembl():
    ''' 
    '''
    def __init__(self, query, hg_version):
        self.query = query.replace("chr","")
        self.hg_version = ScrapeEnsembl.genome.get(hg_version) # convert to ensembl release
        self.hg = EnsemblRelease(self.hg_version) # convert to ensembl release object

    
    genome = {"hg19": 75, "hg38": 83}
    
    def get_gene_info(self):
        ''' Get the gene information at a given genomic position
        '''
         
        # check if the input is a genomic position or genomic range
        if re.search(r"[-:]", self.query) and self.query.replace(":","").isdigit():

            chrom = int(self.query.split(":")[0])
            pos = int(self.query.split(":")[1])
            gene_name = self.hg.gene_names_at_locus(contig=chrom, position=pos)
            if not gene_name:
                msg = " ".join(("No gene found at",self.query,"for genome version",
                                str(self.hg_version)))
                return msg 
            
            gene_info = self.hg.genes_by_name(gene_name[0])
            # gene_info[0].loaction doesn't work, hence the mess below
            gene_location = str(gene_info[0]).split(",")[-1][:-1].split("=")[1]

            gene_info = (gene_info[0].name, gene_info[0].id, 
                         gene_info[0].biotype, gene_location)
            
            return(gene_info)
    
    
    def get_canonical_transcript(self, gene_name):
        ''' Determine and return the canonical transcript of the given gene
        '''
        all_transcripts = self.hg.transcript_ids_of_gene_name(gene_name)
        all_transcript_details = [self.hg.transcript_by_id(x) for x in all_transcripts]
        protein_coding_transcripts = []
        for x in all_transcript_details:
            split_transcript_info = re.split(r"[=,]",str(x))
            transcript = split_transcript_info[1]
            transcript_type = split_transcript_info[9]
            location = split_transcript_info[-1][:-1]
            start = re.split(r"[:-]", location)[1]
            stop = re.split(r"[:-]", location)[2]
            size = int(stop) - int(start)
            if transcript_type == "protein_coding":
                protein_coding_transcripts.append((size,transcript,transcript_type)) 
        
        # sort by size and return the largest protein coding transcript
        if protein_coding_transcripts:    
            canonical_transcript = sorted(protein_coding_transcripts)[-1][1]
            return canonical_transcript



class ExonInfo():

    def __init__(self, transcript, hg_version, pos):
        self.transcript = transcript
        self.hg_version = hg_version
        self.pos = pos
    

    def request_ensembl(self):
        ''' Retrieve all expanded exon information associated with a given transcript ID
            in JSON format
        '''
        try:
            if self.hg_version == "hg19":
                grch = "grch37."
            elif self.hg_version == "hg38":
                grch = ""
            else:
                print("incompatible human genome version")
                sys.exit() 

            url = "".join(("http://",grch,"rest.ensembl.org/overlap/id/",self.transcript,
                           "?feature=exon;content-type=application/json;expand=1"))
            req = requests.get(url)
            req.raise_for_status()
            exon_dics = req.json()
            return exon_dics
        
        except requests.exceptions.RequestException as e:
            return e


    def all_exon_regions(self, exon_dics):
        '''Get every exon ID and its start and stop position which resides within the
           given transcript
        '''
        exon_region = set()
        for exon_dicts in exon_dics:
            if exon_dicts.get("Parent") == self.transcript:
                rank = exon_dicts.get("rank")
                start = exon_dicts.get("start")
                end = exon_dicts.get("end")
                exon = exon_dicts.get("exon_id")
                exon_start_end = (rank, exon, start, end)
                exon_region.add(exon_start_end)

        return exon_region
        

    def all_intron_regions(self,sorted_exon_regions):
        ''' Get all intron numbers, start and stop positions
        ''' 
        intron_region = set()
        
        for i in range(0,len(sorted_exon_regions)):
            for x in range(0,len(sorted_exon_regions)):
                if x == i +1:
                    intron_number = str(sorted_exon_regions[i][0])
                    end_pos_previous_exon = str(sorted_exon_regions[i][3])
                    start_pos_next_exon = str(sorted_exon_regions[x][2])
                    intron_info = (int(intron_number),end_pos_previous_exon,\
                    start_pos_next_exon,str(int(end_pos_previous_exon)-
                                            int(start_pos_next_exon)))
                    intron_region.add(intron_info)
        
        return sorted(intron_region)


    def intron_number(self,intron_region):
        ''' Returns the intron number in which the variant is within
        '''
        
        for intron in intron_region:
            if int(intron[3]) > 0:
             # no idea why the start and end pos are swapped, hence the need for if/else
                for x in range(int(intron[2]),int(intron[1])):          
                    if x == int(self.pos.split(":")[1]):
                        intron_num= intron[0]
                        return intron_num
            else:
                for x in range(int(intron[1]),int(intron[2])):
                    if x == int(self.pos.split(":")[1]):
                        intron_num= intron[0]
                        return intron_num  

                            
    def get_exon_id(self,exon_region):  
        ''' Get the exon id for which the variant position is within.
        
            Use the exon_region tuple to iterate between the start and stop 
            positions of each exon within the transcript until one number 
            matches the variant position number, then return the exon_id 
            associated with that position.
        '''  
        exon_id = []      
        for exons in sorted(exon_region):
            for x in range(exons[2],exons[3]):
                if x == int(self.pos.split(":")[1]):
                    exon_id.append(exons[1])
        
        return exon_id
        



    def exon_number(self,exon_dics,exon_id):
        ''' Returns the exon number in which the variant is within.
        
            use the exon_id to retrieve the matching dictionary and extract
            the rank (exon number) from that dictionary.
        '''
        for y in exon_dics:
            if y.get("exon_id") == exon_id[0] and y.get("Parent") == self.transcript:
                exon_number = y.get("rank")
       
        return exon_number

            
        
    def total_exons(self,exon_dics):
        ''' Get the last exon number in the transcript
        
            List comprehension which extracts all exon numbers found in all the 
            dictionarys with the transcript id. Extract and return the last element
            as the total_exons in the generated list.
        '''
        all_exon_ranks = [i.get("rank") for i in exon_dics 
                          if i.get("Parent") == self.transcript]
        total_exons = all_exon_ranks[-1]
        return total_exons













