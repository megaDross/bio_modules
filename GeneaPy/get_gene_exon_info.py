import re
import Ensembl

def gene_transcript_exon(position, hg_version, transcript=None):
    # defualt values
    gene_info = exon_info = transcript = "-"


    # convert a genomic region to a position in the middle of the region
    if "-" in position and ":" in position:
        split_region = re.split(r'[:-]', position)
        chrom = split_region[0]
        middle_region_position = str(round((int(split_region[2])+
                                           int(split_region[1])) / 2))
        position = "".join((chrom, ":", middle_region_position))
    
    ensembl = Ensembl.ScrapeEnsembl(position, hg_version)
    gene_info = ensembl.get_gene_info()
    
    
    if isinstance(gene_info, tuple):
        gene_name, gene_id, gene_type, gene_range = gene_info
        transcript = ensembl.get_canonical_transcript(gene_name)

        if transcript:
            exon_info = get_exon_number(transcript, hg_version, position)
            exon_id, intron, exon = exon_info

    
    return (gene_info, transcript, exon_info)


def get_exon_number(transcript, hg_version, pos):
    ''' From a transcript and genomic position, determine the exon number
        using the ExonInfo class
    '''
    
    ensembl = Ensembl.ExonInfo(transcript, hg_version, pos)

    exon_dics = ensembl.request_ensembl()
    exon_region = ensembl.all_exon_regions(exon_dics)
    exon_id = ensembl.get_exon_id(exon_region) # filter for exon id in which variant is within
    if not exon_id:
        sorted_exon_regions = sorted(exon_region)
        intron_region = ensembl.all_intron_regions(sorted_exon_regions)
        intron_num = ensembl.intron_number(intron_region)
        
        if intron_num is None:
            return (pos,transcript,"NO INTRON/EXON MATCHED")
        else:
            last_exon = ensembl.total_exons(exon_dics)
            last_intron = int(last_exon)-1
            return ("-",str(intron_num)+"/"+str(last_intron),"-")
    else:
        exon_num = ensembl.exon_number(exon_dics,exon_id)    # use th exon_id to get exon number
        last_exon = ensembl.total_exons(exon_dics)           # get total exons of transcript
        
        
        return (exon_id[0], "-" ,str(exon_num)+"/"+str(last_exon))



