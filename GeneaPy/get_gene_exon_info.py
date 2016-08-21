import re, traceback, logging, sys
from GeneaPy import Ensembl


def gene_transcript_exon(position, hg_version, transcript=None):
    ''' Use the Ensembl module to return all gene, transcript and exon infomation
        relating to the given position

        returns a tuple of tuples:
            ((gene_name, gene_id, gene_type, gene_range), transcript,
             (exon_id, intron, exon))
    '''
    try:
        # defualt values
        gene_info = exon_info = ("-", "-", "-")
        transcript = "-"


        # convert a genomic region to a position in the middle of the region
        if "-" in position and ":" in position:
            position = create_position(position)

        ensembl = Ensembl.ScrapeEnsembl(position, hg_version)
        gene_info = ensembl.get_gene_info()
        
        
        if isinstance(gene_info, tuple):
            gene_name, gene_id, gene_type, gene_range = gene_info
            transcript = ensembl.get_canonical_transcript(gene_name)

            if transcript:
                exon_info = get_exon_number(transcript, hg_version, position)
                exon_id, intron, exon = exon_info

            else:
                transcript = "-"

        return (gene_info, transcript, exon_info)

    except TypeError as e:
        # get traceback from a string so one can identify what is causing the error
        e_type, e_value, e_traceback = sys.exc_info()
        sam = traceback.format_exception(e_type, e_value, e_traceback)
        if "exon_info" in sam[1]:
            msg = " ".join(("ERROR: No exon information found for",
                           position,"in",hg_version))
            return msg


def create_position(position):
    ''' From a genomic range, produce a position which is directly in 
        the middle of the given range
    '''
    
    split_region = re.split(r'[:-]', position)
    chrom = split_region[0]
    middle_region_position = str(round((int(split_region[2])+
                                        int(split_region[1])) / 2))
    position = "".join((chrom, ":", middle_region_position))
    return position




def get_exon_number(transcript, hg_version, pos):
    ''' From a transcript and genomic position, determine the exon number
        using the ExonInfo class within the Ensmebl module
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


#print(gene_transcript_exon("15:48762884", "hg19"))
#print(gene_transcript_exon("15:48762884", "hg38"))
#print(gene_transcript_exon("16:15812194", "hg19"))
#print(gene_transcript_exon("16:15812194", "hg38"))
#print(gene_transcript_exon("3:123418990", "hg19"))
#print(gene_transcript_exon("3:123418990", "hg38"))
