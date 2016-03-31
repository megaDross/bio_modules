import openpyxl, os
from openpyxl.drawing.image import Image
from get_variant_information import *

#1- add better comments
#2- filenames that are the same and if true then add a "gene2" some where


def produce_variant_report(variant_alias_list):
    ''' For each variant alias, extract the approriate variant and mutation
        information and append them to the variant confirmation template
        
            variant_alias_list: should be in a text file where each new line
                                contains a different variant alias or a single
                                variant alias.
    '''

    template = openpyxl.load_workbook("VariantConfirmationReport Template_BU.xlsx")
    templates = template.get_sheet_by_name('Sheet1')
    
    variants = open_variant_validation_spreadsheets()[0]
    mutations = open_variant_validation_spreadsheets()[1]
    
    if variant_alias_list.endswith(".txt"):
        var_alias_list = [ var.rstrip() for var in open(variant_alias_list)]
    else:
        var_alias_list = [variant_alias_list]
    
    for var_alias in var_alias_list: 
        variant_row = get_row(var_alias,variants)
        variant_info = get_variant_info(variant_row,variants,variant_header)

        if variant_header.get("Mutation_ID") != "-":  
            mut_row = get_row(variant_header.get("Mutation_ID"),mutations,2)
            mut_info = get_variant_info(mut_row,mutations,mutation_header)
        
        fill_report(template,templates,var_alias)
        
        transcript_id = variant_header.get("HGVSc")
        hgvs = transcript_id.split(":")[1]
        

        if (variant_header.get("Category") == "ClinVarPathogenic") or ("HGMD"):
            hgmd_clinvar_comment(templates)
            templates["B8"] = "Exon"
        if variant_header.get("Category") == "Gly-X-Y":
            glyxy_comment(templates)
            templates["B8"] = "Exon"
        if variant_header.get("Category") == "LOF":
            lof_comment(templates)
            templates["B8"] = "Exon"
            
        if "-" in hgvs or "+" in hgvs:
            templates["E16"] = "Splicing variant. This will require further investigation"
            templates["B8"] = "Intron"
        
        template.save(variant_header.get("Sample_Name")+"_"+variant_header.get("Variant_Alias")+"_"+"VariantConfirmationReport"+".xlsx")
        



def fill_report(template,templates,var_alias):
    ''' Append the information placed into the variant_header dict into the 
        template report
    '''
    templates["N12"]= var_alias
    templates["E4"] =variant_header.get("Sample_Name")
    templates["E7"] = variant_header.get("Gene")
    templates["E8"] = variant_header.get("Exon_No.")
    templates["E9"] = variant_header.get("HGVSc")
    templates["E10"] =variant_header.get("HGVSp")
    templates["E11"] =variant_header.get("Variant_Position")
    templates["E12"] =variant_header.get("Allele_Balance")
    templates["E13"] =str(variant_header.get ("Allele_Depth(REF)"))+","+str(variant_header.get("Allele_Depth(ALT)"))
    templates["E14"] =str(variant_header.get("Allele_Frequency_ESP(%)"))+"       "+\
    str(variant_header.get("Allele_Frequency_ExAC(%)"))+"       "+str(variant_header.get("Allele_Frequency_dbSNP(%)"))
    templates["E15"] =variant_header.get("Variant_Found")
    
    for sequences in os.listdir("C:\\cygwin64\\home\\dross11\\Playground\\sequences"):
        if var_alias in sequences:
            sequences = "C:\\cygwin64\\home\\dross11\\Playground\\sequences\\"+sequences
            img = Image(sequences)
            templates.add_image(img,"B36")
    
    template.save(variant_header.get("Sample_Name")+"_"+variant_header.get("Variant_Alias")+"_"+"VariantConfirmationReport"+".xlsx")



def hgmd_clinvar_comment(templates):
    ''' Add a comment associated with the HGMD or ClinVar accession number
        found in the database to the template report
    '''    
    comment = mutation_header.get("Report Variant class field")+"\n"+"\n"+\
                mutation_header.get("First Published")
    
    if mutation_header.get("Variant Class") == "DM?":
        full_comment = comment + "\n"+"\n"+"Date of Variant Class Change From DM to DM?: "+\
                        str(mutation_header.get("Date of variant class change"))+"\n"+\
                        mutation_header.get("Reason for Variant class change")
        templates["E16"] = full_comment
    else:
        templates["E16"] = comment



def glyxy_comment(templates):
    ''' Add a comment associated with GLY-X-Y variants to the template report
    '''
    comment = "This mutation is predicted to disrupt the collagen triple helical structure and is therefore likely to be pathogenic"
    templates["E16"] = comment


def lof_comment(templates):
    '''Add a comment associated with LOF variants to the template report. The 
       particular comment added is dependant upon the exon number in which the 
       variant lies within
    '''

    
    exon = variant_header.get("Exon_No.")
    intron = variant_header.get("Intron_No.")
        
    if exon != "-" or exon is not None:
        exon_num = int(exon.split("/")[0])
        exon_total = exon.split("/")[1]
        
        if exon_num == int(exon_total)-2 or exon_num == int(exon_total)-1 or exon_num == int(exon_total):
            templates["E16"] = "This mutations is expected to produce a truncated product"
        elif exon_num < int(exon_total)-2:
            templates["E16"] = "This mutation introduces a premature stop codon and is likely to be pathogenic"
        else:
            templates["E16"] = "LOF mutation present, but the outcome cannot be determined without exon numbering information"
    
    elif intron != "-" or intron is not None:
        templates["B8"] = "Intron"
        templates["E16"] = "This mutation affects the intron"
        
    else:
        templates["E16"] = "ERROR"
        
        
    


produce_variant_report("report_in.txt")

#print os.getcwd()
#print os.listdir("C:\\cygwin64\\home\\dross11\\Playground\\sequences")


### LOF requires EXON info or it throws a split error, try and get the 
### exon info going
    


