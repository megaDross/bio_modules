import openpyxl

variant_header = {
"Sample_Name": '', "Gene": '',
"Exon_No.": '', "HGVSc": '',
"HGVSp": '', "Variant_Position": '',
"Allele_Balance": '', "Allele_Depth(REF)": '',
"Allele_Depth(ALT)": '', "Allele_Frequency_ESP(%)": '',
"Allele_Frequency_ExAC(%)": '', "Allele_Frequency_dbSNP(%)": '',
"Variant_Found": '', "Mutation_ID": '',"Reason for Variant class change": '',
"Category":'', "Variant_Alias":''
}

mutation_header = {
"Report Variant class field": '', "First Published": '',
"Reason for Variant class change": '', "Date of variant class change": '',
"Variant Class": ''
}


def produce_variant_report(variant_alias_list):
    ''' For each variant alias, extract the approriate variant and mutation
        information and append them to the variant confirmation template
        
            variant_alias_list: should be in a text file where each new line
                                contains a different variant alias or a single
                                variant alias.
    '''

    open_spreadsheets()
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
        
        if (variant_header.get("Category") == "ClinVarPathogenic") or ("HGMD"):
            hgmd_clinvar_comment(templates)
        if variant_header.get("Category") == "Gly-X-Y":
            glyxy_comment(templates)
        if variant_header.get("Category") == "LOF":
            lof_comment(templates)
        
        template.save(variant_header.get("Sample_Name")+"_"+variant_header.get("Variant_Alias")+"_"+"VariantConfirmationReport"+".xlsx")
        

def open_spreadsheets():
    '''Open the template report workbook, variant database workbook and
        the relevant variant and mutation sheets
    '''
    global template
    global templates
    global variants
    global mutations
    
    template = openpyxl.load_workbook("VariantConfirmationReport Template_BU.xlsx")
    templates = template.get_sheet_by_name('Sheet1')
    
    wb = openpyxl.load_workbook("All_Yale_&_UK_Variants.xlsx",data_only=True)
    variants = wb.get_sheet_by_name('All_Variants') # or wb["All_variants"]
    mutations = wb["Mutations ID"]



def get_row(query,sheet,column_num=1):
    ''' Search the database/sheet and match with the query/variant_alias
        if found, output its row number in the database
    '''
    for row_number in range(1,500):
        if query == sheet.cell(row=row_number,column=column_num).value:
            return row_number
    

def get_variant_info(row,sheet,dic):
    ''' Extract information associated with the query inputted in get_row() 
        from the database and append it to the items in variant_header dict
    '''
    for i in range(1,100):
        column = sheet.cell(row=2,column=i).value
        get_info =sheet.cell(row=row,column=i).value
        if column in dic:
            dic[column]=get_info
            if dic.get(column) is None:
                dic[column] = "-"
    return dic


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
    
    template.save(variant_header.get("Sample_Name")+"_"+variant_header.get("Variant_Alias")+"_"+"VariantConfirmationReport"+".xlsx")



def hgmd_clinvar_comment(templates):
    ''' Add a comment associated with the HGMD or ClinVar accession number
        found in the database to the template report
    '''    
    comment = mutation_header.get("Report Variant class field")+"\n"+"\n"+\
                variant_header.get("Category")+" ID: "+mutation_header.get("Variant Class")+"\n"+\
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
    if exon != "-" or exon is not None:
        exon_num = exon.split("/")[0]
        exon_total = exon.split("/")[1]
        if (exon_num == int(exon_total)-2) or int(exon_total)-1 or int(exon_total):
            templates["E16"] = "This mutations is expected to produce a truncated product"
        elif exon_num < int(exon_total)-2:
            templates["E16"] = "This mutation introduces a premature stop codon and is likely to be pathogenic"
    else:
        templates["E16"] = "LOF mutation present, but the outcome cannot be determined without exon numbering information"




    

produce_variant_report("test_in.txt")


### LOF requires EXON info or it throws a split error, try and get the 
### exon info going
    


