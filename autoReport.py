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



def do_it(template_report,template_sheet,variant_database,variant_sheet, mutation_sheet,variant_alias_list):
    ''' Open the relevant template report, variant database, mutation database and list of variant aliases to be
        transformed into a variant confirmation report.
    '''
    ### THIS MAYBE A GOOD CANDIDATE FOR OOP
    template = openpyxl.load_workbook(template_report)
    templates = template.get_sheet_by_name(template_sheet)
    
    wb = openpyxl.load_workbook(variant_database,data_only=True)
    variants = wb.get_sheet_by_name(variant_sheet) # or wb["All_variants"]
    mutations = wb[mutation_sheet]
    
    var_alias_list = [ var.rstrip() for var in open(variant_alias_list)]
    
    for var_alias in var_alias_list:
        # extract the varinat info and place in variant_header dic  
        variant_row = get_row(var_alias,variants)
        variant_info = get_variant_info(variant_row,variants,variant_header)
        
        # extract the mutation info and place in mutation_header dic
        # if the mutation ID is a HGMD or ClinVar accession number
        if variant_header.get("Mutation_ID") != "-":  
            mut_row = get_row(variant_header.get("Mutation_ID"),mutations,2)
            mut_info = get_variant_info(mut_row,mutations,mutation_header)
        
        # append variant information to the report
        fill_report(template,templates,var_alias)
        
        # add a the approriate comment associated with the vars category
        if (variant_header.get("Category") == "ClinVarPathogenic") or ("HGMD"):
            hgmd_clinvar_comment(templates)
        if variant_header.get("Category") == "Gly-X-Y":
            glyxy_comment(templates)
        if variant_header.get("Category") == "LOF":
            lof_comment(templates)
        
        template.save(variant_header.get("Sample_Name")+"_"+variant_header.get("Variant_Alias")+"_"+"VariantConfirmationReport"+".xlsx")
        
    

do_it("VariantConfirmationReport Template_BU.xlsx",'Sheet1',"All_Yale_&_UK_Variants.xlsx",'All_Variants',"Mutations ID","test_in.txt")


### LOF requires EXON info or it throws a split error, try and get the 
### exon info going
    


