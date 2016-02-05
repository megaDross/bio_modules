import openpyxl

# open the variant validation template workbook and worksheet
template = openpyxl.load_workbook("VariantConfirmationReport Template_BU.xlsx")
templates = template.get_sheet_by_name('Sheet1')

# open variant validation database
wb = openpyxl.load_workbook("All_Yale_&_UK_Variants.xlsx",data_only=True)
variants = wb.get_sheet_by_name('All_Variants') # or wb["All_variants"]
mutations = wb["Mutations ID"]






# search the database and find the inputted variant alias in database. If found,
# output its row number
def get_row(query,sheet,column_num=1):
    for row_number in range(1,200):
        if query == sheet.cell(row=row_number,column=column_num).value:
            return row_number


# Use a dictionary where the keys will be used to extract column information
# from the above row number
variant_header = {
"Sample_Name": '', "Gene": '',
"Exon_No.": '', "HGVSc": '',
"HGVSp": '', "Variant_Position": '',
"Allele_Balance": '', "Allele_Depth(REF)": '',
"Allele_Depth(ALT)": '', "Allele_Frequency_ESP(%)": '',
"Allele_Frequency_ExAC(%)": '', "Allele_Frequency_dbSNP(%)": '',
"Variant_Found": '', "Mutation_ID": '',"Reason for Variant class change": '',
"Category":''
}

mutation_header = {
"Report Variant class field": '', "First Published": '',
"Reason for Variant class change": '', "Date of variant class change": '',
"Variant Class": ''
}

def get_variant_info(row,sheet,dic):
    for i in range(1,100):
        column = sheet.cell(row=2,column=i).value
        get_info =sheet.cell(row=row,column=i).value
        if column in dic:
            dic[column]=get_info
            
            if dic.get(column) is None:
                dic[column] = "-"
    return dic


var_alias_list = [ var.rstrip() for var in open("test_in.txt")]

for var_alias in var_alias_list:
    # extract the varinat info and place in variant_header dic  
    variant_row = get_row(var_alias,variants)
    variant_info = get_variant_info(variant_row,variants,variant_header)
    
    # extract the mutation info and place in mutation_header dic  
    mut_row = get_row(variant_header.get("Mutation_ID"),mutations,2)
    mut_info = get_variant_info(mut_row,mutations,mutation_header)
    
    
    # add this information to the template file
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
    
    
    # add a comment about the mutation
    if (variant_header.get("Category") == "ClinVarPathogenic") or ("HGMD"):
        comment = mutation_header.get("Report Variant class field")+"\n"+"\n"+\
                  variant_header.get("Category")+" ID: "+mutation_header.get("Variant Class")+"\n"+\
                  mutation_header.get("First Published")
        templates["E16"] = comment
        
        if mutation_header.get("Variant Class") == "DM?":
            full_comment = comment + "\n"+"\n"+"Date of Variant Class Change From DM to DM?: "+\
                           str(mutation_header.get("Date of variant class change"))+"\n"+\
                           mutation_header.get("Reason for Variant class change")
            templates["E16"] = full_comment
        
        
    
    
    
    template.save(variant_header.get("Sample_Name")+".xlsx")
    


    


