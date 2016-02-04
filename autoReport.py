import openpyxl

# open the variant validation template workbook and worksheet
template = openpyxl.load_workbook("VariantConfirmationReport Template_BU.xlsx")
templates = template.get_sheet_by_name('Sheet1')

# open variant validation database
wb = openpyxl.load_workbook("All_Yale_&_UK_Variants.xlsx",data_only=True)
variants = wb.get_sheet_by_name('All_Variants') # or wb["All_variants"]
mut_ids = wb["Mutations ID"]

# list of variant aliases that will be used as input for the test
variant_list = [variants.cell(row=i,column=1).value for i in range(3,13)]
test_var = variant_list[5]

# write the variant alias to the report
templates["N12"] = test_var

# search the database and find the inputted variant alias in database. If found,
# output its row number
for row_number in range(1,200):
    if test_var == variants.cell(row=row_number,column=1).value:
        row = row_number

# use the above row number to retrieve basic specific indformation related 
# to this variant
for i in range(1,100):
    column = variants.cell(row=2,column=i).value
    get_info =variants.cell(row=row,column=i).value
     
    if "Sample_Name" == column:
        sam = get_info
    if "Gene" == column:
        gene = get_info
    if "Exon_No." == column:
        exon = get_info
    if "HGVSc" == column:
        HGVSc = get_info
    if "HGVSp" == column:
        HGVSp = get_info
    if "Variant_Position" == column:
        var_pos = get_info
    if "Allele_Balance" == column:
        AB = get_info
    if "Allele_Depth(REF)" == column:
        AD_REF = get_info
    if "Allele_Depth(ALT)" == column:
        AD_ALT = get_info
    if "Allele_Frequency_ESP(%)" == column:
        AF_ESP = get_info
        if AF_ESP is None:
            AF_ESP = "-"
    if "Allele_Frequency_ExAC(%)" == column:
        AF_ExAC = get_info
        if AF_ExAC is None:
            AF_ExAC = "-"
    if "Allele_Frequency_dbSNP(%)" == column:
        AF_dbSNP = get_info
        if AF_dbSNP is None:
            AF_dbSNP = "-"
    if "Variant_Found" == column:
        found = get_info
        if "Y" or "y" in found:
            found = "Y"
        else:
            found = "N"
                    
    if "Variant_Class" == column:
        comment = get_info
    if "First_Publication" == column:
        pub = get_info
    if "Category" == column:
        category = get_info 
        
   

# Handeling comments here          
if category == "HGMD" and "DM?" in comment:
    full_comment = comment+"\n"+pub+"\n"+"\n"
elif category == "HGMD" or category == "ClinVarPathogenic":
    full_comment = comment+"\n"+pub
    print full_comment

    
    

       
                      
# add this information to the template file
templates["E4"] =sam
templates["E7"] = gene
templates["E8"] = exon
templates["E9"] = HGVSc
templates["E10"] = HGVSp
templates["E11"] = var_pos
templates["E12"] = AB
templates["E13"] = str(AD_REF)+","+str(AD_ALT)
templates["E14"] = str(AF_ESP)+"\t"+str(AF_ExAC)+"\t"+str(AF_dbSNP)
templates["E15"] = found
templates["E16"] = full_comment
        
        

        
    



        

    
        

template.save("crazy_fun_yo.xlsx")