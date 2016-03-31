from get_variant_information import *

output = open("liverpool_variant_info.txt","w")
output.write("Sample_Name"+"\t"+"Category"+"\t"+"Gene"+
             "\t"+"Exon_No."+"\t"+"GC%"+"\t"+"AD(ALT)"
             +"\t"+"AB"+"\t"+"Variant_Found"+"\n")

variants = open_variant_validation_spreadsheets()[0]

answer = []
sams = []

for i in range(1,428):
    get_variant_info(i,variants,variant_header)
    sam = str(variant_header.get("Sample_Name"))
    if sam[:2] == "26":
        UID = variant_header.get("UID/Index")
        sams.append(UID)
    
for sam in sams:
    row = get_row(sam,variants,3)
    info = get_variant_info(row,variants,variant_header)
    answer.append(variant_header.get("Variant_Alias") +"\t "+variant_header.get("Category")+"\t "+
                    str(variant_header.get("Gene"))+"\t "+
                    str(variant_header.get("Exon_No."))+"\t "+
                    str(round(variant_header.get("GC%_Amplicon"),2))+"%"+"\t "+
                    str(variant_header.get("Allele_Depth(ALT)"))+"\t"+
                    str(variant_header.get("Allele_Balance"))+"\t "+
                    variant_header.get("Variant_Found"))

number_sams = len(answer)
print number_sams
for i in answer:
    print i
    output.write(i+"\n")
output.close()
    
    
    

    
    