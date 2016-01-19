from __future__ import division
import csv
import re
import itertools

# db SNP testing ground

##INFO=<ID=CAF,Number=.,Type=String,Description="An ordered, comma 
##delimited list of allele frequencies based on 1000Genomes, starting 
##with the reference allele followed by alternate alleles as ordered in 
##the ALT column. Where a 1000Genomes alternate allele is not in the 
##dbSNPs alternate allele set, the allele is added to the ALT column.  
##The minor allele is the second largest value in the list, and was 
##previuosly reported in VCF as the GMAF.  This is the GMAF reported on 
##the RefSNP and EntrezSNP pages and VariationReporter">


HOME = "C:\cygwin64\home\dross11\\"
DB = "C:\cygwin64\home\dross11\VCF\DB\\"

var_file = csv.reader(open\
(HOME+"All_Yale_&_UK_Variants_160106.txt","rb")\
,delimiter="\t")

ExAC_test = csv.reader(open(DB+"output_ExAC.txt","rb").readlines()[1:]\
,delimiter="\t")

ESP_test = open(DB+"output_ESP.txt","rb")

db_test = open(DB+"output_dbSNP_all.txt","rb")

test = csv.reader(open\
(DB+"Small_test.vcf","rb"),delimiter="\t")

output = open(DB+"output_mutation_frequency.txt","w")

simplified = []

# reapproriate the file that contains all TAAD variants so that the
# variant alias, ref allele, alt allele, chromosome of the variant pos
# & position of the variant position are tab delimintated and stored in
# a list named simplified
for foo in var_file:
    CHR = foo[18].split(":")[0]
    pre_POS = foo[18].split(":")
    if len(pre_POS) == 2:
        POS = re.sub(r"chr","",pre_POS[1])
        rearranged_file = foo[0]+"\t"+foo[25]+"\t"+foo[26]+"\t"+re.sub(r"chr","",CHR)+"\t"+POS
        simplified.append(rearranged_file)

#checking the ESP database
for yale_var,db_var in itertools.product(csv.reader(simplified,delimiter="\t"),db_test):
    database_var= db_var.split("\t")
    
    if len(database_var) >2 and yale_var[3]+":"+yale_var[4] == database_var[0]+":"+database_var[1]:
        ALT_database = database_var[4].split(",")
        info = database_var[7]
        search_CAF = re.findall(r"CAF=.*;",info)
            # TAC = [alt allele count,ref allele count] in INFO fiel ESP* VCF files
            # AAC = alt allele count # RAC = ref allele count # AN = Total alleles
        AAC = str(search_CAF)#.split(";")[0][5:].split(",")[0]
        RAC = str(search_CAF)#.split(";")[0][5:].split(",")[1]
        print yale_var[0]+"\t"+yale_var[3]+":"+yale_var[4]+"\t"+AAC
        #AN = int(AAC) + int(RAC)
        
        #if len(ALT_database) == 1:
        #    if yale_var[2] in ALT_database[0]:
        #        ALT_database = ALT_database[0]
        #        frequency = (int(AAC) / int(AN)) * 100
        #        print str(yale_var[0])+"\t"+str(ALT_database) +"\t" + str(frequency) + "\t" + "ESP"
        #        output.write(str(yale_var[0])+"\t"+str(ALT_database) +"\t" + str(frequency) + "\t" + "ESP"+"\n")
        #        
        #if len(ALT_database) > 1:
        #    if yale_var[2] in ALT_database[1]:
        #        ALT_database = ALT_database[1]
        #        AAC = AAC.split(",")[1]
        #        frequency = (int(AAC) / int(AN)) * 100
        #        print str(yale_var[0])+"\t"+str(ALT_database) +"\t" + str(frequency) + "\t" + "ESP"
        #        output.write(str(yale_var[0])+"\t"+str(ALT_database) +"\t" + str(frequency) + "\t" + "ESP"+"\n")
        #        
        #    if yale_var[2] in ALT_database[0]:
        #        ALT_database = ALT_database[0]
        #        AAC = AAC.split(",")[0]
        #        frequency = (int(AAC) / int(AN)) * 100
        #        print str(yale_var[0])+"\t"+str(ALT_database) +"\t" + str(frequency) + "\t" + "ESP"
        #        output.write(str(yale_var[0])+"\t"+str(ALT_database) +"\t" + str(frequency) + "\t" + "ESP"+"\n")
