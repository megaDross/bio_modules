from __future__ import division
import csv
import re
import itertools

# downloaded the ExAC, ESP & dbSNP vcfs. bgzip -c file.vcf > file.vcf.bgzip
# tabix -p vcf file.vcf.bgzip
# created a file that contains the hg19 gene regions (genes in X/Z assay)
# and used that as input for i: tabix file.vcf.gz $i > output_file.txt
# this filters only for the variants described in $i

# the above output_file.txt for ExAC and ESP was used as input for the below
# program

# 20/01/16 SPLIT THIS CONFUSING MESS UP INTO MODULES, ONE FOR FIND AF/CAF IN EACH
# DATABASE

# not tested on 

HOME = "C:\cygwin64\home\dross11\\"
DB = "C:\cygwin64\home\dross11\VCF\DB\\"

var_file = csv.reader(open\
(HOME+"All_Yale_&_UK_Variants_160106.txt","rb")\
,delimiter="\t")

ExAC_test = csv.reader(open(DB+"output_ExAC.txt","rb").readlines()[1:]\
,delimiter="\t")

ESP_test = open(DB+"output_ESP.txt","rb")

db_test = open(DB+"output_dbSNP_all.txt","rb")

output_ESP = open(DB+"output_mutation_frequency_ESP.txt","w")
output_ExAC = open(DB+"output_mutation_frequency_ExAC.txt","w")
output_dbSNP = open(DB+"output_mutation_frequency_dbSNP.txt","w")

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
for yale_var,db_var in itertools.product(csv.reader(simplified,delimiter="\t"),ESP_test):
    database_var= db_var.split("\t")
    if len(database_var) >7 and yale_var[3]+":"+yale_var[4] == database_var[0]+":"+database_var[1]:
        ALT_database = database_var[4].split(",")
        info = database_var[7]
        search_AAC = re.findall(r"TAC=.*;",info)
            # TAC = [alt allele count,ref allele count] in INFO fiel ESP* VCF files
            # AAC = alt allele count # RAC = ref allele count # AN = Total alleles
        AAC = str(search_AAC).split(";")[0][6:].split(",")[0]
        RAC = str(search_AAC).split(";")[0][6:].split(",")[1]
        AN = int(AAC) + int(RAC)
        
        if len(ALT_database) == 1:
            if yale_var[2] in ALT_database[0]:
                ALT_database = ALT_database[0]
                frequency = (int(AAC) / int(AN)) * 100
                print str(yale_var[0])+"\t"+database_var[0]+":"+database_var[1]+"\t"+str(ALT_database) +"\t" + str(frequency) + "\t" + "ESP"
                output_ESP.write(str(yale_var[0])+"\t"+database_var[0]+":"+database_var[1]+"\t"+str(ALT_database) +"\t" + str(frequency) + "\t" + "ESP"+"\n")
                
        if len(ALT_database) > 1:
            if yale_var[2] in ALT_database[1]:
                ALT_database = ALT_database[1]
                AAC = AAC.split(",")[1]
                frequency = (int(AAC) / int(AN)) * 100
                print str(yale_var[0])+"\t"+database_var[0]+":"+database_var[1]+"\t"+str(ALT_database) +"\t" + str(frequency) + "\t" + "ESP"
                output_ESP.write(str(yale_var[0])+"\t"+database_var[0]+":"+database_var[1]+"\t"+str(ALT_database) +"\t" + str(frequency) + "\t" + "ESP"+"\n")
                
            if yale_var[2] in ALT_database[0]:
                ALT_database = ALT_database[0]
                AAC = AAC.split(",")[0]
                frequency = (int(AAC) / int(AN)) * 100
                print str(yale_var[0])+"\t"+database_var[0]+":"+database_var[1]+"\t"+str(ALT_database) +"\t" + str(frequency) + "\t" + "ESP"
                output_ESP.write(str(yale_var[0])+"\t"+database_var[0]+":"+database_var[1]+"\t"+str(ALT_database) +"\t" + str(frequency) + "\t" + "ESP"+"\n")

#checking the ExAC database
for yale_var,database_var in itertools.product(csv.reader(simplified,delimiter="\t"),ExAC_test):
    if len(database_var) >7 and yale_var[3]+":"+yale_var[4] == database_var[0]+":"+database_var[1]:
        ALT_database = database_var[4].split(",")
        info = database_var[7]
        search_AAC = re.findall(r"AC=.*;",info)
        search_AN = re.findall(r"AN=.*;",info)
            # TAC = [alt allele count,ref allele count] in INFO fiel ESP* VCF files
            # AAC = alt allele count # RAC = ref allele count # AN = Total alleles
        AAC = str(search_AAC).split(";")[0][5:]
        AN = str(search_AN).split(";")[0][5:]
                
        if len(ALT_database) == 1:
            if yale_var[2] in ALT_database[0]:
                ALT_database = ALT_database[0]
                frequency = (int(AAC) / int(AN)) * 100
                print str(yale_var[0])+"\t"+database_var[0]+":"+database_var[1]+"\t"+str(ALT_database) +"\t" + str(frequency) + "\t" + "ExAC"
                output_ExAC.write(str(yale_var[0])+"\t"+database_var[0]+":"+database_var[1]+"\t"+str(ALT_database) +"\t" + str(frequency) + "\t" + "ExAC"+"\n")
                
        if len(ALT_database) > 1:
            if yale_var[2] in ALT_database[1]:
                ALT_database = ALT_database[1]
                AAC = AAC.split(",")[1]
                frequency = (int(AAC) / int(AN)) * 100
                print str(yale_var[0])+"\t"+database_var[0]+":"+database_var[1]+"\t"+str(ALT_database) +"\t" + str(frequency) + "\t" + "ExAC"
                output_ExAC.write(str(yale_var[0])+"\t"+database_var[0]+":"+database_var[1]+"\t"+str(ALT_database) +"\t" + str(frequency) + "\t" + "ExAC"+"\n")
                
            if yale_var[2] in ALT_database[0]:
                ALT_database = ALT_database[0]
                AAC = AAC.split(",")[0]
                frequency = (int(AAC) / int(AN)) * 100
                print str(yale_var[0])+"\t"+database_var[0]+":"+database_var[1]+"\t"+str(ALT_database) +"\t" + str(frequency) + "\t" + "ExAC"
                output_ExAC.write(str(yale_var[0])+"\t"+database_var[0]+":"+database_var[1]+"\t"+str(ALT_database) +"\t" + str(frequency) +"%"+ "\t" + "ExAC"+"\n")
                
#checking the dbSNP database
for yale_var,db_var in itertools.product(csv.reader(simplified,delimiter="\t"),db_test):
    database_var= db_var.split("\t")
    
    if len(database_var) >2 and yale_var[3]+":"+yale_var[4] == database_var[0]+":"+database_var[1]:
        ALT_database = database_var[4].split(",")
        info = database_var[7]
        search_CAF = re.findall(r"CAF=.*;",info)
            # TAC = [alt allele count,ref allele count] in INFO fiel ESP* VCF files
            # AAC = alt allele count # RAC = ref allele count # AN = Total alleles
        CAF = ','.join(search_CAF).split(",")
        if len(CAF) > 1:
            AF = re.sub(";","",CAF[1])
            RAF = CAF[0] 
            #print yale_var[0]+"\t"+yale_var[3]+":"+yale_var[4]+"\t"+AF+"%"
            
        
            if len(ALT_database) == 1:
                if yale_var[2] in ALT_database[0]:
                    ALT_database = ALT_database[0]
                    print str(yale_var[0])+"\t"+database_var[0]+":"+database_var[1]+"\t"+str(ALT_database) +"\t" + AF + "\t" + "dbSNP"
                    output_dbSNP.write(str(yale_var[0])+"\t"+database_var[0]+":"+database_var[1]+"\t"+str(ALT_database) +"\t" + AF + "\t" + "dbSNP"+"\n")
                    
            if len(ALT_database) > 1:
                if yale_var[2] in ALT_database[1]:
                    ALT_database = ALT_database[1]
                    AF = re.sub(";","",CAF[2])
                    print str(yale_var[0])+"\t"+database_var[0]+":"+database_var[1]+"\t"+str(ALT_database) +"\t" + AF + "\t" + "dbSNP"
                    output_dbSNP.write(str(yale_var[0])+"\t"+database_var[0]+":"+database_var[1]+"\t"+str(ALT_database) +"\t" + AF + "\t" + "dbSNP"+"\n")
                    
                if yale_var[2] in ALT_database[0]:
                    ALT_database = ALT_database[0]
                    AF = re.sub(";","",CAF[1])
                    print str(yale_var[0])+"\t"+database_var[0]+":"+database_var[1]+"\t"+str(ALT_database) +"\t" + AF + "\t" + "dbSNP"
                    output_dbSNP.write(str(yale_var[0])+"\t"+database_var[0]+":"+database_var[1]+"\t"+str(ALT_database) +"\t" + AF + "\t" + "dbSNP"+"\n")
        else:

            print str(yale_var[0])+"\t"+database_var[0]+":"+database_var[1]+"\t"+str(database_var[4]) +"\t" + "unknown_AF" + "\t" + "dbSNP"
            output_dbSNP.write(str(yale_var[0])+"\t"+database_var[0]+":"+database_var[1]+"\t"+str(database_var[4]) +"\t" + "unknown_AF" + "\t" + "dbSNP"+"\n")
output_ESP.close()
output_ExAC.close()
output_dbSNP.close()
