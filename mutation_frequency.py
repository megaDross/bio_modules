from __future__ import division
import csv
import re
import itertools

HOME = "C:\cygwin64\home\dross11\\"
DB = "C:\cygwin64\home\dross11\VCF\DB\\"

var_file = csv.reader(open\
(HOME+"All_Yale_&_UK_Variants_160106.txt","r+")\
,delimiter="\t")

ExAC_test = csv.reader(open(DB+"Small_ExAC.r0.3.sites.vep.vcf","r+").readlines()[1:]\
,delimiter="\t")

ESP_test = csv.reader(open\
(DB+"ESP6500SI-V2-SSA137.GRCh38-liftover.chr15.snps_indels.vcf","r+"),delimiter="\t")


simplified = []

# reapproriate the file that contains all TAAD variants so that the
# variant alias, ref allele, alt allele, chromosome of the variant pos
# & position of the variant position are tab delimintated and stored in
# a list named simplified
for foo in var_file:
    #if re.findall(r":",foo[18]):
     #   POS = re.sub(r"chr","",foo[18])
      #  rearranged_file = foo[0]+"\t"+foo[25]+"\t"+foo[26]+"\t"+POS
       # simplified_yale.append(rearranged_file)
    CHR = foo[18].split(":")[0]
    pre_POS = foo[18].split(":")
    if len(pre_POS) == 2:
        POS = re.sub(r"chr","",pre_POS[1])
        rearranged_file = foo[0]+"\t"+foo[25]+"\t"+foo[26]+"\t"+CHR+"\t"+POS
        simplified.append(rearranged_file)

# Checking the ESP databases
for yale_var,database_var in itertools.product(csv.reader(simplified,delimiter="\t"),ESP_test):
    if len(database_var) >7:
        # if the variant locationmatch and the alternative alleles match...
        # ONLY DOES EXACT MATCHES, SO IF TWO ALT THEN WONT WORK - FIX THIS
        if yale_var[3]+":"+yale_var[4] == database_var[0]+":"+database_var[1]\
        and yale_var[2] == database_var[4]:
            info = database_var[7]
            search = re.findall(r"TAC=.*;",info)
            # TAC = [alt allele count,ref allele count] in INFO fiel ESP* VCF files
            # AAC = alt allele count # RAC = ref allele count # AN = Total alleles
            AAC = str(search).split(";")[0][6:].split(",")[0]
            RAC = str(search).split(";")[0][6:].split(",")[1]
            AN = int(AAC) + int(RAC)
            frequency = (int(AAC) / AN) * 100
            #print yale_var[0]+"\t"+database_var[4]+"\t"+str(frequency)+"%"


#checking the ExAC database
for yale_var,database_var in itertools.product(csv.reader(simplified,delimiter="\t"),ExAC_test):
    if len(database_var) >7:
        # if the variant locationmatch and the alternative alleles match...
        # ONLY DOES EXACT MATCHES, SO IF TWO ALT THEN WONT WORK - FIX THIS
        if yale_var[3]+":"+yale_var[4] == database_var[0]+":"+database_var[1]\
        and yale_var[2] == database_var[4]:
            info = database_var[7]
            search = re.findall(r"AC=.*;",info)
            search2 = re.findall(r"AN=.*;",info)
            # TAC = [alt allele count,ref allele count] in INFO fiel ESP* VCF files
            # AAC = alt allele count # RAC = ref allele count # AN = Total alleles
            AAC = str(search).split(";")[0][5:]
            AN = str(search2).split(";")[0][5:]
            frequency = (int(AAC) / int(AN)) * 100
            print yale_var[0]+"\t"+database_var[4]+"\t"+str(frequency)+"%"
 

    
    
    
    
    
    
    
