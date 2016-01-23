import os
import itertools
import re

op = open("UK_TAAD_Variants_in_dbSNP.vcf") # tabix <dbSNP.vcf> <positions
out = open("greater_than_0.1_AF.txt","w")
original = open("var.both.TAAD_UK_Nov2015.filters.vcf","r")
#
for i in op:
    info = i.split("\t")[7]
    loc = i.split("\t")[0]+":"+i.split("\t")[1]
    rs = i.split("\t")[2]
    search_CAF = re.findall(r"CAF=.*;",info)
    CAF = ','.join(search_CAF).split(",")
    if len(CAF) > 1:
            AF = re.sub(";","",CAF[1])
            RAF = CAF[0]
            if AF != "." and float(AF) > 0.1:
                #print loc+"\t"+rs+"\t"+str(AF) 
                out.write(loc+"\t"+rs+"\t"+str(AF)+"\t"+i.split("\t")[4]+"\n")
out.close()

out = open("greater_than_0.1_AF.txt","r")
oot = open("positions_in_TAAD_UK.txt","w")
positions_in_z = []



test = open("head.vcf","r")

#for x,z in itertools.product(original,out):
#    if not x.startswith("#") and not z.startswith("#"):
#        cells_x = x.split("\t")
#        pos_x = cells_x[0]+":"+cells_x[1]
#        alt_x = cells_x[4]
#        sample_cells_x = cells_x[16:len(cells_x)]
#        
#        cells_z = z.split("\t")
#        pos_z = cells_z[0]
#        alt_z = cells_z[3]
#        positions_in_z.append(pos_z+"\t"+alt_z)
#        
#        if pos_x+"\t"+alt_x in positions_in_z:
#            
#            if "./." in sample_cells_x:
#                print pos_x + "\t"+ "found in " + \
#                str(len(sample_cells_x) - sample_cells_x.count("./."))\
#                + " samples"
#            else:
#                print pos_x + "\t" + "nothing"
#     

dic = {}
# construct a dictionary where the keys are var pos and the item
# is the alt nucleotide. All of which is derived from the 
# greater_than_0.1_AF.txt
for i in out:
    if not i.startswith("#"):
        cells = i.split("\t")
        pos = cells[0]
        alt = cells[3]
        dic[pos]=alt

for (x),(pos,alt) in itertools.product(original,dic.items()):
    if not x.startswith("#"):
        cells_x = x.split("\t")
        pos_x = cells_x[0]+":"+cells_x[1]
        alt_x = cells_x[4]
        sample_cells_x = cells_x[16:len(cells_x)]
        sample_format = filter(lambda b: b != "./.", sample_cells_x)
        print len(sample_format)
        
        
        #if pos in pos_x and int(DP_x) > 50:
        #    if "./." in sample_cells_x:
        #        print pos_x + "\t"+ "found in " + \
        #        str(len(sample_cells_x) - sample_cells_x.count("./."))\
        #        + " samples"

        if pos in pos_x:
            if "./." in sample_cells_x:
                print pos_x + "\t"+ "found in " + \
                str(len(sample_cells_x) - sample_cells_x.count("./."))\
                + " samples"
        
#print dic.get("18:48584855")

# Look for samples that are of a certain depth only, look at genotype field 
# and VCF Specification list
#   - so only pull samples from the original if > 50 RD & Alt Genotype...
        
