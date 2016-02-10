import urllib2
import re

### WHY ARE YOU DOING THIS ###
# having looked at my validation database for variants AB>0.3 and not validating,
# I wonder whether I can see whether they cluster around the same exons? this would
# indicate whether some genomic areas are harder to sequence than others.

### POTENTIAL PROBLEMS ###
# the variable 'pos' is essentially the coordinates for where the variant occurs along
# with a dash and number reresenting 101bp upstream. This extra 200bp is required to 
# produce a range, as a range is all that can be plugged into the REST API.
# This could be problematic later on!!!!

transcript = "ENST00000316623.5:c.8504dupC"
pos = "15:48703298-48703498"
#pos = "7:140424943-140624564"
name = "LX14-AI-UKB2"


transcript2 = transcript.split(":")[0]

# use the Ensembl code to output transcript page and isolate total number of exons
url =urllib2.urlopen("http://grch37.ensembl.org/Homo_sapiens/Transcript/Exons?t="+transcript2)
search = url.read()
split_url = re.search(r"transcript has.* exons",search).group().split(">")
total_exons = split_url[len(split_url)-1].split(" ")[0]


# prints the end position of the exon from the variant position range given
# planning on using as a search term for the exons_list variable output
exon_url = urllib2.urlopen("http://grch37.rest.ensembl.org/overlap/region/human/"+pos+\
"?feature=exon;content-type=application/json").read()
end_pos_exon = re.search(r"end.*[0-9],",exon_url).group().split(",")[0].split(":")[1]


# prints all exons associated with Ensembl Code, and split them into seperate elements
exons_list = urllib2.urlopen("http://grch37.rest.ensembl.org/overlap/id/"+\
transcript2[:len(transcript2)-2]+"?feature=exon;content-type=application/json;expand=1").read()
seperated_exons_list = re.split(r'[{}]',exons_list)

results = []

# within the seperated exons list, try and find an entry which contains the nucleotide 
# end position and the ENSEMBL transcript ID and place this into an empty list. The empty
# list converts the nested lists into one single list, where every nested list becomes 
# an element.
for x in range(1,len(seperated_exons_list)):
    split_all_exons = re.split(r'[{}]',exons_list)[x] # = seperated_exons_list[x]
    search_all_exons = re.findall(transcript2[:len(transcript2)-2]+".*"+end_pos_exon+".*",split_all_exons)
    if search_all_exons is not None:
        results.append(search_all_exons)

# this gets rid of the empty lists ([]) within the list variable named 'empty'
matches = [z for z in results if z]
print matches
# This extracts the exon ID from the matched exon 
ID = str(matches[0]).split(",")[6]
almost = ID.split(":")[1]
exon_ID = re.sub("\"","",almost)

# This extracts the exon number from the matched exon
exon_number = str(matches[0]).split(",")[9].split(":")[1]


print exon_ID +":"+"\t"+exon_number+"/"+total_exons



#print exons_list
#print re.search(end_pos_exon+".*}",exons_list).group()
