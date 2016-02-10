import requests, json

transcript_and_hgvs = "ENST00000316623.5:c.8504dupC"
transcript = transcript_and_hgvs.split(":")[0][:-2]
#pos = "15:48703298-48703498"
pos = "7:140424943-140624564"
name = "LX14-AI-UKB2"

url = "http://grch37.rest.ensembl.org/overlap/id/"
ext = "?feature=exon;content-type=application/json;expand=1"
req = requests.get(url+transcript+ext)
req.raise_for_status()

all_data = json.loads(req.text)

# get the exon id and the start and end positions in which they occupy
exon_region = set()
for i in all_data:
    if i.get("Parent") == transcript:
        start = i.get("start")
        end = i.get("end")
        exon = i.get("exon_id")
        exon_start_end = (exon,start,end)
        exon_region.add(exon_start_end)
        # i.get("rank"), last entry is last exon of transcript
        
        
### NEXT TIME ON GET EXON NUMBERING... 
# 1 - use exon_region in a loop to match the variant position with an exon_id
# 2 - use this exon_id to get from all_data and extract rank whic equals exon number



