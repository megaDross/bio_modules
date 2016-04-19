import requests,bs4,re

web = "http://genome.ucsc.edu/cgi-bin/das/hg19/dna?segment="

r = requests.get(web+"chr1:169314404,169314444")
r.raise_for_status()

url = bs4.BeautifulSoup(r.text,"html.parser")
refined_search = re.findall(r"[tacg{5}].*",url)
print(refined_search)
