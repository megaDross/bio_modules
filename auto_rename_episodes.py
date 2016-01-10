import os
import re
from glob import glob




for Pre_Title in glob('*/'):
    Almost_Title = re.sub("/","",Pre_Title)
    Title = re.sub(" ","_",Almost_Title)
    #print Title

for Season_no in glob('*/*/'):
    pre = Season_no.split("/")[1]
    Season = re.split(r"[ _]", pre)[1]
    #print Season

first_episode = 1

