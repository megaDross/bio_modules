### HOW IT WORKS

1 - a vcf file containing variants which need to be validated by Sanger sequencing is converted to a tsv file using the --vcf option. All PCR products submitted for Sanger sequencing should be named as they are in the first column of the converted tsv. This ensures the variant referenced in the tsv is correctly matched with the intended AB1 file. 

2 - Type of mutation (SNP, insertion or deletion) is assessed

3 - Each ID (column 1 entries) in the tsv file is used to pair a specific variant with the correct AB1 files. 

4 - AB1 files are converted to .tab, .seq & .poly files using ttuner. 

--- for each matched file

5 - Using the upstream, downstream and hg_version arguments, a portion of the reference sequence of interest is scrapped (via UCSC or locally).

6 - If an _R_, -R-, _R. or -R. is in the AB1 file name, then reverse complement the sequence string.

7 - Use the reference sequence preceding and proceding the variant to find the variants position (index number) within the AB1 file.  

8a - If we are looking for a SNP, use the index number to extract the string at the position in which one expects to find the variant. Extract the two most prominant nucleotides at this position from the generated poly file and concatenate them.

8b - If we are looking for an INDEL, get every het call within the sequence procedding the variant of interest (length of which specified by the downstream variable) by extracting said calls from the tab file. Check that the alternate calls are within the sequence extracted from the AB1 file. If so, return indel as being real. 

9 - Assess whether the het call, if found, is the same or different when compared to the reference base at said position. Return this for every matched file found.

---

10 - Filter out any of the returned values where no AB1 file was matched

11 - Determine which returned values contains the variant called within the vcf file. 

12 - Communicate to the user whether the variant was found or not.


