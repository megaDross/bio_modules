# GeneaPy

A set of python3 tools that allows one to manipulate sequences and PCR primers 

## get_seq
Basic usage will return a sequence from the human genome reference (default: hg19) from a given genomic position and given numbers upstream and downstream. A genomic range can be given instead, which renders the upstream and downstream options redundant. This reference sequence can then optionally be compared to a .seq sanger sequence file to determine whether a variant is within the given genomic position in said .seq file. Gene, transcript and exonic information can optionally be scrapped from ensembl for the said genomic position given.

### Example Usage
```bash
# return a reference sequence 20 bp upstream and downstream from a given position. Scrapped from the defaulted human genome version (hg19)
python3 get_seq.py chr1:169314424

# as above but also return a header containing sequence metadata, including gene information
python3 get_seq.py chr1:169314424 --header --ensembl

# as above but also compare the contents of a given .seq file to the scrapped reference to determine whether the given position within the .seq file contains a variant.
python3 get_seq.py chr1:169314424 --seq_file NME7.seq --header --ensembl

# return a reference sequence from a given genomic range and provide metadata in a header (without gne information) 
python3 get_seq.py chr1,169315000,169316550 --header

# return a reference sequence 100b upstream and 200bp downstream from a given position scrapped from human genome version 38
python3 get_seq.py 1:169314424 --upstream 100 --downstream 200 --hg_version hg38

# return a reference sequence for each psition in a given input file, generate genomic/transcript/exonic infomation and compare to an automatically paired .seq file from a given directory. Finally, output all scrapped data to file.
python3 get_seq.py in.txt --header --ensembl --seq_dir ../seq_files/ --output_file out.txt
```


## primer_finder
Takes variant position(s) as input and matches it with an appropriate primer in a given file containing primer information.

### Example Usage
```bash
# find if a primer pair is available for a single genomic position
python3 primer_finder.py 15:48729400 --primer_database primer_database.txt

# filter out primer pairs where the forward or reverse primer is less than 100bp from the given genmic position
python3 primer_finder.py 15:48729400 --primer_database primer_database.txt --distance 100

# filter out primer pairs where the amplcion generated is more than 300bp
python3 primer_finder.py 15:48729400 --primer_database primer_database.txt --size 300

# filter out primer pairs where the amplicon generated is more than 50% GC rich
python3 primer_finder.py 15:48729400 --primer_database primer_database.txt --gc 50

# find if primer pairs are available for multiple positions in a given input file and ouput to file
python3 primer_finder.py --input_file in.txt --primer_database primer_database.txt --output_file out.txt

```


## unknown_primer
Take primer sequences as input and output the amplicon sequence information

### Example Usage
```bash
# find the amplicon information for a given primer pair for human genome version 19
python3 unknown_primer.py --primers TAACAGATTGATGATGCATG CCCATGAGTGGCTCCTAAA

# return amplicon infomration for every primer pair in a given input file for human genome version 38 and output the results to a file
python3 unknown_primer.py --input_file in.txt --hg_version hg38 --output_file out.txt
```

## Testing
To perform testing execute the following command within the test directory:
```bash  
python3 -m unittest *.py
```
