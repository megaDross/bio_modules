# GeneaPy

A collection of scripts to help aid PCR-based genetic variant validations.

## Install
To install and test GeneaPy:
```bash
git clone https://github.com/superDross/GeneaPy
pip3 install -r GeneaPy/requirements.txt
export PYTHONPATH=$PYTHONPATH:/path/to/GeneaPy/
cd GeneaPy/test
python3 -m unittest *.py
```

## get_locus_metadata
Scrape a genomic positions metadata from Ensembl and UCSC.

#### Example
```
$ python3 get_locus_metadata.py --position chr15:48729400 --genome_version hg19
--------------------------------------------------
chr15:48729400 in hg19(75)
--------------------------------------------------
Gene
name: FBN1
id: ENSG00000166147
loc: 15:48700503-48938046
type: protein_coding

Transcript
name: FBN1-001
id: ENST00000316623
type: protein_coding
canon: True

Intron
id: N/A
no: 52/65
strand: -
--------------------------------------------------
>15:48729350-48729450
tttaatcaaaagattatctaataatgcaataatataattgctatctaaat
Gaagggacaaaaaagtagcacttaattttccaagatagatggagaaaaat
a
--------------------------------------------------
```

## get_seq
Scrapes a DNA sequence covering a given genomic range from the UCSC DAS server.

#### Example

Scrape sequence 10bp upstream and downstream from the given position from the UCSC DAS server
```
$ python3 get_seq.py chr1:169314424 --upstream 10 --downstream 10 --header --genome_version hg38
>1:169314414,169314434 hg38
tagtttctagTgtaaatacta
```
Scrape the genomic range's sequence from a local FASTA genome file
```
$ python3 get_seq.py chr1:169314424-169314623 --genome ~/human_genome19.fasta --header
>1:169314424,169314623 hg19
cattcctgattttgataatttgtatcttcagtcttttttcttggtcagtc
taaccaaagttttgcctatgttaattgtttccaataactaatttttagtt
tcattgattttctccattttctattttactgttttctactccaatactta
tctcctttcttttctatttcctttaggtttaatttgttcttatttttctt

```

## unknown_primer
Retrieve genetic metadata from a given pair of primers.

#### Example
From a parsed primer pair:
```
$ python3 unknown_primer.py --primers CTGTTCACAGGGCTTGTTCC CTGGGCAGAGAGTCATTTAAAGT --genome_version hg19
Primer  F_Primer        R_Primer        Genome  Gene    Transcript      Exon    Intron  Product_Size    Primer_Range    GC%
query   CTGTTCACAGGGCTTGTTCC    CTGGGCAGAGAGTCATTTAAAGT hg19    FBN1    FBN1-001        -       41/65   421bp   chr15:48755298-48755718   39.2
```
From a tab delimited input file (e.g. test/expected_output/unknown_primer_in.txt):
```
$ python3 unknown_primer.py --input input_file.txt --output primer_file.txt
```

## primer_finder
Takes variant position(s) as input and matches it with an appropriate primer in a given file containing primer information (primer database).

### Primer Database
A primer database can be created from an existing list of primer pairs provided it is formatted correctly (e.g. test/expected_output/unknown_primer_in.txt):
```
python3 unknown_primer.py --input primer_sequences.txt --output primer_database.txt
```

#### Example
Filter the primer database for primer pairs that are within FBN1 intron 21 human genome version 19 and which produce a 500bp product with a maximum of 60% GC content.
```
$ python3 primer_finder.py --database test/expected_output/primer_database.txt --gene FBN1 --intron 21 --size 500 --gc 60 --genome_version hg19
Primer              F_Primer              R_Primer Genome  Gene Transcript Exon Intron  Product_Size   GC%  Chrom     Start       End                         
HX10  TGCAAAGACCATTGGAGTGG  ATGGAGTCTGCAAGAACAGC   hg19  FBN1   FBN1-001    -  21/65           483  39.1     15  48787270  48787752
```
Filter for primers within a given position which is at least 100bp away from the forward and reverse primer.
```
$ python3 primer_finder.py --database test/expected_output/primer_database.txt --variant 15:48787370 --distance 100
Primer              F_Primer                R_Primer Genome  Gene Transcript   Exon Intron  Product_Size   GC%  Chrom     Start       End  Dist_F  Dist_R
FUK22  AGCTCAGAACCCCTAACCAC  AATGTCAGCTTTTCCTGCAAAA   hg19  FBN1   FBN1-001  22/66      -           582  40.4     15  48787031  48787612     339     242
```
Filter for primers within all positions detailed in a given file (e.g. test/primer_finder_input.txt)
```
$ python3 primer_finder.py --database ../test/expected_output/primer_database.txt --input  ../test/expected_output/primer_finder_input.txt --output filtered_primers.txt
```
