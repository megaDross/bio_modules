## TODO - TOP

- parse every single query from the TAAD spreadsheet to get_seq and check the results (base changes) are the same as those found in the spreadsheet. This is the only way to truly determine the accuracy of get_seq.py

- write tests for the newly implemented features in this branch

- tidy up code

- Write some code that gets gene information at a given position from REST ensembl directly instead of relying on the bloated pyensembl module

- LONG TERM: detect InDels and not just SNVs



## TODO - MINOR

- allow user to input there own transcript for scrapping exon information

- Create a configuration file/dict which has things like path to genome and seq_dir i.e. use configparser

- An option that allows you to use an annotated genome to get gene/transcript/exon info http://www.ensembl.org/info/data/ftp/index.html



## DONE

- Fix the testing outcomes when --input_file and --output_file is use in get_seq.py

- Use locally download genome instead of UCSC: have a script that downloads the genome and indexes it by then scrape from it using pyfasta, BioPythons SeqIO (records[chrom].seq[start:end]), bx-python, pysam or samtools (faidx hg19.fa 15:1566656:15768969)   https://www.biostars.org/p/1638/

- sort matching het calls at the given index in order of quality score (highest to lowest) and return the two highest scoring bases

- Process every match made in CompareSeqs.get_matching_seq_file() until an alternative base is found, instead of just processing the first one in the list. Also write the file used to compare to reference seq in the output file

- remove transcription-translation feature. They are not useful, just some pointless bells and whistles.

- write some code in get_seq to put the reverse complement of the scrapped seq int the mat CompareSanger.match_with_seq_file() method. Reverse complement the preseq then match then reverse complement the full_seq

- before returning het call check if K M etc. so it doesnt get called like K/A

- in output_file, place the name of the seq_file used to compare to the var_base

- a boolean option (--download) which will download human genome, ttuner etc.

- ensure altered program works with existing testing and debug according to errors rasied

## BUGS

- if reverse complement is used then the opposite bases will be called i.e we expect G/T but because we used _R file we get the reverse C/A. Ensure this is actually/ definitely happening and fix accordingly

## TESTING

The following AB1 files have given the expected output:
  B02_FUK20_VG_UKB_FUK20_F_004.ab1    A/G 15:48782270
  G03_HX9_AD_HX9_F_013.ab1            G 15:48787732
  D01_FUK31_LS_UKB_FUK31_F_007.ab1    C 17:48274593
  A01_HUK4_16_23DW_HUK4_              A/G 2:189863418
  B01_CXE_EM_CXE_F_003                W 17:48273298

The following AB1 files is giving 4 bases instead of the expected T/G. Good candidate for filtering bases by quality scores:
  H01_FUK15_RL_UKB_FUK15_F_015.ab1  A/G  15:48737567 

The following is a good match for multiple matching seq files:
  G03_HX9_AD_HX9_F_013.ab1   G 15:48787732

AB1 files where the _R.ab1 has the variant present but the _F.ab1 does not:
  C04_GXYUK2_RD_PP_GXYUK2_RD_R_006.ab1  G/C 2:189861900       post_seq
  HUK3-SC-R.ab1                         G/A    2:189861205

