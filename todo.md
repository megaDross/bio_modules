## TODO

- write some code in get_seq to put the reverse complement of the scrapped seq int the mat CompareSanger.match_with_seq_file() method. Reverse complement the preseq then match then reverse complement the full_seq

- Process every match made in CompareSeqs.get_matching_seq_file() until an alternative base is found, instead of just processing the first one in the list. Also write the file used to compare to reference seq in the output file

- Fix the testing outcomes when --input_file and --output_file is use in get_seq.py

- remove transcription-translation feature. They are not useful, just some pointless bells and whistles.

- allow user to input there own transcript for scrapping exon information

- Use locally download genome instead of UCSC: have a script that downloads the genome and indexes it by then scrape from it using pyfasta, BioPythons SeqIO (records[chrom].seq[start:end]), bx-python, pysam or samtools (faidx hg19.fa 15:1566656:15768969)   https://www.biostars.org/p/1638/

- Create a configuration file/dict which has things like path to genome and seq_dir i.e.
      configuration = {"hg19": "path/to/hg19.fa", "hg38": "path/to/hg38.fa",
                       "seq_dir": "path/to/seq_dir/"}






The following AB1 files have given the expected output:
  FUK20-VG  A/G 15:48782270
  HX9-AD  G 15:48787732
  FUK31 C 17:48274593
  A01_HUK4_16_23DW_HUK4_ A/G 2:189863418
  B01_CXE_EM_CXE_F_003 W 17:48273298

The following AB1 files is giving 4 bases instead of the expected T/G. Good candidate for filtering bases by quality scores:
  15:48737567  H01_FUK15_RL_UKB_FUK15_F_015.ab1

