from Bio import SeqIO
import os, re, subprocess
import useful


file_path = useful.cwd_file_path(__file__)

class CompareSeqs(object):
    ''' A collection of methods used to compare a reference sequence 
        with a sanger sequence contained within a .seq file
    '''
    def __init__(self, upstream, downstream, seq_file=None, seq_dir=None):
        self.seq_file = open(seq_file).read().replace("\n","")
        self.upstream = upstream
        self.downstream = downstream
        self.seq_dir = seq_dir
    
    UIPAC = {"A":"A", "C":"C", "G":"G", "T":"T",
             "R":"A/G", "Y":"C/T", "S":"G/C",
             "W":"A/T", "K":"G/T", "M":"A/C",
             "B":"C/G/T", "D":"A/G/T", "H":"A/C/T",
             "V":"A/C/G", "N":"N"}

    @staticmethod
    def get_matching_seq_file(query, directory):
        ''' Find a file name that best matches given query and return the
            sequence string in matched file
        '''
        # appending is required in case there is more than one match, returns first match
        store_matches = []
        for f in os.listdir(directory):
            if query in f:
                file_match = directory+f
                store_matches.append(file_match)
                
        sorted_matches = sorted(store_matches)
        return sorted_matches[0]

    @staticmethod
    def convert_ab1_to_seq(ab1):
        ''' Extract the sequence string from .ab1 file and return

            WARNING: This method needs ALOT of testing
        '''
        handle = open(ab1, "rb")
        new_file_name = ab1.replace("ab1", "seq")
        new_file = open(new_file_name, "w")

        # Use SeqIO and re to get the sequence from the ab1 file
        for record in SeqIO.parse(handle, "abi"):
                raw_data = record.annotations.get("abif_raw")
                sequence = raw_data.get("PBAS2")
                new_file.write(sequence)
                new_file.close()
                return(new_file_name)


    def match_with_seq_file(self,sequence):
        ''' Find part of a given sequence in a given seq_file and return the equivalent 
            full sequence deriving from the .seq file

            returns a tuple containing the matched sanger sequence and the variant 
            position base and the index of the var_pos_seq
        '''
        if self.seq_file:
            preseq = sequence[:self.upstream].upper()
            ref_seq = sequence[self.upstream].upper()
            postseq = sequence[self.upstream:].upper()[1:]
            
            
            if re.search(preseq, self.seq_file):    
                start, end, matched_seq = CompareSeqs.get_start_end_indexes(preseq, 
                                                                            self.seq_file)
                var_pos_seq = CompareSeqs.UIPAC.get(self.seq_file[end])
                downstream_seq = self.seq_file[end+1:end+self.downstream+1]
                full_seq = "".join((matched_seq.lower(),var_pos_seq.upper(),
                                     downstream_seq.lower()))
                
                return (full_seq,ref_seq,var_pos_seq.upper(), end)

            elif re.search(postseq, self.seq_file):
                #print("POST")
                start, end, matched_seq = CompareSeqs.get_start_end_indexes(postseq, 
                                                                            self.seq_file)
                var_pos_seq = CompareSeqs.UIPAC.get(self.seq_file[start-1])
                # below may not work great if it produces a negative number for indexing
                # i.e if start index = 6, self.downstream = 20
                upstream_seq = self.seq_file[(start-1)-self.downstream:start-1]
                full_seq = "".join((upstream_seq.lower(),var_pos_seq.upper(),
                                    matched_seq.lower()))
                return (full_seq,ref_seq,var_pos_seq.upper(), start-1)
            
            
            else:
                pass
        else:
            pass

    @staticmethod
    def get_start_end_indexes(seq, seq_file):
        ''' Find a given string in a given file and get the indexes of said
            string in the file
        '''
        find = [(m.start(0), m.end(0)) for m in re.finditer(seq, seq_file)][0]
        start = find[0]
        end = find[1]
        matched_seq = seq_file[start:end]
        return (start, end, matched_seq)

    @staticmethod
    def compare_nucleotides(base_1, base_2):
            ''' Compare two nucleotides
            '''
            if base_1 != base_2:
                return("the nucleotides given are DIFFERENT",2)

            elif base_1== base_2:
                return("the nucleotides given are the SAME",1)
            
            else:
                return "not found"


    @staticmethod
    def get_het_call(ab1_file, var_index):
        ''' Get the two bases being called as N
        '''
        ttuner = "/home/david/bin/tracetuner_3.0.6beta/rel/Linux_64/ttuner"
        out_dir = filepath.split("/")[:-1]+"test/test_files/"
        subprocess.call([ttuner, "-tabd", out_dir, "-het", "-id", self.seq_dir])
        subprocess.call([ttuner, "-s", "-id", self.seq_dir])

        ttuner_seq_file = open(out_dir+ab1_file+".seq")
        tab_file = open(out_dir+ab1_file+".tab")
        split_tab = [x.rstrip("\n").split(" ") for x in tab_file 
                     if not x.startswith("#")]
      
        het_base1 = "".join([x.replace("\n", "") for x in ttuner_seq_file 
                         if not x.startswith(">")])
       
        for line in split_tab:
            het_base2, qual, scan_pos, index = (line[3], line[9], line[17], line[24])
        
        return het_base1+"/"+het_base2    




