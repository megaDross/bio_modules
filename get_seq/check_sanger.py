import os, re


class CompareSeqs(object):
    def __init__(self, upstream, downstream, seq_file):
        self.seq_file = seq_file
        self.upstream = upstream
        self.downstream = downstream
    
    UIPAC = {"A":"A", "C":"C", "G":"G", "T":"T",
             "R":"A/G", "Y":"C/T", "S":"G/C",
             "W":"A/T", "K":"G/T", "M":"A/C",
             "B":"C/G/T", "D":"A/G/T", "H":"A/C/T",
             "V":"A/C/G", "N":"N"}

    @staticmethod
    def get_matching_seq_file(query, directory):
        ''' find a file name that best matches given query
        '''
        # appending is required in case there is more than one match, returns first match
        store_matches = []
        for f in os.listdir(directory):
            if query in f:
                file_match = directory+f
                store_matches.append(file_match)
                
         
        sorted_matches = sorted(store_matches)
        return sorted_matches[0]

    def match_with_seq_file(self,sequence):
        ''' search for the sequence output from 
            get_region_info() in a given 
            .seq file and output it

            returns a tuple containing the sanger
            sequence and the var_pos nucelotide

            # NEEDS MUCH MORE TESTING
        '''
        if self.seq_file:
            # get the sequence preceding the var_pos (preseq) and the var_pos sequence (ref_seq) 
            # from the returned get_region_info() value 
            preseq = sequence[:self.upstream].upper()
            ref_seq = sequence[self.upstream].upper()
            seq_file = open(self.seq_file, "r").read()
            seq_file = seq_file.replace("\n","")
            
            # find the preseq in the seq_file string and output the indexes where the match occurred within the 
            # seq_file as a tuple
            if re.search(preseq, seq_file):
                find = [(m.start(0), m.end(0)) for m in re.finditer(preseq, seq_file)][0]
                start = find[0]
                end = find[1]
                
                # get the full sequence of interest from the seq_file
                matched_seq = seq_file[start:end]
                var_pos_seq = CompareSeqs.UIPAC.get(seq_file[end])  # convert the UIPAC to bases
                downstream_seq = seq_file[end+1:end+self.downstream+1]
                full_seq = "".join((matched_seq.lower(),var_pos_seq.upper(),
                                     downstream_seq.lower()))
                return(full_seq,ref_seq,var_pos_seq.upper())
            
            else:
                pass
        else:
            pass


    @staticmethod
    def compare_nucleotides(base_1, base_2):
            '''compare two nucleotides
            '''
            
            # assess whether a variant is present in the sanger sequence in the given proposed variant position
            if base_1 != base_2:
                return("the nucleotides given are DIFFERENT",2)

            elif base_1== base_2:
                return("the nucleotides given are the SAME",1)
            
            else:
                return "not found"


