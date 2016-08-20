from Bio import SeqIO
import os, re, subprocess
import GeneaPy.useful as useful


file_path = useful.cwd_file_path(__file__)

class CompareSeqs(object):
    ''' A collection of methods used to compare a reference sequence 
        with a sanger sequence contained within a .seq file
    '''
    def __init__(self, upstream, downstream, seq_file=None, seq_dir=None):
        if seq_file.endswith(".seq"):
            self.seq_file = open(seq_file).read().replace("\n","")
        else:
            self.seq_file = seq_file
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
        for f in [x for x in os.listdir(directory) if x.endswith(".ab1")]:
            if query.replace("-","_") in f or query.replace("_","-") in f:
                file_match = directory+f
                store_matches.append(file_match)
        
        # trying to find matches that have both _ and -
        sorted_matches = sorted(store_matches)
        return sorted_matches

    @staticmethod
    def convert_ab1_to_seq(ab1):
        ''' Extract the sequence string from .ab1 file and return
            
            INCOMPATIBLE WITH TTUNER!!!!!!!
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



    def match_with_seq_file(self,sequence, num=2):
        ''' Find part of a given sequence in a given seq_file and return the equivalent 

            returns a tuple containing the matched sanger sequence and the variant 
            position base and the index of the var_pos_seq
        '''
        if int(num) < 1:
            return None
        
        elif self.seq_file:
            if self.seq_file.endswith("ab1"):
                ttuner = "/home/david/bin/tracetuner_3.0.6beta/rel/Linux_64/ttuner"
                #subprocess.call([ttuner, "-sd", self.seq_dir[:-1], "-id", self.seq_dir[:-1]])
                seq_file = "".join([x.rstrip("\n") for x in open(self.seq_file+".seq")
                                            if not x.startswith(">")])
            else:
                seq_file = self.seq_file
            
            
            preseq = sequence[:self.upstream].upper()
            ref_seq = sequence[self.upstream].upper()
            postseq = sequence[self.upstream:].upper()[1:]
            
            
            if re.search(preseq, seq_file):    
                start, end, upstream_seq = CompareSeqs.get_start_end_indexes(preseq, 
                                                                            seq_file)
                # for Reverse sequence seq files, get the actual index
                if num == 1:
                    end = len(seq_file) - end -1

                var_pos_seq = CompareSeqs.UIPAC.get(seq_file[end])
                downstream_seq = seq_file[end+1:end+self.downstream+1]
                full_seq = "".join((upstream_seq.lower(),var_pos_seq.upper(),
                                     downstream_seq.lower()))
                
                return (upstream_seq.lower(), downstream_seq.lower(),
                        ref_seq, var_pos_seq.upper(), end)


            elif re.search(postseq, seq_file):
                start, end, downstream_seq = CompareSeqs.get_start_end_indexes(postseq, 
                                                                            seq_file)
                # for Reverse sequence files, get the actual index
                if num == 1:
                    start = len(seq_file) - start -1

                var_pos_seq = CompareSeqs.UIPAC.get(seq_file[start-1])
                # below may not work great if it produces a negative number for indexing
                # i.e if start index = 6, self.downstream = 20
                upstream_seq = seq_file[(start-1)-self.downstream:start-1]
                full_seq = "".join((upstream_seq.lower(),var_pos_seq.upper(),
                                    downstream_seq.lower()))
                return (upstream_seq.lower(), downstream_seq.lower(), 
                        ref_seq,var_pos_seq.upper(), start-1)

             
            else:
                update_object = CompareSeqs(self.upstream, self.downstream, 
                                            useful.reverse_complement(seq_file), 
                                            self.seq_dir)
                return update_object.match_with_seq_file(sequence, num-1)

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


    def get_het_call(self, var_index):
        ''' Get the two bases being called at a given variant position
        '''
        if self.seq_file.endswith(".ab1"):
            ttuner = "/home/david/bin/tracetuner_3.0.6beta/rel/Linux_64/ttuner"

            # if no tab or seq files found then create them
            if not os.path.isfile(self.seq_file+".tab"):
                subprocess.call([ttuner, "-tabd", self.seq_dir, "-id", 
                                 self.seq_dir, "-mix"])

            if not os.path.isfile(self.seq_file+".seq"):
                subprocess.call([ttuner, "-sd", self.seq_dir, "-id", self.seq_dir])

            # open tab file and split by space and filter out empty strings
            tab_file = open(self.seq_file+".tab")
            split_tab = [x.rstrip("\n").split(" ") for x in tab_file 
                         if not x.startswith("#")] 
            filtered_tab = useful.filter_list_of_lists(split_tab)

            # if an index in the tab file matches the variants index, then append the associated base to an empty list
            all_matches = []
            for line in filtered_tab:
                het_base2, qual, scan_pos, index = (line[0], line[1], line[2], line[3])
                if int(index) == int(var_index):
                    all_matches.append((qual, het_base2))

            # sort matches by quality and return highest bases as het call if more than one call is found for the given position/index
            sorted_matches = sorted(all_matches) 
            
            return sorted_matches

    
    def base_caller(self, sorted_matches, ref_base):
        ''' use the sorted matches from get_het_call() to call the approiriate base
        '''
        non_singular_bases = ["Y", "R", "W", "S", "K", "M"]
        
        # if only one base, check if its a non singular base
        if len(sorted_matches) ==  1 and sorted_matches[0] not in non_singular_bases:
            het_call = "/".join((ref_base, sorted_matches[0][1]))
        
        elif len(sorted_matches) == 1:
            het_call = sorted_matches[0]
        
        # if more than one base 
        elif len(sorted_matches) > 1:
            first_call, second_call = (sorted_matches[0][1], sorted_matches[1][1])
            
            # if more any calls have a non singular base within them
            if [True for x in [first_call, second_call] if x in non_singular_bases]:
                
                clean_calls = [x for x in [first_call, second_call] if x not in 
                               non_singular_bases]
                
                if len(clean_calls) == 1 and ref_base != clean_calls[0]:
                    het_call = "/".join((ref_base, clean_calls[0]))
                else: 
                    het_call = ref_base

            else:
                het_call = "/".join((first_call,second_call))

        else:
            het_call = None

        return het_call




