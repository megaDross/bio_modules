from Bio import SeqIO
import os, re, subprocess
import itertools
import GeneaPy.useful as useful


file_path = useful.cwd_file_path(__file__)

class CompareSeqs(object):
    ''' A collection of methods used to compare a reference sequence 
        with a sanger sequence contained within a .seq file
    '''
    def __init__(self, upstream, downstream, path_to_seq_file=None, seq_dir=None,
                 seq_filename = None):
        self.seq_file = CompareSeqs.handle_seq_file(path_to_seq_file, seq_dir)
        # if statement helps with the recursive function in match_with_seq_file(), otherwise the actual reverse complemented sequence is used as the seq_filename, which causes errors in get_het_calls when it tries to open a string instead of an actual file
        self.seq_filename = path_to_seq_file if not seq_filename else seq_filename
        self.upstream = upstream
        self.downstream = downstream
        self.seq_dir = seq_dir
    
    UIPAC = {"A":"A", "C":"C", "G":"G", "T":"T",
             "R":"A/G", "Y":"C/T", "S":"G/C",
             "W":"A/T", "K":"G/T", "M":"A/C",
             "B":"C/G/T", "D":"A/G/T", "H":"A/C/T",
             "V":"A/C/G", "N":"N"}

    non_singular_bases = ["Y", "R", "W", "S", "K", "M"]

    @staticmethod
    def handle_seq_file(seq_file, seq_dir):
        ''' Decide how to process seq_files
        
            -id is used instead of -if as the -if only works on a file containing paths to the ab1 files requiring processing; cannot directly select a file in ttuner.
        '''
        ttuner = "/home/david/bin/tracetuner_3.0.6beta/rel/Linux_64/ttuner"

        if seq_file.endswith(".seq"):
            seq_file = "".join([x.strip("\n") for x in open(seq_file) 
                                              if not x.startswith(">")])

        elif seq_file.endswith("ab1") and os.path.isfile(seq_file+".seq") \
            and os.path.isfile(seq_file+".tab"):
            seq_file = "".join([x.rstrip("\n") for x in open(seq_file+".seq")
                                            if not x.startswith(">")]) 
 
        elif seq_file.endswith("ab1") and not os.path.isfile(seq_file+".tab"):
            subprocess.call([ttuner, "-tabd", seq_dir, "-id", 
                             seq_dir, "-mix"], stdout=open(os.devnull, 'wb'),
                             stderr=open(os.devnull, 'wb'))
            subprocess.call([ttuner, "-sd", seq_dir, "-id", seq_dir], 
                            stdout=open(os.devnull, 'wb'), stderr=open(os.devnull, 'wb'))

        else:
            seq_file = seq_file

        return seq_file

    @staticmethod
    def get_matching_seq_file(query, directory):
        ''' Find a file name that best matches given query and return the
            sequence string in matched file
        '''
        length = len(query)
        cut = query.replace("_","-").split("-")
        reversed_query = "_".join((cut[1],cut[0])) if len(cut) > 1 else query
        

        if length < 6:
            return []

        # appending is required in case there is more than one match, returns first match
        store_matches = []
        for f in [x for x in os.listdir(directory) if x.endswith(".ab1")]:
            if query.replace("-","_") in f or query.replace("_","-") in f or \
               reversed_query in f:
                file_match = directory+f
                store_matches.append(file_match)
        
        # trying to find matches that have both _ and -
        sorted_matches = sorted(store_matches)

        if not sorted_matches:
            return CompareSeqs.get_matching_seq_file(query[:-1], directory)

        print(sorted_matches)
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
            
            preseq = sequence[:self.upstream].upper()
            ref_seq = sequence[self.upstream].upper()
            postseq = sequence[self.upstream:].upper()[1:]
            

            if re.search(preseq, self.seq_file):    
                start, end, upstream_seq = CompareSeqs.get_start_end_indexes(preseq, 
                                                                            self.seq_file)
                # for Reverse sequence seq files, get the actual index
                if num == 1:
                    end = len(self.seq_file) - end -1
                
                # check if there is an insertion
                insertion = self.check_if_insertion(postseq, end)
                if insertion:
                    start_insert, end_insert, extra = insertion
                    postseq = sequence[self.upstream:].upper()[1+extra:]
                    downstream_seq = self.seq_file[end+1:end+self.downstream+1]
                    return (upstream_seq.lower(), downstream_seq.lower(),
                           ref_seq, "-", insertion)
                

                # else check if it's a point mutation
                else:
                    var_pos_seq = CompareSeqs.UIPAC.get(self.seq_file[end])
                    downstream_seq = self.seq_file[end+1:end+self.downstream+1]
                    full_seq = "".join((upstream_seq.lower(),var_pos_seq.upper(),
                                         downstream_seq.lower()))
                    
                    return (upstream_seq.lower(), downstream_seq.lower(),
                            ref_seq, var_pos_seq.upper(), end)


            elif re.search(postseq, self.seq_file):
                start, end, downstream_seq = CompareSeqs.get_start_end_indexes(postseq, 
                                                                            self.seq_file)
                # for Reverse sequence files, get the actual index
                if num == 1:
                    start = len(self.seq_file) - start -1
            
                var_pos_seq = CompareSeqs.UIPAC.get(self.seq_file[start-1])
                # below may not work great if it produces a negative number for indexing
                # i.e if start index = 6, self.downstream = 20
                upstream_seq = self.seq_file[(start-1)-self.downstream:start-1]
                full_seq = "".join((upstream_seq.lower(),var_pos_seq.upper(),
                                    downstream_seq.lower()))
                return (upstream_seq.lower(), downstream_seq.lower(), 
                        ref_seq,var_pos_seq.upper(), start-1)

             
            else:
                update_object = CompareSeqs(self.upstream, self.downstream, 
                                            useful.reverse_complement(self.seq_file), 
                                            self.seq_dir, self.seq_filename)

                return update_object.match_with_seq_file(sequence, num-1)

        else:
            pass



    def check_if_insertion(self, postseq, var_index, num=1):
        ''' Determine whether an insertion occurs within a given variant position
        '''
        start_index_postseq = var_index + 1

        # stops recursive function if the num gets too high
        if float(num) > float(len(postseq)/4):
            return None

        seq_het_calls = self.index_basecall_dictionary(postseq,
                                                       start_index_postseq, num)
        matched_seq, seq_len = self.get_index_range(postseq, seq_het_calls,
                                                    start_index_postseq, num)
        
        # ensure the difference in length is greater than 80% then return the base range at which the insertion covers (seq_len will decrease as the number of recursions increases).
        if len(matched_seq)/seq_len > 0.8:
            return (var_index, var_index+num, num-1)
        else:
            return self.check_if_insertion(postseq, var_index, num+1)



    def index_basecall_dictionary(self, sequence, start_index_seq, num):
        ''' Create a dictioney where the key is the seq_file index of the 
            given sequence and the item is a tuple containing every called 
            base at said index
        '''
        seq_het_calls = {}

        # entire index (index from seq_file) range of the given sequence
        indexes_sequence = range(start_index_seq, start_index_seq+self.upstream)
        count = 0

        # get all base calls for every index in the given sequence & place in dict
        for index in indexes_sequence:
            hets = self.get_het_call(index)
            seq_het_calls[index] = ''

            if hets:
                convert = CompareSeqs.UIPAC.get(hets[0][1])
                if "/" in convert:
                    seq_het_calls[index] = (convert.split("/")[0], convert.split("/")[1])
                    
                elif hets[1][1] in CompareSeqs.non_singular_bases:
                    seq_het_calls[index] = (sequence[count], hets[0][1])

                else:
                    seq_het_calls[index] = (hets[0][1], hets[1][1])

            else:
                seq_het_calls[index] = tuple(sequence[count])
                
            count += 1

        return seq_het_calls

    
    def get_index_range(self, sequence, het_call_dict, start_index_seq, num):
        ''' Compare every base in a given sequence with the calls in the 
            dictionary generated by index_basecall_dictionary() starting 
            from the first index key in the dict plus the given num

            i.e. keys range from 91, 111 and num = 2 then:
                compare every sequence[0] with items in 93, then the
                sequence[1] with items in 94 etc.
        '''
        # entire index (from seq_file) range of a given sequence it starts from the first index plus num
        end_index_seq = start_index_seq + self.upstream
        indexes_sequence = range(start_index_seq + num, end_index_seq)
        counts = range(0, self.upstream)

        # the more recursion, the smaller this will become, hence why you can't just use self.downstream to divide the len(master_seq)
        len_indexes_sequence = 0
        matched_seq = []

        for count, index in zip(counts, indexes_sequence):
            len_indexes_sequence += 1
            if sequence[count] in het_call_dict.get(index):
                matched_seq.append(sequence[count])

        return (matched_seq, len_indexes_sequence-1) 


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
        # open tab file and split by space and filter out empty strings
        tab_file = open(self.seq_filename+".tab")
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
        sorted_matches = sorted(all_matches, key=lambda x: int(x[0]), reverse=True)

        return sorted_matches

    
    def base_caller(self, sorted_matches, ref_base):
        ''' use the sorted matches from get_het_call() to call the approiriate base
        '''
        # if only one base, check if its a non singular base
        if len(sorted_matches) ==  1 and sorted_matches[0] not in CompareSeqs.non_singular_bases:
            het_call = "/".join((ref_base, sorted_matches[0][1]))
        
        elif len(sorted_matches) == 1:
            het_call = sorted_matches[0]
        
        # if more than one base 
        elif len(sorted_matches) > 1:
            first_call, second_call = (sorted_matches[0][1], sorted_matches[1][1])
            
            if first_call in CompareSeqs.non_singular_bases and \
               second_call not in CompareSeqs.non_singular_bases and \
               int(sorted_matches[0][0]) - int(sorted_matches[1][0]) > 10 and \
               second_call in CompareSeqs.UIPAC.get(first_call):
                het_call = CompareSeqs.UIPAC.get(first_call)

            # if more any calls have a non singular base within them
            elif [True for x in [first_call, second_call] if x in CompareSeqs.non_singular_bases]:
                
                # get rid of any no singular bases
                clean_calls = [x for x in [first_call, second_call] if x not in 
                               CompareSeqs.non_singular_bases]
                
                if len(clean_calls) == 1 and ref_base != clean_calls[0]:
                    het_call = "/".join((ref_base, clean_calls[0]))

                elif not clean_calls and first_call in CompareSeqs.non_singular_bases:
                    het_call = CompareSeqs.UIPAC.get(first_call)

                else: 
                    het_call = ref_base

            else:
                het_call = "/".join((first_call,second_call))

        else:
            het_call = None

        return het_call




