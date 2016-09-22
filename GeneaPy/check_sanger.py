import re
import GeneaPy.useful as useful
import get_AB1_file

# draw a map or mind map which visually shows how all this works. Perhaps read the OOP book to ensure best practices are adhered to.

# Perhaps have a function within this file which takes much of the actual basecalling out of the get_seq script (lines 187-214). life_saver()

# PLAYING WITH ACTUAL STRING: match_with_seq_file(), check_if_insertion(), check_if_deletion()

# HET CALLS: the rest

file_path = useful.cwd_file_path(__file__)

UIPAC = {"A":"A", "C":"C", "G":"G", "T":"T",
         "R":"A/G", "Y":"C/T", "S":"G/C",
         "W":"A/T", "K":"G/T", "M":"A/C",
         "B":"C/G/T", "D":"A/G/T", "H":"A/C/T",
         "V":"A/C/G", "N":"N"}

non_singular_bases = ["Y", "R", "W", "S", "K", "M"]



class CompareSeqs(object):
    ''' A collection of methods used to compare a reference sequence 
        with a sanger sequence contained within a .seq file
    '''
    def __init__(self, upstream, downstream, alt_answer=None, mut_type=None,
                 path_to_seq_file=None, seq_dir=None, seq_filename = None):
        self.seq_file = get_AB1_file.handle_seq_file(path_to_seq_file, seq_dir)
        # if statement helps with the recursive function in match_with_seq_file(), otherwise the actual reverse complemented sequence is used as the seq_filename, which causes errors in get_het_calls when it tries to open a string instead of an actual file
        self.seq_filename = path_to_seq_file if not seq_filename else seq_filename
        self.upstream = upstream
        self.downstream = downstream
        self.seq_dir = seq_dir
        self.alt_answer = alt_answer
        self.mut_type = mut_type

    def match_with_seq_file(self, sequence, num=2):
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
                
                # check if there is an insertion or deletion 
                if self.mut_type not in ("snp", "d"):
                    print("looking for insertion\n")
                    insertion = self.check_if_insertion(postseq, end) 
                else:
                    insertion = None

                if self.mut_type not in ("snp", "i"):
                    print("looking for deletion\n")
                    deletion = self.check_if_deletion(postseq, end) 
                else:
                    deletion = None

                #print(deletion)
                indel = insertion if insertion else deletion

                if indel:
                    print("Indel\n")
                    start_insert, end_insert, extra, indel_type = indel
                    postseq = sequence[self.upstream:].upper()[1+extra:]
                    downstream_seq = self.seq_file[end+1:end+self.downstream+1]
                    return (upstream_seq.lower(), downstream_seq.lower(),
                           ref_seq, "-", indel)
                

                # else check if it's a point mutation
                else:
                    print("looking for SNP\n")
                    var_pos_seq = UIPAC.get(self.seq_file[end])
                    downstream_seq = self.seq_file[end+1:end+self.downstream+1]
                    full_seq = "".join((upstream_seq.lower(),var_pos_seq.upper(),
                                         downstream_seq.lower()))
                    
                    return (upstream_seq.lower(), downstream_seq.lower(),
                            ref_seq, var_pos_seq.upper(), end)


            elif re.search(postseq, self.seq_file):
                start, end, downstream_seq = CompareSeqs.get_start_end_indexes(postseq, 
                                                                            self.seq_file)
                #deletion = self.check_if_deletion(sequence, preseq, postseq, start, end, "postseq")

                # for Reverse sequence files, get the actual index
                if num == 1:
                    start = len(self.seq_file) - start -1
                 

                var_pos_seq = UIPAC.get(self.seq_file[start-1])
                # below may not work great if it produces a negative number for indexing
                # i.e if start index = 6, self.downstream = 20
                upstream_seq = self.seq_file[(start-1)-self.downstream:start-1]
                full_seq = "".join((upstream_seq.lower(),var_pos_seq.upper(),
                                    downstream_seq.lower()))
                return (upstream_seq.lower(), downstream_seq.lower(), 
                        ref_seq,var_pos_seq.upper(), start-1)

             
            else:
                update_object = CompareSeqs(self.upstream, self.downstream, 
                                            self.alt_answer, self.mut_type, 
                                            useful.reverse_complement(self.seq_file), 
                                            self.seq_dir, self.seq_filename)

                return update_object.match_with_seq_file(sequence, num-1)

        else:
            pass


    def check_if_deletion(self, postseq, var_index, num=1, matched_seq_list=[]):
          ''' Determine whether an insertion occurs within a given variant position
          '''
          start_index_postseq = var_index + 1
  
          # stops recursive function if the number of bases is more than the postseq len
          if int(num) == self.upstream or \
             re.search(postseq, self.seq_file):
              return None

          # get het calls for all within the indexs that cover the postseq 
          seq_het_calls = self.index_basecall_dictionary(postseq,
                                                         start_index_postseq, num)
              
          matched_seq, seq_len = self.compare_ref_het_calls(postseq, seq_het_calls,
                                                      start_index_postseq, num, "deletion")

          # a means by which to get all matche_seqs despite the number of recursions
          matched_seq_list.append((len(matched_seq)/seq_len, num))
          
          # get tuple that has a matched_seq len closest to 1.0
          get_min = min(matched_seq_list, key=lambda x:abs(x[0]-1))
                                                           
          
          # TEMP
          #print("\nSEQ_FILE:\t"+ str(self.seq_file[start_index_postseq + 1:            start_index_postseq+self.upstream]))
          #print("-"*40+"\nLENGTH DIFF:\t"+str(len(matched_seq)/seq_len)+"\n"+"-"*40)
                  
          # if number of recursions is more than half
          if float(num) > float(len(postseq)/2)-1:  
              #print(get_min)
              num = get_min[1]
              return (var_index, var_index+num, num-1, "d")
          else:
              return self.check_if_deletion(postseq, var_index, num+1, matched_seq_list )

           
    def check_if_insertion(self, postseq, var_index, num=1):
        ''' Determine whether an insertion occurs within a given variant position
        '''
        start_index_postseq = var_index + 1

        # stops recursive function if the num gets too high
        if float(num) > float(len(postseq)/4):
            return None

        seq_het_calls = self.index_basecall_dictionary(postseq,
                                                       start_index_postseq, num)
        matched_seq, seq_len = self.compare_ref_het_calls(postseq, seq_het_calls,
                                                    start_index_postseq, num, "insertion")
        
        # ensure the the number of bases that match is at least 80% of the actual sequence used to compare with the hets, then return the base range at which the insertion covers (seq_len will decrease as the number of recursions increases).
        if len(matched_seq)/seq_len > 0.8:
            return (var_index, var_index+num, num-1, "i")
        else:
            return self.check_if_insertion(postseq, var_index, num+1)


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
        print(sorted_matches)
        return sorted_matches


    def base_caller(self, sorted_matches, ref_base):
        ''' use the sorted matches from get_het_call() to call the approiriate base
        ''' 
        if self.alt_answer:
            print("Alt answer:\t "+self.alt_answer)
            ppp = [UIPAC.get(x[1]) for x in sorted_matches]
            if self.alt_answer in ppp or self.alt_answer+"/"+ref_base in ppp\
               or ref_base+"/"+self.alt_answer in ppp:
                het_call =  ref_base+"/"+self.alt_answer
                print("Het Call:\t"+het_call)
            else:
                het_call = None

        elif not het_call:
            # if only one base, check if its a non singular base
            if len(sorted_matches) ==  1 and sorted_matches[0] not in non_singular_bases:
                het_call = "/".join((ref_base, sorted_matches[0][1]))
            
            elif len(sorted_matches) == 1:
                het_call = sorted_matches[0]
            
            # if more than one base 
            elif len(sorted_matches) > 1:
                first_call, second_call = (sorted_matches[0][1], sorted_matches[1][1])
                
                if first_call in non_singular_bases and \
                   second_call not in non_singular_bases and \
                   int(sorted_matches[0][0]) - int(sorted_matches[1][0]) > 10 and \
                   second_call in UIPAC.get(first_call):
                    het_call = UIPAC.get(first_call)

                # if more any calls have a non singular base within them
                elif [True for x in [first_call, second_call] if x in non_singular_bases]:
                    
                    # get rid of any no singular bases
                    clean_calls = [x for x in [first_call, second_call] if x not in 
                                   non_singular_bases]
                    
                    if len(clean_calls) == 1 and ref_base != clean_calls[0]:
                        het_call = "/".join((ref_base, clean_calls[0]))

                    elif not clean_calls and first_call in non_singular_bases:
                        het_call = UIPAC.get(first_call)

                    else: 
                        het_call = ref_base

                else:
                    het_call = "/".join((first_call,second_call))


        return het_call


    def index_basecall_dictionary(self, sequence, start_index_seq, num):
        ''' Create a dictioney where the key is the seq_file index of the 
            given sequence and the item is a tuple containing every called 
            base at said index
        '''
        seq_het_calls = {}

        # get the entire index (index from seq_file) range of the given sequence
        indexes_sequence = range(start_index_seq, start_index_seq+self.upstream)
        count = 0

        # get all base calls for every index in the given sequence & place in dict
        for index in indexes_sequence:
            hets = self.get_het_call(index)
            seq_het_calls[index] = ''

            if hets:
                convert = UIPAC.get(hets[0][1])
                if "/" in convert:
                    seq_het_calls[index] = (convert.split("/")[0], convert.split("/")[1])
                    
                elif hets[1][1] in non_singular_bases:
                    seq_het_calls[index] = (sequence[count], hets[0][1])

                else:
                    seq_het_calls[index] = (hets[0][1], hets[1][1])

            else:
                seq_het_calls[index] = tuple(sequence[count])
                
            count += 1

        return seq_het_calls

    
    def compare_ref_het_calls(self, sequence, het_call_dict, start_index_seq, num, indel):
        ''' Compare every base in a given sequence with the calls in the 
            dictionary generated by index_basecall_dictionary() starting 
            from the first index key in the dict plus the given num

            i.e. keys range from 91, 111 and num = 2 then:
                compare every sequence[0] with items in 93, then the
                sequence[1] with items in 94 etc.
        '''
        # entire index (from seq_file) range of a given sequence it starts from the first index plus num
        end_index_seq = start_index_seq + self.upstream

        if indel == "insertion":
            indexes_sequence = range(start_index_seq + num, end_index_seq)
            counts = range(0, self.upstream)

        elif indel == "deletion":
            indexes_sequence = range(start_index_seq + 1, end_index_seq)
            counts = range(0 + num - 1, self.upstream)
        
        #print("\n\nNUM:\t"+str(num))
        #print("\nREF_SEQ:\t"+str(sequence[0+num-1:self.upstream]))


        # the more recursion, the smaller this will become, hence why you can't just use self.downstream to divide the len(master_seq)
        len_indexes_sequence = 0
        matched_seq = []

        for count, index in zip(counts, indexes_sequence):
            len_indexes_sequence += 1
            if sequence[count] in het_call_dict.get(index):
                matched_seq.append(sequence[count])

        return (matched_seq, len_indexes_sequence-1) 

    

