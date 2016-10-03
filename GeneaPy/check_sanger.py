import re, regex
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
    def __init__(self, upstream, downstream, ref_base=None,  alt_answer=None,
                 mut_type=None, path_to_seq_file=None, seq_dir=None, 
                 seq_filename = None):
        self.seq_file = get_AB1_file.handle_seq_file(path_to_seq_file, seq_dir)
        # if statement helps with the recursive function in match_with_seq_file(), otherwise the actual reverse complemented sequence is used as the seq_filename, which causes errors in get_het_calls when it tries to open a string instead of an actual file
        self.seq_filename = path_to_seq_file if not seq_filename else seq_filename
        self.upstream = upstream
        self.downstream = downstream
        self.seq_dir = seq_dir
        self.mut_type = mut_type
        self.alt_answer = alt_answer
        self.ref_base = ref_base

    def match_with_seq_file(self, sequence, num=2):
        ''' Find part of a given sequence in a given seq_file and return the equivalent 

            returns a tuple containing the matched sanger sequence and the variant 
            position base and the index of the var_pos_seq
        '''
        print(num)

        if int(num) < 1:
            return None
        
        elif self.seq_file:
            
            preseq = sequence[:self.upstream].upper()
            ref_seq = sequence[self.upstream].upper()
            postseq = sequence[self.upstream:].upper()[1:]

            print(preseq)
            print(postseq)
            print(self.seq_filename) 
            
            print(self.seq_file)
            #if re.search(preseq, self.seq_file):    
            if regex.search(r'(?:'+preseq+'){s<=2}', self.seq_file):    
                print("\nPRESEQ MATCH\n")

                start, end, upstream_seq = CompareSeqs.get_start_end_indexes(preseq, 
                                                                            self.seq_file)

                # for Reverse sequence seq files, get the actual index
                if num == 1:
                    end = len(self.seq_file) - end -1
                
                # check if there is an insertion or deletion 
                if self.mut_type != "snp":
                    print("looking for indel\n")
                    indel = self.check_if_indel(postseq, end) 
                else:
                    indel = None
                               

                if indel:
                    downstream_seq = postseq[len(self.ref_base)-1:]
                    return (upstream_seq.lower(), downstream_seq.lower(),
                            self.ref_base, "-", indel)
                           

                # else check if it's a point mutation
                if self.mut_type not in ("d", "i"):
                    print("looking for SNP\n")
                    var_pos_seq = UIPAC.get(self.seq_file[end])
                    downstream_seq = self.seq_file[end+1:end+self.downstream+1]
                    full_seq = "".join((upstream_seq.lower(),var_pos_seq.upper(),
                                         downstream_seq.lower()))
                    
                    return (upstream_seq.lower(), downstream_seq.lower(),
                            ref_seq, var_pos_seq.upper(), end)


            # elif re.search(postseq, self.seq_file):
            elif regex.search(r'(?:'+postseq+'){s<=2}', self.seq_file):
                print("\nPOSTSEQ MATCH\n")

                start, end, downstream_seq = CompareSeqs.get_start_end_indexes(postseq, 
                                                                            self.seq_file)

                # for Reverse sequence files, get the actual index
                if num == 1:
                    start = len(self.seq_file) - start -1
               
                # check if there is an insertion or deletion 
                if self.mut_type != "snp":
                    print("looking for indel\n")
                    indel = self.check_if_indel(postseq, end) 
                else:
                    indel = None
  

                if indel:
                    upstream_seq = postseq[:len(self.ref_base)-1]
                    return (upstream_seq.lower(), downstream_seq.lower(),
                            self.ref_base, "-", indel)
                           

                # else check if it's a point mutation
                if self.mut_type not in ("d", "i"):
                    print("looking for SNP\n")
                     
                    var_pos_seq = UIPAC.get(self.seq_file[start-1])
                    # below may not work great if it produces a negative number for indexing
                    # i.e if start index = 6, self.downstream = 20
                    upstream_seq = self.seq_file[(start-1)-self.downstream:start-1]
                    full_seq = "".join((upstream_seq.lower(),var_pos_seq.upper(),
                                        downstream_seq.lower()))
                    
                    return (upstream_seq.lower(), downstream_seq.lower(), 
                            ref_seq,var_pos_seq.upper(), start-1)

             
            else:
                # reverse complement the sequence and use as input recursively
                print("\nREVERSE\n")
                update_object = CompareSeqs(self.upstream, self.downstream, 
                                            self.ref_base, 
                                            self.alt_answer, self.mut_type, 
                                            useful.reverse_complement(self.seq_file), 
                                            self.seq_dir, self.seq_filename)

                return update_object.match_with_seq_file(sequence, num-1)

        else:
            pass


    def check_if_indel(self, postseq, var_index, num=1, matched_seq_list=[]):
        ''' Determine whether a deletion occurs within a given variant position
        '''
        if self.ref_base and self.alt_answer:
            start_index_postseq = var_index + 1

            # stops recursive function if the number of bases is more than the postseq len
            if int(num) == self.upstream or \
             re.search(postseq, self.seq_file):
                return None

            # get every het call for each position within the postseq 
            seq_het_calls = self.index_basecall_dictionary(postseq,
                                                           start_index_postseq, num)
            
            # use the number of bases to construct a sequence that depicts the resulting string produced from the deletion (al_answer_deletion)
            num_bases = len(self.ref_base) - 1
            alt_answer_deletion = self.alt_answer + postseq[num_bases:num_bases+num_bases]
            
            # check that the deleted base is within the positions het calls dict, if so add 1 to score and return the het call if all are found in the dict
            score = 0
            for answer, key_item  in zip(alt_answer_deletion[1:], 
                                         sorted(seq_het_calls.items())[:num_bases]):
                pos, calls = key_item
                print((pos, calls, answer))
                if answer in calls:
                    score += 1

            if score == len(self.ref_base) - 1:
                return self.ref_base+"/"+self.alt_answer

        
                 
    
    @staticmethod
    def get_start_end_indexes(seq, seq_file):
        ''' Find a given string in a given file and get the indexes of said
            string in the file
        '''
        find = [(m.start(0), m.end(0)) for m in regex.finditer(r'(?:'+seq+'){s<=2}'
                                                               , seq_file)][0]
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
        print("SORTED "+str(sorted_matches))
        return sorted_matches


    def base_caller(self, sorted_matches, ref_base):
        ''' use the sorted matches from get_het_call() to call the approiriate base
        ''' 

        if self.alt_answer:
            print("Alt answer:\t "+self.alt_answer)
            
            # complement bases
            comp_alt_answer = "".join([useful.complement.get(x) for x in self.alt_answer])
            comp_ref_base = useful.complement.get(ref_base)
            
            # get all het calls in a list
            ppp = [UIPAC.get(x[1]) for x in sorted_matches]

            if self.alt_answer in ppp or self.alt_answer+"/"+ref_base in ppp\
               or ref_base+"/"+self.alt_answer in ppp:
                het_call = ref_base + "/" + self.alt_answer
                print("Het Call:\t"+het_call)

            # check complements in case seq_file is reverse complemented
            elif comp_alt_answer in ppp or comp_alt_answer+"/"+comp_ref_base in ppp \
               or comp_ref_base+"/"+comp_alt_answer in ppp:
                het_call = ref_base + "/" + self.alt_answer
                print("Comp Het Call:\t"+het_call)

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
            given base in the sequence and the item is a tuple containing every called 
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

        print("HET CALL DICT\t"+str(seq_het_calls))
        return seq_het_calls

    
   
