from Bio import SeqIO
import os, re, subprocess
import itertools
import GeneaPy.useful as useful


file_path = useful.cwd_file_path(__file__)

class CompareSeqs(object):
    ''' A collection of methods used to compare a reference sequence 
        with a sanger sequence contained within a .seq file
    '''
    def __init__(self, upstream, downstream, seq_file=None, seq_dir=None):
        self.seq_file = CompareSeqs.handle_seq_file(seq_file, seq_dir)
        self.seq_filename = seq_file
        self.upstream = upstream
        self.downstream = downstream
        self.seq_dir = seq_dir
    
    UIPAC = {"A":"A", "C":"C", "G":"G", "T":"T",
             "R":"A/G", "Y":"C/T", "S":"G/C",
             "W":"A/T", "K":"G/T", "M":"A/C",
             "B":"C/G/T", "D":"A/G/T", "H":"A/C/T",
             "V":"A/C/G", "N":"N"}

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
        #print(reversed_query)
        

        if length < 6:
        #    print("nothing found for\t"+query)
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
        #    print(query+"         "+str(len(sorted_matches)))
            return CompareSeqs.get_matching_seq_file(query[:-1], directory)

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
                ## DELETION
                #dictionary = {}

                #non_singular_bases = ["Y", "R", "W", "S", "K", "M"]
                #
                #n = 0
                #for index in range(end+1, end+self.upstream):
                #    hets = self.get_het_call(index)
                #    dictionary[index] = ''
                #    if hets:
                #        uip = CompareSeqs.UIPAC.get(hets[0][1])
                #        if "/" in uip:
                #            dictionary[index] = (uip.split("/")[0], uip.split("/")[1])
                #            
                #        elif hets[1][1] in non_singular_bases:
                #            dictionary[index] = (postseq[n], hets[0][1])

                #        else:
                #            dictionary[index] = (hets[0][1], hets[1][1])

                #    else:
                #        dictionary[index] = tuple(postseq[n])
                #    n += 1

                #            
                #print(dictionary)


                # END OF DELETION #

                # for Reverse sequence seq files, get the actual index
                if num == 1:
                    end = len(self.seq_file) - end -1

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




