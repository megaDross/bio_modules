import os, sys

   
def write_to_output(self,seq_list,header=None):
    ''' creates an output file containing all scrapped data
    '''
    # create a list comprehension of all input and open output file
    seq_list = [x for x in seq_list]
    output = open(self.output_file,"w")

    # if header parsed then write to output file
    if header:
        output.write(header)

    # write each element of list to the file and close
    for seq in seq_list:
        output.write(seq+"\n")
    output.close()
    return output


