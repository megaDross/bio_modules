import os, sys, io

   
def write_to_output(seq_list, output_name, header=None):
    ''' creates an output file containing all scrapped data
    '''
    print("Writing all findings to {} ....".format(output_name))

    # create a list comprehension of all input and open output file
    seq_list = [x for x in seq_list]
    output = io.open(output_name,"w", encoding='utf-8')

    # if header parsed then write to output file
    if header:
        output.write(header)

    # write each element of list to the file and close
    for seq in seq_list:
        output.write(seq+"\n")
    output.close()
    return output


