import os, sys

class ProcessIO(object):
    ''' Process input and output, specifically designed to work with Click
    '''

    def __init__(self,input_file,output_file=None):
        
        self.input_file = input_file
        self.output_file = output_file

    def process_input(self):      
        ''' Determine if input is a string, file or piped from another program
            and return it in the apropriate format..
        '''
        # allows one to pipe in an argument at the commandline, requires required=False
        # in @click.argument()
        if not self.input_file:
            input_file = input()
        
        # if input is a tuple and contains more than one element, create a list 
        # THIS NEEDS TO BE SIMPLFIED (USED BY UNKNWON_PRIMER & EXTRACTINFO)
        elif len(tuple(self.input_file)) > 1:  # was len(tuple(self.input_file[0])) > 1
            #input_file = ["query"+"\t"+self.input_file[0]+"\t"+self.input_file[1]]
            input_file = "query"
            for i in range(len(tuple(self.input_file))):
                input_file += '\t'
                input_file += self.input_file[i]
            return [input_file]

        # if input is a tuple and contains only one element
        elif len(tuple(self.input_file)) == 1:
            input_file = ["query"+"\t"+self.input_file[0]]
            return input_file

        # if input is not a file, create a list
        elif os.path.isfile(self.input_file) is False:
            input_file = ["query"+"\t"+self.input_file]
            return input_file 
        
        # if input is a file, open file
        elif os.path.isfile(self.input_file) is True:
            input_file = open(self.input_file,"r+")
            return input_file

        else:
            print("Invalid Input!")
           

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


