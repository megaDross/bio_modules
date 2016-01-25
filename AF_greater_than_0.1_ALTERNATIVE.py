import csv, re

def filter_matching_alt_and_pos(first_file, second_file):
    # csv allows one to deal with rows only. nested loop means every row in file1
    # will be iterated along with every row in file2, allowing one to easyily filter.
    for x in csv.reader(open(first_file, 'rb'), delimiter='\t'):
        for y in csv.reader(open(second_file, 'rb'), delimiter='\t'):
            if x[3] == y[4] and x[0] in y[0]+":"+y[1]:
                z = filter(lambda a: a !="./." and "0/0" not in a, y)
                print z
                        


            

#def match_datestamp_and_alt_and_pos(first_file, second_file):
#    
#    for z in filter_matching_alt_and_pos(first_file, second_file):
#        z = filter(lambda a: a !="./.", z)
#        
                 #if not int(chunk[2])> 50:
                 #   # continue automatically skips to the next iteration on element
                 #   continue
                 #if not chunk[:1] == "1/1":
                 #   continue
                 #yield z

            
                
                


first_file = "greater_than_0.1_AF.txt"
second_file = "test.txt"
#print list(match_datestamp_and_alt_and_pos(first_file,second_file))
print filter_matching_alt_and_pos(first_file,second_file)