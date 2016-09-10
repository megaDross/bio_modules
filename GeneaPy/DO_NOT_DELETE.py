vcf = "/home/david/Downloads/FB_0031_GTGAAA_L001_pe.vcf"


file_of_interest = ["test1\t15:17383686\tA/G", 
                    "test2\t17:18975678\tACGT/T",
                    "test3\t20:17896459\tG/TGAT"]

def vcf2validationlist(vcf, filename, starting_num=1):
    ''' Create a variant validation tsv file that can be used
        as input to get_seq.py
    '''
    counter = starting_num  
    out_file = open(filename, "w")

    for line in open(vcf, "r"):
        if not line.startswith("#"):
            split_line = line.split("\t")
            name = "mutation_{}".format(counter)
            var_pos = "{}:{}".format(split_line[0], split_line[1])
            ref = split_line[3]
            alt = split_line[4]

            new_line = "{}\t{}\t{}\t\n".format(name, var_pos, ref+"/"+alt)
            # append each line to a file and use as input to the above
            out_file.write(new_line)
            counter += 1
    out_file.close()

# convert a vcf file to a varinat validation list
vcf2validationlist(vcf, "test_out.tsv", 4)
file_of_interest = open("test_out.tsv")

for line in file_of_interest:
    split_line = line.split("\t")
    name = split_line[0]
    var_pos = split_line[1]

    try:
        mut = split_line[2]
    except IndexError:
        # have a func here that trys and autodetects the mutation
        continue

    ref, alt = tuple(mut.split("/"))
    # determine whether mut is a SNP, deletion or insertion
    if len(ref) == len(alt):
        print("We are looking for a SNP in {}".format(name))
    elif len(ref) > len(alt):
        print("We are looking for a DELETION in {}".format(name))
    elif len(ref) < len(alt):
        print("we are looking for an INSERTION in {}".format(name))



