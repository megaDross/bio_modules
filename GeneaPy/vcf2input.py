vcf = "/home/david/Downloads/FB_0031_GTGAAA_L001_pe.vcf"

def vcf2input(vcf, filename, starting_num=1):
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
#vcf2validationlist(vcf, "test_out.tsv", 4)

