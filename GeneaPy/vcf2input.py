vcf = "/home/david/Downloads/FB_0031_GTGAAA_L001_pe.vcf"

def vcf2input(vcf, filename):
    ''' Create a variant validation tsv file that can be used
        as input to get_seq.py. This will create a unique id
        for each sample/variant.
    '''
    counter = 0  
    out_file = open(filename, "w")

    # store the VCF header here
    header = ""
   
    # nested loop for every sample and every line
    for line in open(vcf, "r"):
        if not line.startswith("##"):
            if line.startswith("#"):
                header = line.strip("\n")
            else:
                for num in range(9,len(header.split("\t"))):         
                    split_line = line.split("\t")
                    
                    # get the ID, var position and ref/alt
                    var_pos = "{}:{}".format(split_line[0], split_line[1])
                    variant = "/".join((split_line[3], split_line[4]))
                    ID = split_line[2]

                    # assuming this is a multi sample vcf then the sample name will be in the header
                    SAM = header.split("\t")[num]

                    # get a unique ID for each sample/variant in the VCF
                    name = "{}-{}-{}".format(ID, SAM, counter)
                    
                    # append each line to a file
                    new_line = "{}\t{}\t{}\t\n".format(name, var_pos, variant)
                    out_file.write(new_line)

                    counter += 1
    out_file.close()


# convert a vcf file to a varinat validation list
#vcf2validationlist(vcf, "test_out.tsv", 4)

