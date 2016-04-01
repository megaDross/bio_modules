sam_file = "mixup_sample_investigation.txt"
vcf = "var.both.TAAD_UK_Nov2015.filters.vcf"    # doesn't contain anything like 23 etc

# this is designed to work with VCF headers containing sample names
# if a sample name in a given list matches the one in the VCF header then,
# the index of the split header is returned. 

sam_names = ["24BO-C0911","24SY0793",
            "22GK1188","24HS1554",
            "21-GW-0002","23ND1907",
            "23DD1910","29XX1914",
            "23JM1912","23KM1908"]

def get_index_number(f):
    with open(f) as input_file:
        row = [i for i in input_file if "##" not in i] # skip unwanted headers
        head = next(iter(row))    # generator to deal with the header line only. 
        split_head = head.split("\t")
        for sam_name in sam_names:
                if sam_name in split_head:
                    index = split_head.index(sam_name) # get the split_head index of the match
                    return index


index = get_index_number("test.txt")
print index