# GeneaPy

A set of python3 tools that allows one to manipulate sequences and PCR primers 

## primer_finder
Takes variant position(s) as input and matches it with an appropriate primer in a given file containing primer information

### Example Usage
```bash
# find if a primer pair is available for a single genomic position
python3 primer_finder.py 15:48729400 --primer_database primer_database.txt

# find if primer pairs are available for multiple positions in a given input file and ouput to file
python3 primer_finder.py --input_file in.txt --primer_database primer_database.txt --output_file out.txt
```

## Testing
To perform testing execute the following command within the test dir:
  python3 -m unittest *.py
