import pandas as pd
import argparse
import warnings
import logging
import GeneaPy.modules.custom_exceptions as ex
from GeneaPy.modules.common import correct_hg_version

logging.basicConfig(filename='primer_finder.error.log',
                    format="%(asctime)s:%(levelname)s:%(message)s")

def primer_finder(db, variant=None, input_file=None, size=None, distance=None, gc=None, 
                  hg=None, gene=None, exon=None, intron=None, output=None):
    ''' Find primers which fulfill the parsed conditions.
    
    Args:
        db: a primer database (tsv format)
        input_file: variants to find primers for (tsv format)
        size: maximum product size of the resulting amplicon
        distance: minimum distance from F/R primers
        gc: maximum GC% of resulting amplicon
        hg: human genome filter
        gene: gene to filter for
        exon: exon number to filter for
        intron: intron number to filter for
        ouput: output file name
    
    Returns:
        a DataFrame containing primer pairs which passed
        the given filters
        
    Notes:
        primer database should look like pd_example.tsv
        input_file should look like:
            variant_name    15:48765543
    '''
    # Only way to get rid of the SettingWithCopyWarning
    warnings.filterwarnings('ignore')
    hg = correct_hg_version(hg) if hg else None
    db_df = database2df(db)
    primers = filter_for_variants(db_df, input_file, variant)
    primers = extra_filters(primers, size, distance, gc,
                            gene, exon, intron, hg)
    if output:
        primers.to_csv(output, sep='\t', index=False)
    return primers

def database2df(db):
    ''' Transform primer database to a DataFrame.'''
    db = pd.read_csv(db, delimiter='\t')
    db['Chrom'], db['Start'], db['End'] = db.Primer_Range.str.split('[-:]').str
    db['Chrom'] = db.Chrom.str.replace('chr', '')
    db = db.drop('Primer_Range', axis=1)
    db['Product_Size'] = db['Product_Size'].str.replace('bp', '')
    nums = ['Chrom', 'Start', 'End', 'Product_Size', 'GC%']
    db = convert2numeric(db, nums)
    return db

def filter_for_variants(db, input_file, variant):
    ''' Filter primer database for given variants'''
    if input_file:
        var_df = input2df(input_file)
        primers = get_variant_file_primers(db, var_df)
        report_unmatched_variants(var_df, primers)
    elif variant:
        primers = get_variant_primers(db, variant)
    else:
        primers = db
    return primers

def input2df(input_file):
    ''' Transform input file to a DataFrame.'''
    df = pd.read_csv(input_file, delimiter='\t', header=None)
    df.columns = ['Variant', 'Variant_Position']
    df['Chrom'], df['Pos'] = df['Variant_Position'].str.split(':').str
    df = df.drop('Variant_Position', axis=1)
    
    df = convert2numeric(df, ['Chrom', 'Pos'])
    return df

def convert2numeric(df, cols):
    ''' Convert the given DataFrame columns to numeric type'''
    df[cols].apply(pd.to_numeric, errors='coerce', axis=1)
    for c in cols:
        df[c] = pd.to_numeric(df[c], errors='coerce')
    return df

def get_variant_file_primers(db, var_file):
    ''' Returns a DataFrame of variants and their matching primers.'''
    m = pd.merge(var_file, db)
    within_primers = ((m['Start'] < m['Pos']) & (m['Pos'] < m['End']))
    output = m[within_primers]
    output['Dist_F'] = output['Pos'] - output['Start']
    output['Dist_R'] = output['End'] - output['Pos']
    output = output.drop(['Start', 'End'], axis=1)
    return output

def report_unmatched_variants(before, after):
    ''' Report variants which have no primer match'''
    try:
        no_match = before[~before['Variant'].isin(after['Variant'])]['Variant'].tolist()
        if no_match:
            raise ex.UnmatchedVariants(no_match)
    except ex.UnmatchedVariants as e:
        for unmatched in e.unmatched:
            logging.info('Cannot find primer for {}'.format(unmatched))

def get_variant_primers(db, var):
    ''' Returns a DataFrame of primers matching a given variant '''
    chrom, pos = var.replace('chr', '').split(":")
    within_primers = ((db['Start'] < int(pos)) & (int(pos) < db['End']))
    output = db[within_primers]
    output['Dist_F'] = int(pos) - output['Start']
    output['Dist_R'] = output['End'] - int(pos)
    return output

def extra_filters(primer, size, distance, gc, gene, exon, intron, hg):
    ''' Filter primers on given argument values.'''
    try:    
        if size:
            primer = primer[primer.Product_Size < size]
        if distance:
            primer = primer[((primer.Dist_F > distance) &
                             (primer.Dist_R > distance))]
        if gc:
            primer = primer[primer['GC%'] < gc]
        if gene:
            primer = primer[primer['Gene'] == gene]
        if exon:
            exon_prep = primer[primer['Exon'] != '-']
            primer['exon_no'], exon_total = exon_prep.Intron.str.split("/").str
            primer = primer[(primer.exon_no == str(exon))]
            primer = primer.drop('exon_no', axis=1)
        if intron:
            intron_prep = primer[primer['Intron'] != '-']
            primer['intron_no'], intron_total = intron_prep.Intron.str.split("/").str
            primer = primer[(primer.intron_no == str(intron))]
            primer = primer.drop('intron_no', axis=1)
        if hg:
            primer = primer[primer['Genome'] == hg]
        if primer.empty:
            raise ex.EmptyDataFrame()
        return primer

    except (ValueError, ex.EmptyDataFrame):
        logging.error('No primers found for the specific filter arguments parsed')
    
def get_parser():
    parser = argparse.ArgumentParser(description='Find primers for a set of given variants.')
    parser.add_argument('-d', '--database', type=str, required=True, help='primer database')
    parser.add_argument('-i', '--input', type=str, help='input file with variant positions')
    parser.add_argument('-v', '--variant', type=str, help='only amplicons that cover the given variant position')
    parser.add_argument('-s', '--size', type=int, help='maximum size of the resulting amplicon')
    parser.add_argument('-k', '--distance', type=int, help='minimum distance from forward/reverse primer')
    parser.add_argument('-c', '--gc', type=int, help='maximum GC%% of the resulting amplicon')
    parser.add_argument('-hg', '--genome', type=str, help='only amplicons produced in specific genome version')
    parser.add_argument('-g', '--gene', type=str, help='only amplicons lying within the given gene')
    parser.add_argument('-e', '--exon', type=str, help='only amplicons produced in a given exon numer')
    parser.add_argument('-n', '--intron', type=int, help='only amplicons lying within the given intron number')
    parser.add_argument('-o', '--output', type=str, help='output the results of processing --input', default='output_primer_finder.txt')
    return parser

def cli():
    parser = get_parser()
    args = vars(parser.parse_args())
    primers = primer_finder(db=args['database'], variant=args['variant'], 
                            input_file=args['input'], size=args['size'], 
                            distance=args['distance'], gc=args['gc'], 
                            hg=args['genome'], gene=args['gene'], 
                            exon=args['exon'], intron=args['intron'], 
                            output=args['output'])
    print(primers)
    

if __name__ == '__main__':
    cli()
