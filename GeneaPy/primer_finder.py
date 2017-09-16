import pandas as pd
import warnings

def main(db, i, size=None, distance=None, gc=None):
    ''' Find primers for a set of given variants.
    
    Args:
        db: a primer database (tsv format)
        i: variants to find primers for (tsv format)
        size: maximum product size of the resulting amplicon
        distance: minimum distance from F/R primers
        gc: maximum GC% of resulting amplicon
    
    Returns:
        A tsv file containing variants and their matching 
        primers. 
        
    Notes:
        primer database should look like pd_example.tsv
        i should look like:
            variant_name    15:48765543
    '''
    # Only way to get rid of the SettingWithCopyWarning
    warnings.filterwarnings('ignore')

    # change to DataFrame, find matching primers and apply filters
    db_df, var_df = files2df(db, i)
    primers = get_primers(db_df, var_df)
    primers = extra_filters(primers, size, distance, gc)
    
    # list of variants with no primer match
    no_match = var_df[~var_df['Variant'].isin(primers['Variant'])]['Variant'].tolist()
    if no_match:
        err = 'No primer pairs could be found for {}'.format(', '.join(no_match))
        print(err)
    
    return primers


def files2df(db, i):
    ''' Transform primer database and variant files to DataFrames.
    '''
    # prepare primer database
    db = pd.read_csv(db, delimiter='\t')
    db['Chrom'], db['Start'], db['End'] = db.Primer_Range.str.split('[-:]').str
    db = db.drop('Primer_Range', axis=1)

    # prepare input file
    i = pd.read_csv(i, delimiter='\t', header=None)
    i.columns = ['Variant', 'Variant_Position']
    i['Chrom'], i['Pos'] = i['Variant_Position'].str.split(':').str
    i = i.drop('Variant_Position', axis=1)

    return (db, i)


def get_primers(db_df, var_df):
    ''' Returns a DataFrame of variants and their matching primers.

    Args:
        db_df: primer database DataFrame
        var_df: variants of interest DataFrame
    '''
    # merge input and database by Chrom
    m = pd.merge(var_df, db_df)

    # convert some columns to numeric type
    m['GC%'] = m['GC%'].str.replace('%', '')
    m['Product_Size'] = m['Product_Size'].str.replace('bp', '')
    nums = ['Chrom', 'Pos', 'Start', 'End', 'Product_Size', 'GC%',
            'Number_Amplicons']
    m[nums].apply(pd.to_numeric, errors='coerce', axis=1)
    for c in nums:
        m[c] = pd.to_numeric(m[c], errors='coerce')

    # filter for those within primers
    within_primers = ((m['Start'] < m['Pos']) & (m['Pos'] < m['End']))
    output = m[within_primers]

    # variant distance from start and end of primer
    output['Dist_F'] = output['Pos'] - output['Start']
    output['Dist_R'] = output['End'] - output['Pos']
    output = output.drop(['Start', 'End'], axis=1)
    
    return output


def extra_filters(primer, size, distance, gc):
    ''' Filter primers on product size, distance and
        product GC%.

    Args: 
        primers:DataFrame containing variants and matching primers
        size: maximum product size
        distance: minimum distance from F/R primers
        gc: maximum product GC%
    '''
    
    if size:
        primer = primer[primer.Product_Size < 500]
    if distance:
        primer = primer[((primer.Dist_F > distance) &
                         (primer.Dist_R > distance))]
    if gc:
        primer = primer[primer['GC%'] < gc]
        
    return primer


if __name__ == '__main__':
    i = '/home/david/projects/GeneaPy/test/primer_finder_in.txt'
    db = '/home/david/projects/GeneaPy/test/primer_database.txt'
    p = main(db, i, size=500, gc=50, distance=100)
    print(p.shape)
    print(p)
