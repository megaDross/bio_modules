''' Functions used across multiple scripts'''

def correct_hg_version(hg_version):
    ''' Change genome build name'''
    hg_dict = {'grch37': 'hg19', 'grch38': 'hg38'}
    if hg_version.lower().startswith('g'):
        hg_version = hg_dict.get(hg_version.lower())
    return hg_version


