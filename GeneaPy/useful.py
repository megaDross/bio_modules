from __future__ import division
import re, os, sys

complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A',
              'Y': 'R', 'R': 'Y', 'W': 'W', 'S': 'S', 'K': 'M', 'M': 'K',
              'D': 'H', 'V': 'B', 'H': 'D', 'B': 'V',
              'N': 'N'}


class conditional_decorator(object):
    ''' Allows own to only use a decorator if a certain condition has been met.
        
        Example usage:
        @conditional_decorator(transcribe, over_10bp)
        def get_seq(.....)
    '''
    def __init__(self, dec, condition):
        self.decorator = dec
        self.condition = condition

    def __call__(self, func):
        if not self.condition:
            # return the function unchanged; not decorated
            return func
        return self.decorator(func)


def sorted_nicely(l):
    """ Sorts the given iterable in the way that is expected (alphanumerically).
 
    Required arguments:
    l -- The iterable to be sorted.
 
    """
    convert = lambda text: int(text) if text.isdigit() else text
    alphanum_key = lambda key: [convert(c) for c in re.split('([0-9]+)', key)]
    return sorted(l, key = alphanum_key)


def cwd_file_path(x):
    ''' Return the absolute path for the file which is currently 
        being executed
        
        Required arguments:
        x -- always = __file__
    '''
    file_path = os.path.dirname(os.path.abspath(x))
    
    if sys.platform == "win32":
        file_path = file_path + "\\"
        return file_path
    
    if sys.platform == "cygwin" or sys.platform == "linux2" or sys.platform == "linux":
        file_path = file_path + "/"
        return file_path

def gc_content(sequence):
    ''' calculates the GC content of a given DNA sequence
    '''
    gc_percent = round(((sequence.count("G")+sequence.count("C")+sequence.count("c")+\
                    sequence.count("g"))/ len(sequence)) * 100,2)
    return gc_percent

def chunks(l, n):
    ''' Yield successive n-sized list from list(l)

        To produce a nested list: 
            print(list(chunks(l,n)))
    '''
    for i in range(0, len(l), n):
        yield l[i:i+n]

def filter_list_of_lists(list_of_lists):
    ''' Filter out empty strings from all lists in a list of lists
    '''
    filtered_lists = []
    for row in list_of_lists:
        filtered = []
        for element in row:
            if element:
                filtered.append(element)
        filtered_lists.append(filtered)

    return filtered_lists


def reverse_complement(seq):
    ''' reverse complement a given DNA sequence
    '''
    reverse_complement = "".join(complement.get(base, base) for base in reversed(seq))
    return(reverse_complement)



