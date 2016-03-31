from __future__ import division
import re, os, sys

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
    
    if sys.platform == "cygwin" or sys.platform == "linux2":
        file_path = file_path + "/"
        return file_path
        
def gc_content(sequence):
    gc_percent = round(((sequence.count("G")+sequence.count("C")+sequence.count("c")+\
                    sequence.count("g"))/ len(sequence)) * 100,2)
    return gc_percent
    
    

    