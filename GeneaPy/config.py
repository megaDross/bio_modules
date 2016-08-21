import configparser
import os, sys

home = os.path.expanduser("~") 

default_config = { 

    "ttuner_path": home+"/bin/tracetuner_3.0.6beta/rel/Linux_64/ttuner",
    "default": "hg19",
    "hg19_path": home+"/.config/genome/hg19.fa",
    "hg38_path": home+"/.config/genome/hg38.fa"
}


def config():
    ''' Ask the user if they wish to use the default config
    '''
    answer = input("Do you want to use the default configuration file (y or n)?\n")

    if answer in ("y", "Y", "yes", "Yes", "YES"):
        write_config(default_config.get("ttuner_path"), default_config.get("hg19_path"),
                     default_config.get("hg38_path"), default_config.get("default"))
    else:
        write_config()



def write_config(ttuner_path=None, hg19_path=None, hg38_path=None, default=None):
    ''' write configuration file for get_seq.py
    '''
    # ask user for variables if they are empty
    if None in (ttuner_path, hg19_path, hg38_path, default):

        ttuner_path = input("Please specify the path to ttuner:\n")
        hg19_path = input("Please specify the path to the human genome GRcH 37 (hg19) fasta file:\n")
        hg38_path = input("Please specify the path to the human genome GRcH 38 (hg38) fasta file:\n")
        default = input("Which genome is preffered as the default (hg19 or hg38)?\n")
    
    # write the file with the given variables
    config = configparser.ConfigParser()
    config['TTUNER'] = {'Path': ttuner_path}
    config['GENOME'] = {'DefaultGenome': default, 'PathHG19': hg19_path, 
                        'PathHG38': hg38_path}
    with open(home + '/.config/GeneaPy.config', 'w') as configfile:
        config.write(configfile)




def read_config(default=None):
    ''' get the items from the config file and return them
    '''
    # read the config file
    config = configparser.ConfigParser()
    config.read(home + "/.config/GeneaPy.config")
    
    # get the variables from the config file and return them
    ttuner = config["TTUNER"]["Path"]
    if not default:
        default = config["GENOME"]["DefaultGenome"]
    local_genome = config["GENOME"]["Path" + default.upper()]

    return (ttuner, default, local_genome)




