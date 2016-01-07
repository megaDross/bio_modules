import os
import re

def sorted_nicely( l ):
    """ Sorts the given iterable in the way that is expected (alphanumerically).
 
    Required arguments:
    l -- The iterable to be sorted.
 
    """
    convert = lambda text: int(text) if text.isdigit() else text
    alphanum_key = lambda key: [convert(c) for c in re.split('([0-9]+)', key)]
    return sorted(l, key = alphanum_key)


def rename_episodes():
    # Ask user to type a PATH for the directory in which this program will be executed.
    # Defaults to working directory
    Dir = raw_input("Type the path in which you wish this to take place in (defaults to working dir): " )
    if Dir is "":
        Dir = os.getcwd()
    
    
    # Ask for the SHOW TITLE input. Defaults to the name of the directory which is above
    # the directory stored in Dir. Changes all white spaces to _.
    Title = raw_input("Type the TV shows title (defaults to directory above the path given previously):")
    if Title is "":
        Dirs = Dir.split("/")
        Title = Dirs[len(Dirs)-2]
    
    if " " in Title:
        Title = re.sub(" ","_",Title)
    
    
    # Season and Episode input. While loop if a letter is found, only accepts until an
    # integer is given by the user.
    season = str(raw_input("Type the season number: " ))
    while season.isdigit() == False :
        season =raw_input("No Letters!\ngive me a number: " )
        
    first_episode = eval(raw_input("Type the first episode number in the series: " ))
        
    last_episode = (len(Dir)+first_episode) - 2
    
    
    
    # give the file extension
    extension = raw_input("Give the desired file extension: " )
    
    # defines the new episode range
    new = range(first_episode,last_episode-1)
    
    
    # Warn the user what each file will be renamed to
    for (current),(want) in zip(sorted_nicely(os.listdir(Dir)),new):
        new_filename = str(Title+"_"+str(season)+"x"+str(want)+extension)
        print "This file: "+current+"\n"+"Will be renamed to: "+new_filename+"\n"

    # Ask user if they wish to proceed
    while raw_input("Proceed with renaming(enter Y for yes)?") != "Y":
        raw_input("Proceed with renaming(enter Y for yes or ctrl+C to terminate)?")
    
    # Rename files
    for (current),(want) in zip(sorted_nicely(os.listdir(Dir)),new):
        new_filename = str(Title+"_"+str(season)+"x"+str(want)+extension)
        os.rename(current,new_filename)
        
    