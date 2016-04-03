#!/home/david/bin/Canopy/User/bin/python
from __future__ import division
import requests, bs4
from useful_tools import useful
from Tkinter import *
import tkMessageBox

# create a blank window and give it a title
window = Tk()
window.title("Unkown Primer")

# assign string variables and set defaults
f = StringVar()
r = StringVar()
option = StringVar()

f.set("TAACAGATTGATGATGCATG")
r.set("CCCATGAGTGGCTCCTAAA")
option.set("hg38")

# create labels and entry boxes for forward primer input
textFprimer = Label(window, text="Enter forward primer sequence:",padx=30,pady=30)  # pad defines empty space around label
textFprimer.grid(row=1,column=1)    # position within empty window
boxFprimer = Entry(window,bd=5, textvariable=f)   # bd is border depth, assign variable name
boxFprimer.grid(row=1,column=2)

# create labels and entry boxes for reverse primer input
textRprimer = Label(window,text="Enter reverse primer sequence:",padx=30,pady=30)
textRprimer.grid(row=2,column=1)
boxRprimer = Entry(window,bd=5,textvariable=r)
boxRprimer.grid(row=2,column=2)

# create a label and drop down menu for picking human genome version
textOptions = Label(window,text="Pick a human genome version:",padx=50)
textOptions.grid(row=3,column=1)
option_choice = OptionMenu(window,option,"hg16","hg17","hg18","hg19","hg38")
option_choice.configure(width=10)
option_choice.grid(row=3,column=2)

# create a button that intiates the below function
button = Button(window,text="Go")
button.grid(row=4,column=1)

def get_unknown_primer():
    ''' Generate an amplicon sequence from inputted primer sequences, which
        is further manipulated t derive inofrmation from the sequence.
    '''
    # store variables from entry box and option box
    f_primer = f.get()
    r_primer = r.get()
    hg_version = option.get()
    

    # generate amplicon sequence using isPCR tool and the above variables
    req = requests.get("https://genome.ucsc.edu/cgi-bin/hgPcr?hgsid=483751629_"+\
                        "vuLjoO4UVF9h4vF4TEp9U8OQiFd7&org=Human&db="+hg_version+\
                        "&wp_target=genome&wp_f="+f_primer+"&wp_r="+r_primer+"&Submit=submit"+\
                        "&wp_size=4000&wp_perfect=15&wp_good=15&boolshad.wp_flipReverse=0")
    req.raise_for_status()                                # get the request error code if failed
    entire_url = bs4.BeautifulSoup(req.text,"html.parser")
    pre_elements = entire_url.select('pre')               # get all <pre> elements on webpage
    
    # if nothing between pre-elements, raise error
    isPCR = pre_elements[0].getText()                     # get text for first <pre> element
    amplicon = "\n".join(isPCR.split("\n")[1:])           # split newlines, take 2nd to last and join back together
    
    # amplicon sequence returned in a new dialog box 
    tkMessageBox.showinfo("Amplicon for "+hg_version,amplicon.replace("\n",""))


# configure the command which intiates upon pressing the button
button.configure(command=get_unknown_primer)

# not sure.... kind of intiates the whole thing
mainloop()      

