import os
import tkinter as tk
import tkinter.scrolledtext as scrolledtext
from dataclasses import dataclass
import sys

def resource_path(relative_path):
    try:
        base_path = sys._MEIPASS
    except Exception:
        base_path = os.path.abspath(".")

    return os.path.join(base_path, relative_path)

#Automate the path changing to the directory where the .csv is
#dir_path = os.path.dirname(os.path.realpath(__file__))
#os.chdir(dir_path)

class baseGUI: 
    def __init__(self):
        self.window = tk.Tk()
        self.window.geometry("1100x500")

        #domain variables
        self.keywords = ""  # Initialize keywords variable
        self.subjects = ""

        #label of the textbox
        self.research = tk.Label(self.window, text = "Keywords input \n Separate by a comma to differentiate them \n Currently keywords setting is set to OR", font=("Calibri", 11))        
        self.research.pack(padx=0, pady = 5)        
        
        #textbox from which we want to get the input
        self.textbox_research = tk.Text(self.window, height = 1, font=("Calibri", 12))
        self.textbox_research.pack(padx=0, pady=0)
        
        #logs label
        self.logs_label = tk.Label(self.window, text = "Results", font=("Calibri", 11))        
        self.logs_label.pack(padx=0, pady = 0)  

        #Logs box
        self.logs = scrolledtext.ScrolledText(self.window, height = 15, font=("Calibri", 12), yscrollcommand=True)
        self.logs.pack(padx=0, pady=0)

        #start the research process
        self.button = tk.Button(self.window, text= "Start researching", command= self.confirm)
        self.button.place_configure(relx = 0.85, rely= 0.05)
            
        self.window.mainloop()
        
    def confirm(self):
        pass 

    def keywords_treatment(self):
        pass
    
    def update(self):
        pass


class iGEM__team_search_GUI(baseGUI):
    def confirm(self):
        self.button.config(text='Running...')
        self.logs.insert(tk.END, "The program started, please wait...\n")
        self.button.after(1000, self.keywords_treatment)

    def keywords_treatment(self):
        self.keywords = str(self.textbox_research.get(1.0, "end-1c")).split(',')
        print(self.keywords)

        results = subject_finder(self.keywords)

        self.subjects = ['\t'.join(line.split('\t')[0:2]) for line in results]

        self.logs.insert(tk.END, '\n'.join(set(self.subjects)))
        self.logs.insert(tk.END, "\nResearch complete !\n")

        self.button.config(text="Start research")

@dataclass
class team:
    url : str
    village : str
    abstract : str

def subject_finder(keywords):
    results = []
    data = open('data\\all_team_data.txt', 'r', encoding = 'utf-8')
    
    all_abstract = [team(url = line.split('\t')[0],
                         village = line.split('\t')[1],
                         abstract = line.split('\t')[-1].strip('\n').lower()) for line in data.readlines() if line.split('\t')[-1].strip('\n').lower() != '']
    
    n = 1
    for i in all_abstract:
        n+=1
        for word in keywords :
            if word.lower() in i.abstract.lower():
                results.append(i.url + '\t' + i.village + '\t' + i.abstract + '\n')
    data.close()
    return results         

if __name__ == '__main__':
    iGEM__team_search_GUI()
