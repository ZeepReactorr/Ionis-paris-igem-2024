import os 
from dataclasses import dataclass
import sys 

@dataclass
class team:
    url : str
    village : str
    abstract : str

dir_path = os.path.dirname(os.path.realpath(__file__))
os.chdir(dir_path)

def subject_finder(keywords):
    main_dir = os.getcwd()
    filename = 'all_team_data.txt'
    results = []
    data = open(os.path.join(main_dir, filename), 'r', encoding = 'utf-8')
    
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
    
#print(subject_finder(['plant virus']))