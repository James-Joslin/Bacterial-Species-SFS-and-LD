# -*- coding: utf-8 -*-
"""
Created on Tue Oct 20 11:01:00 2020

@author: James
"""

# =============================================================================
#     Current_Strains = [ele for ele in cog if all(ch in ele for ch in Species[b])]
#     Filtered_Strains = []
#     for string in Current_Strains:
#         new_string = string.replace("ATGC002.COG00001"+Species[b], "")
#         Filtered_Strains.append(new_string)
#     Strains = [ x for x in Filtered_Strains if "ATGC" not in x ]
#     Final_Strains = []
#     for string in Strains:
#         new_string = string.replace(".fasta", "")
#         Final_Strains.append(new_string)
# =============================================================================

import gc
import re
import glob
import os
import pandas as pd
import numpy as np
from tqdm import tqdm
from time import sleep
  

Nucleotide_Table = {
    'A':'Adenine','T':'Thymine',
    'G':'Guanine','C':'Cytosine',
    '-' : 'Indel'}
#LoopCount = 1
cog = os.listdir("C:/Users/james/Desktop/Post_Grad_Files/ATGC002/Fasta_Rows2")
Species = []
IDs = []

for g in range(len(cog)):
    cog_temp = cog[g]
    cog_temp = cog_temp.split(".")
    cog_temp = cog_temp[1]
    cog_temp_number = re.findall(r'\d+', cog_temp)
    cog_temp_number = cog_temp_number[0]
    Species_temp = cog_temp.strip("COG"+cog_temp_number)
    IDs.append(cog_temp_number)
    Species.append(Species_temp)
IDs = np.unique(IDs)
IDs = IDs.tolist()
Species = np.unique(Species)
Species = Species.tolist()


for a in range (len(Species)):
    Species_temp = Species[a]
    Species_temp = Species_temp.split("_")
    Species_temp = Species_temp[0] + "_" + Species_temp[1]
    Species[a] = Species_temp
Species = np.unique(Species)
Species = Species.tolist()
Species_Extra = Species[2]

# pneumoniae duplicate below
Species.append(Species_Extra)


CleanedSpecies = []

New_Species = glob.glob("C:/Users/james/Desktop/Post_Grad_Files/ATGC002/Full_Strains/*.txt")
for w in range(len(New_Species)):
    List_NewSpecies = open(New_Species[w])
    List_NewSpecies = List_NewSpecies.read()
    List_NewSpecies = List_NewSpecies.split("\n")
    CleanedSpecies.append(List_NewSpecies)

with tqdm(total=100) as pbar:
    for b in range(1):
        b = 3
        Selected_Species = CleanedSpecies[b]
        for c in range(2249, len(IDs), 1):
            sleep(0.0001)
            pbar.update(100/(len(IDs)-2249))
            gc.collect()
            rSquared_Parameters = {'RSquared' : [],
                 'Distance' : [],
                 'Species' : []}
            rSquared_Results = pd.DataFrame(rSquared_Parameters, columns = ['RSquared',
                 'Distance',
                 'Species'])
            sequencesCombined = []
            Poly_Points_List = []
            for d in range(len(Selected_Species)):
                gene = 'pass'
                filename = "ATGC002.COG"+IDs[c]+Species[b]+Selected_Species[d]+".fasta"
                #print(filename)
                seqdirectory = "C:/Users/james/Desktop/Post_Grad_Files/ATGC002/Fasta_Rows2/"+str(filename) 
                infile = open(seqdirectory,'r')
                sequences = infile.read()
                sequences_split = sequences.split()
                sequencesCombined.append(sequences_split[1])   
            NumberHyphens = sequences.count('-')
            #print(NumberHyphens/len(sequencesCombined[0]))
            if (NumberHyphens/len(sequencesCombined[0])) >= 0.1:
                gene = 'fail'
            else:
                for i in range(0, len(sequencesCombined[0]), 1):
                    Nucleotides_Dict = {}
                    for j in range(len(sequencesCombined)): 
                        if gene == 'pass':
                            Nucleotide = sequencesCombined[j][i:i+1] 
                            if Nucleotide not in Nucleotides_Dict:
                                Nucleotides_Dict[Nucleotide] = 1 
                            else: 
                                Nucleotides_Dict[Nucleotide] += 1
                        else:
                            pass
                    if len(Nucleotides_Dict) == 2: 
                             Nucleotides_List = []
                             for x in list(Nucleotides_Dict.keys()):
                                 Nucleotides_List.append(Nucleotide_Table[x])
                             if Nucleotides_List[0] == Nucleotides_List[1]: 
                                 pass
                             else:
                                 if Nucleotides_List[0] == "Indel" or Nucleotides_List[1] == "Indel":
                                     pass
                                 else:
                                     Poly_Point = i
                                     Poly_Points_List.append(Poly_Point)
                    else:
                        pass
                for h in range(0, (len(Poly_Points_List)-1), 1):
                    for q in range(h + 1, len(Poly_Points_List)):
                        dinucleotide_Dict = {}
                        rSquared_Stats = [0,0]
                        for e in range(len(sequencesCombined)):
                            dinucleotide = sequencesCombined[e][Poly_Points_List[h]:Poly_Points_List[h]+1]+sequencesCombined[e][Poly_Points_List[q]:Poly_Points_List[q]+1]
                            if dinucleotide not in dinucleotide_Dict:
                                dinucleotide_Dict[dinucleotide] = 1
                            else:
                                dinucleotide_Dict[dinucleotide] += 1
                        dinucleotide_list = list(dinucleotide_Dict.keys())
                        total = len(sequencesCombined)
                        
                        #First Site
                        firstNucleotide = dinucleotide_list[0][0]
                        firstCount = 0
                        
                        for dinucleotide in dinucleotide_Dict:
                            if dinucleotide[0] == firstNucleotide:
                                firstCount += dinucleotide_Dict[dinucleotide]
                        firstFrequency = firstCount/total
                        
                        #Second Site
                        secondNucleotide = dinucleotide_list[0][1]
                        secondCount = 0
                        
                        for dinucleotide in dinucleotide_Dict:
                            if dinucleotide[1] == secondNucleotide:
                                secondCount += dinucleotide_Dict[dinucleotide]
                        secondFrequency = secondCount/total
                           
                        haplotypeFrequency = dinucleotide_Dict[dinucleotide_list[0]]/total
                        
                        D = (haplotypeFrequency - firstFrequency * secondFrequency)
                        
                        rSquared = D ** 2 / (firstFrequency*(1-firstFrequency)*secondFrequency*(1-secondFrequency))
                        
                        dist =  Poly_Points_List[q] - Poly_Points_List[h]
                        
                        RSquared_Values = {'RSquared' : [rSquared],
                                           'Distance' : [dist],
                                           'Species' : [Species[b]]}
                        
                        RSquared_Row = pd.DataFrame(RSquared_Values, columns = ['RSquared',
                                                                                'Distance',
                                                                                'Species'])
                                                    
                        rSquared_Results = pd.concat([rSquared_Results, RSquared_Row], axis=0, ignore_index = True)
    
                rSquared_Results.to_csv("C:/Users/james/Desktop/Post_Grad_Files/ATGC002/9_Gene_Disequilibrium_Cleaned_2/" + "COG"+str(IDs[c])+str(Species[b]) + "RSquared_Results_Cleaned" + ".csv", 
                                 index = False) 