# -*- coding: utf-8 -*-
"""
Created on Tue Feb  2 19:42:22 2021

@author: james
"""
from tqdm import tqdm
from time import sleep
import re
import glob
import os
import pandas as pd
import numpy as np

AminoAcid_Table = {
    'TTT':'Phe','TCT':'Ser','TAT':'Tyr','TGT':'Sys', 
    'TTC':'Phe','TCC':'Ser','TAC':'Tyr','TGC':'Sys',
    'TTA':'Leu','TCA':'Ser','TAA':'Stop','TGA':'Stop',
    'TTG':'Leu','TCG':'Ser','TAG':'Stop','TGG':'Trp',
    'CTT':'Leu','CCT':'Pro','CAT':'His','CGT':'Arg',
    'CTC':'Leu','CCC':'Pro','CAC':'His','CGC':'Arg',
    'CTA':'Leu','CCA':'Pro','CAA':'Gln','CGA':'Arg',
    'CTG':'Leu','CCG':'Pro','CAG':'Gln','CGG':'Arg',
    'ATT':'Ile','ACT':'Thr','AAT':'Asn','AGT':'Ser',
    'ATC':'Ile','ACC':'Thr','AAC':'Asn','AGC':'Ser',
    'ATA':'Ile','ACA':'Thr','AAA':'Lys','AGA':'Arg',
    'ATG':'Met','ACG':'Thr','AAG':'Lys','AGG':'Arg',
    'GTT':'Val','GCT':'Ala','GAT':'Asp','GGT':'Gly',
    'GTC':'Val','GCC':'Ala','GAC':'Asp','GGC':'Gly',
    'GTA':'Val','GCA':'Ala','GAA':'Glu','GGA':'Gly',
    'GTG':'Val','GCG':'Ala','GAG':'Glu','GGG':'Gly'}

Folder = glob.glob("C:/Users/james/Desktop/Post_Grad_Files/ATGC002/Full_Strains/" + "/*txt")

substringa = "_large"
substringb = "_small"


with tqdm(total=100) as pbar:
    for w in range(len(Folder)):
        open_file = open(Folder[w], 'r')
        read_file = open_file.read()
        Strains = read_file.split("\n")
        All_Clades = []
        SFSn_All = []
        SFSs_All = []
        df_name = Folder[w]
        df_name = df_name.split("\\")
        df_name = df_name[1].split(".")
        df_name = df_name[0]
        
        Gene_Str = []
        sequencesCombined = []
        Gaps_All = ["Gap"]
        Complex_All = ["ComplexCodons"]
         
        fullstring = df_name
        if substringa in fullstring:
            species = fullstring.replace("_large", "", 1)
        elif substringb in fullstring:
            species = fullstring.replace("_small", "", 1)
        else:
            species = fullstring
        
        cog = os.listdir("C:/Users/james/Desktop/Post_Grad_Files/ATGC002/Fasta_Rows2")
        for g in range(len(cog)):
            cog_temp = cog[g]
            cog_temp = cog_temp.split(".")
            cog_temp = cog_temp[1]
            cog_temp = re.findall(r'\d+', cog_temp)
            cog_temp = cog_temp[0]
            cog[g] = cog_temp
        cog = np.unique(cog)
        cog = cog.tolist()
    
        for l in range(len(cog)):
            sleep(0.0001)
            pbar.update(100/(len(cog)*len(Folder)))
            SFSn = [0]*int(len(Strains)/2)  #Your annotation - you have to be very careful in constructing lists in this way; sometimes 
            SFSs = [0]*int(len(Strains)/2)  #inadvertently construct references to lists and niot the lists themselves
            sequencesCombined = []
            Complex_Codons = 0
            NumGaps = 0
            TotalNumGaps = 0
            gene = 'pass'
    
            for z in range(len(Strains)):
                filename = "ATGC002.COG"+cog[l]+str(species)+str(Strains[z])+'.fasta' 
                seqdirectory = "C:/Users/james/Desktop/Post_Grad_Files/ATGC002/Fasta_Rows2/"+str(filename) 
                isdir = os.path.isfile(seqdirectory)
                #print(isdir)
                if isdir == True:
                    infile = open(seqdirectory,'r')
                    sequences = infile.read()
                    sequences_split = sequences.split()
                    sequencesCombined.append(sequences_split[1])   
                                    #I now realise that instead of doing this manually I could've multiplied the length of the zeros, by half the length of the strain list
                                    #And if it was a decimal, then I know that this decimal would be ".5" then perhaps I could've written a couple of lines to round this down to the closest
                                    #whole number, not sure how to do this though, but I reckon it would be doable.
            if isdir == True:  
                NumberHyphens = sequences.count('-')
                #print(NumberHyphens/len(sequencesCombined[0]))
                if (NumberHyphens/len(sequencesCombined[0])) >= 0.1: #Check if number of Gaps makes up less than 10%
                    gene = 'fail' #If the gene fails the loop continues but nothing gets written and so the gene is ignored in the final output
                else:
                    for i in range(0, len(sequencesCombined[0]), 3): #repeat process for the length of the sequence moving up in intervals of 3
                        codons = {} #Build an empty dictionary ready for populating
                        NumGaps = 0
                        for j in range(len(sequencesCombined)): #Repeat loop for the number of sequences in my list, in this case 10
                            if gene == 'pass':
                                codon = sequencesCombined[j][i:i+3] #slice every 3 elements of the sequence to form the codon
                                if (codon.count('-') == 1) or (codon.count('-') == 2): #check that slice doesn't have any gaps that are't a factor of 3
                                    gene = 'fail' #if the current slice has gaps that aren't a factor of 3 then fail the entire gene
                                else:
                                    if '---' in codon: #for counting the number of gaps in the passed sequences
                                        NumGaps += 1
                                    else:   #what codons do we have?             
                                        if codon not in codons:
                                            codons[codon] = 1 #if this codon isn't present in the dictionary, add it to the dictionary, if it is present
                                        else: #add one onto the total
                                            codons[codon] += 1
                            else:
                                pass
                        if (NumGaps == 0): # If there weren't gaps present in the slice check the SFS of that slice, if there were gaps, then don't write the SFS for that slice and tally the number of gaps instead
                            if gene == 'pass':
                              if len(codons) == 2: #if there are two codons within the dictionary 
                                  amino_acid = [] #form an empty table called amino_acid
                                  for x in list(codons.keys()): # add the amino acids from the amino acid dictionary, that the codons of the codon dictionary 
                                      amino_acid.append(AminoAcid_Table[x]) #code for to the amino acid_table
                                  if amino_acid[0] == amino_acid[1]: #check if these two amino acids are the same
                                      sub = min(codons, key=codons.get) # if they are the same do this 
                                      #print('Syn - codon: ' +  str(sub)) # sub is simply a string which tells me the least common of the amino acids               
                                      sub = codons[sub]
                                      SFSs[sub-1] += 1
                                  else:
                                      sub = min(codons, key=codons.get) # if they are not the same do this:
                                      #print('Nsyn - codon: ' +  str(sub))
                                      sub = codons[sub]
                                      SFSn[sub-1] += 1
                              elif len(codons) > 2:
                                  Complex_Codons += 1
                              else:
                                  pass
                            else:
                                pass
                        else:
                            TotalNumGaps += 1 #the total number from the cumulation of gaps gets tallied here instead
                if gene == 'pass': #if the gene passed all of the tests repeated for each slice, the cumulative information over all the slices is appended into a single line within a list
                    Gene_Str.append("COG" + str(cog[l]))
                    SFSn_str = str(SFSn)
                    SFSn_str = SFSn_str.strip('[]')
                    SFSn_All.append(SFSn_str)
                    SFSs_str = str(SFSs)
                    SFSs_str = SFSs_str.strip('[]')        
                    SFSs_All.append(SFSs_str)
                    Gaps_Str = str(TotalNumGaps)
                    Gaps_Str = Gaps_Str.strip('[]')
                    Gaps_All.append(Gaps_Str)
                    Complex_Codons_Str = str(Complex_Codons)
                    Complex_Codons_Str = Complex_Codons_Str.strip('[]')
                    Complex_All.append(Complex_Codons_Str)
                else:
                    pass
            else:
                pass
        SeriesID = pd.Series(Gene_Str)     
        SeriesPN = pd.Series(SFSn_All)     
        SeriesPS = pd.Series(SFSs_All)     
        Final_df = pd.concat([SeriesID, SeriesPN, SeriesPS], axis=1)    
        Final_df.to_csv("C:/Users/james/Desktop/Post_Grad_Files/ATGC002/5_SFS/" + df_name + "_SFS.csv", 
                                 index = False, header = False, sep = ",")

