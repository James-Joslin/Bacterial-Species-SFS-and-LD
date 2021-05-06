# -*- coding: utf-8 -*-
"""
Created on Mon Nov 16 19:05:21 2020

@author: James

Mean Distance.             No. Of values             Mean rSquared
1                          500                       0.6
(80*2+40*3)/120            (80+40) = 120             (80*0.3+40*0.4)/120
(30*4+50*5+30*6)/110       (40+30+50+30) = 110       (30*0.2+50*0.2+30*0.1)/110

"""


import pandas as pd
import numpy as np
import glob



RollingCount = 0
RollingRSquared = 0
RollingDistances = 0

#df = pd.read_csv("F:/PostGrad_Research/Polymorphic_Data/Clade_SFSs/OxytocaAll.csv")
Threshold = 50000

path = "C:/Users/james/Desktop/Post_Grad_Files/ATGC002/9_Gene_Disequilibrium_Cleaned" # use your path
all_files = glob.glob(path + "/*.csv")

for a in range(len(all_files)):
    df = pd.read_csv(all_files[a])
    rSquaredSums = pd.Series([])
    Counts = pd.Series([])
    Distances = pd.Series([])
    SeriesSpecies = []
    SeriesSFS = []
    RollingCount = 0
    RollingDistances = 0
    RollingRSquared = 0    
    
    for i in range(len(df)):
        count = df["Counts"].iloc[i]
        distance = df["Distance"].iloc[i]
        rSquared = df["RSquared"].iloc[i]
        species = df["Species"].iloc[i]
        SFS = df["SFS"].iloc[i]
        
        RollingCount += count
        
        RollingDistanceTemp = distance*count
        RollingDistances += RollingDistanceTemp  
        
        RollingRSquaredTemp = rSquared*count
        RollingRSquared += RollingRSquaredTemp 
        
        if RollingCount >= Threshold:
            Counts = Counts.append(pd.Series([RollingCount]), ignore_index = True)
            
            DistancesTemp = RollingDistances/RollingCount
            Distances = Distances.append(pd.Series([DistancesTemp]), ignore_index = True)
            
            rSquaredTemp = RollingRSquared/RollingCount
            rSquaredSums = rSquaredSums.append(pd.Series([rSquaredTemp]), ignore_index = True)
            
            SeriesSpecies.append(species)
            
            SeriesSFS.append(SFS)
            
            RollingCount = 0
            RollingDistances = 0
            RollingRSquared = 0
        else:
            pass
        
        # if i == (len(df)-1):
        #     Counts = Counts.append(pd.Series([RollingCount]), ignore_index = True)
            
        #     DistancesTemp = RollingDistances/RollingCount
        #     Distances = Distances.append(pd.Series([DistancesTemp]), ignore_index = True)
            
        #     rSquaredTemp = RollingRSquared/RollingCount
        #     rSquaredSums = rSquaredSums.append(pd.Series([rSquaredTemp]), ignore_index = True)
            
        #     SeriesSpecies.append(species)
            
        #     SeriesSFS.append(SFS)
            
        #     RollingCount = 0
        #     RollingDistances = 0
        #     RollingRSquared = 0
        # else:
        #     pass        
        
    SeriesSpecies = pd.Series(SeriesSpecies) 
    SeriesSFS = pd.Series(SeriesSFS) 
    
    Final_df = pd.concat([Counts, Distances, rSquaredSums, SeriesSpecies, SeriesSFS], axis=1)    
    Final_df.columns=["Counts", "MeanDistances", "MeanRSquared", "Species", "SFS"]
    Final_df.to_csv("C:/Users/james/Desktop/Post_Grad_Files/ATGC002/10_Gene_Disequilibrium_Smoothed/" + species + "_" + SFS + ".csv", 
                             index = False) 