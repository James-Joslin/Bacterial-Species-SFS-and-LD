# -*- coding: utf-8 -*-
"""
Created on Sun Nov  1 10:47:54 2020

@author: James
"""

import pandas as pd
import glob
from tqdm import tqdm
from time import sleep

path = "C:/Users/james/Desktop/Post_Grad_Files/ATGC002/9_5_Adjusted_LD" # use your path
all_files = glob.glob(path + "/*.csv")

rSquared_Parameters = {'RSquared' : [],
             'Distance' : [],
             'Species' : []}

rsquared_Means = pd.DataFrame(rSquared_Parameters, columns = ['RSquared',
             'Distance',
             'Species'])
with tqdm(total=100) as pbar:
    for filename in all_files:
        sleep(0.0001)
        pbar.update(100/(len(all_files)))
        df = pd.read_csv(filename, index_col=None, header=0)
        try:
            df['Distance'] = df['Distance'].astype('category')
            Counts = df.groupby('Distance').count()
            Counts = Counts.rename(columns={'RSquared': 'Count'})
            Counts = Counts.reset_index(drop=True)
            
            df_means = df.groupby(['Species','Distance'], as_index=False)['RSquared'].mean()
            df_means['Counts'] = Counts['Count']
            rsquared_Means = pd.concat([rsquared_Means, df_means], axis=0, ignore_index = True)
        except:
            pass
    
rsquared_MeansFinal = rsquared_Means.groupby(['Species','Distance'], as_index=False)['RSquared'].mean()
rsquared_FinalCounts = rsquared_Means.groupby(['Species','Distance'])['Counts'].sum()
rsquared_FinalCounts = rsquared_FinalCounts.reset_index(drop=True)
rsquared_MeansFinal['Counts'] = rsquared_FinalCounts

rsquared_MeansFinal.to_csv("C:/Users/james/Desktop/Post_Grad_Files/ATGC002/9_5_Adjusted_LD/ParsedResults/MeanRSquared.csv", 
                             index = False)