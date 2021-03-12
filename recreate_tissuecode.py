import pandas as pd
import os
import numpy as np

os.chdir("C:/Users/Pavan Nayak/Documents/UCI/Winter 2021 Q11/ECO_EVO_283")
shortRNA = pd.read_table('RNAprefixlist.txt',header=None,names=['FullSampleName'])
shortRNA["TissueCode"] = np.nan

for i in range(0,len(shortRNA)):
    pre = str(shortRNA.iloc[i,0])
    tcode=pre[5]
    shortRNA.iloc[i,1] = tcode

shortRNA.to_csv('shortRNAlist1.txt', sep='\t')
