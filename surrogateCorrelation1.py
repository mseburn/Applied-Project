##########################################################
## Consistilator:  consistilatorV2.py                   ##
##  ______     ______     __  __                        ##
## /\  __ \   /\  ___\   /\ \/\ \                       ##
## \ \  __ \  \ \___  \  \ \ \_\ \                      ##
##  \ \_\ \_\  \/\_____\  \ \_____\                     ##
##   \/_/\/_/   \/_____/   \/_____/                     ##
## @Developed by: Plaisier Lab                          ##
##   (https://plaisierlab.engineering.asu.edu/)         ##
##   Arizona State University                           ##
##   242 ISTB1, 550 E Orange St                         ##
##   Tempe, AZ  85281                                   ##
## @Author:  Chris Plaisier                             ##
## @License:  GNU GPLv3                                 ##
##                                                      ##
## If this program is used in your analysis please      ##
## mention who built it. Thanks. :-)                    ##
##########################################################

from subprocess import *
from multiprocessing import Pool, cpu_count, Manager
from scipy.stats import pearsonr
import numpy as np
import pandas as pd
import json

# Structure of Studer datasets
with open('series.json','r') as inFile:
    subsets = json.load(inFile)

# Read in expression data
expression = {}
for gse in subsets:
    expression[gse] = pd.read_csv('../'+gse+'/'+gse+'_genesExpMatrix_all.csv', header=0, index_col=0).transpose()

# Load up surrogate variables
metaData = pd.read_csv('metaData_StuderExperiments.csv', header=0, index_col=0)

# Load up all TFs to use
inputTfs = []
for i in ['GSE20573','GSE26867','GSE32658','GSE45223','GSE51533']:
    for j in subsets[i].keys():
        with open(i+'/correlated_tf_regulators_'+i+'.csv','r') as inFile:
            header = inFile.readline().strip().split(',')
            #print header
            while 1:
                inLine = inFile.readline()
                if not inLine:
                    break
                splitUp = dict(zip(header,inLine.strip().split(',')))
                inputTfs += splitUp['MEME_'+j].split(';')+splitUp['WEEDER_'+j].split(';')+splitUp['TFBS_DB_'+j].split(';')

inputTfs = [i for i in list(set(inputTfs)) if not i=='']

# Calculate correlation coefficient for all
rhos = {}
pvalues = {}
for tf in inputTfs:
    for pert in list(metaData.columns.values)[2:]:
        for gse in subsets:
            for subset in subsets[gse]:
                if not pert+'_'+gse+'_'+subset in rhos:
                    rhos[pert+'_'+gse+'_'+subset] = {}
                    pvalues[pert+'_'+gse+'_'+subset] = {}
                surrogate = metaData[[pert]].loc[subsets[gse][subset]]
                if int(tf) in list(expression[gse].columns.values) and not (int(np.sum(surrogate)==0) or int(np.sum(surrogate))==int(surrogate.size)):
                    tfExp = expression[gse][[int(tf)]].loc[subsets[gse][subset]]
                    tmp = pearsonr(tfExp,surrogate)
                    rhos[pert+'_'+gse+'_'+subset][tf] = tmp[0][0]
                    pvalues[pert+'_'+gse+'_'+subset][tf] = tmp[1][0]
                    #print (tf, pert, gse, subset, rhos[pert+'_'+gse+'_'+subset][tf], pvalues[pert+'_'+gse+'_'+subset][tf])
                else:
                    rhos[pert+'_'+gse+'_'+subset][tf] = np.nan
                    pvalues[pert+'_'+gse+'_'+subset][tf] = np.nan
                    #print (tf, gse, subset, rhos[pert+'_'+gse+'_'+subset][tf], pvalues[pert+'_'+gse+'_'+subset][tf])


#pd.DataFrame(results).to_csv('results.csv')
pd.DataFrame(rhos).to_csv('rhos.csv')
pd.DataFrame(pvalues).to_csv('pvalues.csv')

#Load in Gene2entrez Id
import csv 
with open('gene2entrezId.csv', mode='r') as infile:
    reader = csv.reader(infile)
    with open('gene2entrezId_new.csv', mode='w') as outfile:
        writer = csv.writer(outfile)
        mydict = {rows[1]:rows[0] for rows in reader}

#mydict = pd.read_csv('gene2entrezId.csv',header=0, index_col=0)


# Plot
import seaborn as sns

sns.set(color_codes=True)
dg=pd.DataFrame(rhos)
db=dg.fillna(0)
df = db.transpose()
#df = db.dropna()
#df=dg.fillna(0)
#df = df.drop(df.index[[0, 1, 2, 3, 7, 8, 9, 10, 13, 17, 20, 21, 22, 24, 27, 29, 31, 33, 34, 35, 36, 37, 38, 39, 42, 43, 44, 45, 46, 47, 48, 49, 50, 51, 54, 55, 56, 57, 58, 59, 60, 61, 62, 66, 67, 68, 69, 70, 72, 73, 74, 75, 76, 77, 78, 79, 80, 81, 82, 83, 86, 89, 93, 97, 100, 103, 106, 110, 114]], axis=0)
df = df.loc[[ind1 for ind1 in df.index if sum(df.loc[ind1]!=0)>0]]
keepers = []
rho_cutoff = 0.8
for tf1 in df.columns:
    if max(abs(df[tf1]))>=rho_cutoff:
        keepers.append(tf1)

df2 = df[keepers]

df1 = [mydict[entrez] for entrez in df2.keys() if entrez in mydict]

e1 = [entrez for entrez in df2.keys() if entrez in mydict]

df4 = df2[e1]
df4.columns = df1
#df4=df2
##### Dropped 'all' experiment series, duplicates for combos, KSR and N2.
df4 = df4.drop(df4.index[[0, 1, 2, 3, 7, 8, 9, 10, 13, 17, 20, 21, 22, 24,
                          27, 29, 31, 33, 34, 35, 36, 37, 38, 39, 42, 43,
                          44, 45, 46, 47, 48, 49, 50, 51, 54, 55, 56, 57, 
                          58, 59, 60, 61, 62, 66, 67, 68, 69, 70, 72, 73,
                          74, 75, 76, 77, 78, 79, 80, 81, 82, 83, 84, 85,
                          86, 87, 88, 89, 90, 91, 92, 93, 94, 95, 96, 97,
                          98, 99, 100, 101, 102, 103, 104, 105, 106, 107,
                          108, 109, 110, 111, 112, 113, 114, 115, 116]], axis=0)
##### Renamed Combos 
df4.rename(index = {'Ascorbic_acid_GSE32658_l': 'Combo1_GSE32658_l', 
                     'Ascorbic_acid_GSE32658_lsfc':'Combo1_GSE32658_lsfc',
                     'Ascorbic_acid_GSE32658_lsf':'Combo1_GSE32658_lsf',
                     'F8_GSE32658_lsf':'Combo2_GSE32658_lsf', 
                     'F8_GSE32658_lsfc':'Combo2_GSE32658_lsfc', 
                     'SB431542_GSE45223_nc':'Combo3_GSE45223_nc',
                     'SB431542_GSE45223_be':'Combo3_GSE45223_be',
                     'SB431542_GSE26867_lsb':'Combo3_GSE26867_lsb',
                     'SB431542_GSE26867_lsb3i':'Combo3_GSE26867_lsb3i',
                     'DAPT_GSE26867_lsb3i':'Combo4_GSE26867_lsb3i',
                     'GDNF_GSE26867_lsb':'Combo5_GSE26867_lsb',
                     'GDNF_GSE26867_lsb3i':'Combo5_GSE26867_lsb3i'}, 
                                 inplace = True) 
##### Changed row order
custom_dict = {'LDN193189_GSE32658_l':0, 'LDN193189_GSE32658_lsf':1, 
               'LDN193189_GSE32658_lsfc':2, 'SB431542_GSE32658_l':3, 
               'SB431542_GSE32658_lsf':4, 'SB431542_GSE32658_lsfc':5,
               'Noggin_GSE51533_pip':6, 'Combo3_GSE26867_lsb':7,
               'Combo3_GSE26867_lsb3i':8,  'Combo3_GSE45223_nc':9,
               'Combo3_GSE45223_be':10, 'Combo4_GSE26867_lsb3i':11,
               'CHIR99021_GSE32658_lsfc':12, 'CHIR99021_GSE45223_nc':13,
               'CHIR99021_GSE45223_be':14, 'BMP4_EDN3_GSE45223_be':15,
               'SHH_C25II_GSE20573_nsbs':16, 'Combo2_GSE32658_lsf':17,
               'Combo2_GSE32658_lsfc':18, 'Combo5_GSE26867_lsb':19,
               'Combo5_GSE26867_lsb3i':20,'Combo1_GSE32658_l':21,
               'Combo1_GSE32658_lsf':22, 'Combo1_GSE32658_lsfc':23}
df4 = df4.iloc[df4.index.map(custom_dict).argsort()]
#df2.keys() = df1
#rhos.keys()
#[i.split('_')[1] for i in rhos.keys()]
#gseNums = [j for i in rhos.keys() for j in i.split('_') if 'GSE' in j]
#gseNums = [[j for j in i.split('_') if j.find('GSE')==0][0] for i in df4.index]
##### Added row colorbar based on treatment
gseNums = [i.split('_')[0] for i in df4.index]
#len(rhos.keys())
#def unique(lut):
   # unique_list = [] 
    
    # traverse for all elements 
    #for x in lut: 
        # check if exists in unique_list or not 
        #if x not in unique_list: 
            #unique_list.append(x) 
    # print list 
    #for x in unique_list: 
        #print (x)

lut = dict(zip(list(set(gseNums)),["red", "orange", "yellow", "green", "blue", "purple", "pink", "gray", "saddlebrown", "aqua", "black"]))

#df3 = dict(zip(np.unique(lut), ["red", "orange", "yellow", "green", "blue", "purple", "brown", "silver", "pink"]))
row_colors = [lut[i] for i in gseNums]

g = sns.clustermap(df4.fillna(0), cmap = "bwr", method = 'average',row_cluster =False , figsize = (100,60), linecolor = 'grey', row_colors=row_colors)

g.savefig('test_avg_row.pdf',format='pdf')

# Surrogate to Surrogate correlation
# All GSE
metaData1 = metaData.corr().drop(index='Day', columns='Day')
sa = sns.clustermap(metaData1, cmap ="RdYlBu_r")

sa.savefig('Surrogate_correlation.pdf',format='pdf')

# GSE20573
Subset20 = metaData[metaData['GSE']== 'GSE20573'] 

Subset20 = Subset20.corr()
Subset20 = Subset20.fillna(0)


Subset20 = Subset20.loc[:, (Subset20 != 0).any(axis=0)]
Subset20 = Subset20[(Subset20.T != 0).any()]
Subset20 = Subset20.drop(index='Day', columns='Day')
s = sns.clustermap(Subset20, cmap ="RdYlBu_r")
s.savefig('Surrogate_correlation_20.pdf',format='pdf')

# GSE26867
Subset26 = metaData[metaData['GSE']== 'GSE26867'] 

Subset26 = Subset26.corr()
Subset26 = Subset26.fillna(0)


Subset26 = Subset26.loc[:, (Subset26 != 0).any(axis=0)]
Subset26 = Subset26[(Subset26.T != 0).any()]
Subset26 = Subset26.drop(index='Day', columns='Day')
s1 = sns.clustermap(Subset26, cmap ="RdYlBu_r")
s1.savefig('Surrogate_correlation_26.pdf',format='pdf')

# GSE32658
Subset3 = metaData[metaData['GSE']== 'GSE32658'] 

Subset3 = Subset3.corr()
Subset3= Subset3.fillna(0)


Subset3 = Subset3.loc[:, (Subset3 != 0).any(axis=0)]
Subset3 = Subset3[(Subset3.T != 0).any()]
Subset3 = Subset3.drop(index='Day', columns='Day')
s2 = sns.clustermap(Subset3, cmap ="RdYlBu_r")
s2.savefig('Surrogate_correlation_3.pdf',format='pdf')

# GSE45223

Subset4 = metaData[metaData['GSE']== 'GSE45223'] 

Subset4 = Subset4.corr()
Subset4= Subset4.fillna(0)


Subset4 = Subset4.loc[:, (Subset4 != 0).any(axis=0)]
Subset4 = Subset4[(Subset4.T != 0).any()]
Subset4 = Subset4.drop(index='Day', columns='Day')
s3 = sns.clustermap(Subset4, cmap ="RdYlBu_r")
s3.savefig('Surrogate_correlation_4.pdf',format='pdf')

# GSE51533

Subset5 = metaData[metaData['GSE']== 'GSE51533'] 

Subset5 = Subset5.corr()
Subset5= Subset5.fillna(0)


Subset5 = Subset5.loc[:, (Subset5 != 0).any(axis=0)]
Subset5 = Subset5[(Subset5.T != 0).any()]
Subset5 = Subset5.drop(index='Day', columns='Day')
s4 = sns.clustermap(Subset5, cmap ="RdYlBu_r")
s4.savefig('Surrogate_correlation_5.pdf',format='pdf')