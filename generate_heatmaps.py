#prep the environment and load libraries
import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt

import sys
import time
#

#tracking steps and time
print('Starting run (' + time.strftime("%d %b %Y %H:%M:%S", time.localtime()) + ')')

#read in our depth files and metadata
mapping_data = pd.read_csv('all_depths.txt',sep='\t')
meta_data = pd.read_csv('metadata.txt',sep='\t')

#merged metadata into the main table
mapping_data = mapping_data.merge(meta_data, on='SampleID')

#mapping_data = mapping_data[mapping_data['SampleID'] == sample]

#apply a mapping to generate internal normalization as a column
#also force cast columns into types for safety

#re-assign variable types so they aren't all strings
mapping_data['Position'] = mapping_data['Position'].astype(int)
mapping_data['Depth'] = mapping_data['Depth'].astype(float)

#Replace any zeros with NaN (this is to fix log(0) errors)
mapping_data['Depth'] = mapping_data['Depth'].replace(0, np.nan)
#Generate LogDepths and replace NaNs with zeros
mapping_data['LogDepth'] = np.log10(mapping_data['Depth']).replace(np.nan, 0)
mapping_data['Depth'] = mapping_data['Depth'].replace(np.nan, 0)

#Create groupins by SampeId to allow for ordering later
group1 = mapping_data.groupby('SampleID').Depth
group2 = mapping_data.groupby('SampleID').LogDepth

#Normalize both Depth and Log Depth
mapping_data['NormalizedDepth'] = 100 * (mapping_data.Depth - group1.transform('min'))/group1.transform(np.ptp)
mapping_data['LogNormalizedDepth'] = 100 * (mapping_data.LogDepth - group2.transform('min'))/group2.transform(np.ptp)

mapping_data['NormalizedDepth'] = mapping_data['NormalizedDepth'].astype(float)
mapping_data['LogNormalizedDepth'] = mapping_data['LogNormalizedDepth'].astype(float)

print('Finished Normalizing (' + time.strftime("%d %b %Y %H:%M:%S", time.localtime()) + ')')

#Create an ordering by percent mapping
ordering = mapping_data.groupby('SampleID').Order.agg('mean').sort_values().index

#need to pivot the table to put it in the format seaborn wants
mapping_data['SampleID'] = pd.Categorical(mapping_data['SampleID'], ordered=True, categories=ordering)

#Generate the heatmap data which requires a pivoted table
all_heatmap_data = pd.pivot_table(mapping_data, values='LogNormalizedDepth', index=['SampleID'], columns='Position')

#Setup the basic heatmap plotting components
f, ax = plt.subplots(figsize=(20,15))
sns.heatmap(all_heatmap_data, mask=all_heatmap_data<1, vmin=0, vmax=100, cmap='Blues')

plt.xticks((0,1000000,2000000), (0,1000000,2000000))

plt.savefig('all_samples.png',dpi=300)

#Draw verticles lines at the known 16S coordinates
ax.vlines([co-ordinates of sites of interest],0,1,color="red",transform=ax.get_xaxis_transform())

#print('Plot Generated for ' + sample + ' (' + time.strftime("%d %b %Y %H:%M:%S", time.localtime()) + ')')
print('Plot 1 Generated and Saved (' + time.strftime("%d %b %Y %H:%M:%S", time.localtime()) + ')')

#save figure as a PNG
plt.savefig('all_samples_marked.png',dpi=300)

print('Plot 2 saved ('  + time.strftime("%d %b %Y %H:%M:%S", time.localtime()) + ')')
