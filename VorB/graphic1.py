#make dotplot comparing # of strains gene appears in to average identity across identified gene set. This will be done for host strain marker genes and all potential phage genes. 


#Example: python graphic1.py ./out/ GCA_000963495.1_ASM96349v1_genomic ./VorB/ ./Pseudomonas_simiae/
import matplotlib.pyplot as plt
import os 
import sys
import numpy as np
import pandas as pd
import seaborn as sns
from collections import Counter

path = sys.argv[1]
header = sys.argv[2]
out_path = sys.argv[3]
query = sys.argv[4]

try:
	os.mkdir(out_path + 'graphics/')
except:
	pass

data = os.listdir(out_path + 'summary/')
name = []
x_data = []
y_data = []
cat = []
name2 = []
#king = []
dict1 = {}
dict12 = {}
for d in data:
	if d.endswith('.tsv'):
		if d == 'query_marker_summary.tsv':
			file1 = open(out_path + 'summary/' + d , 'r')
			for line in file1: 
				if not line.startswith('marker'):
					id  = 'Host_COG'
					strains = int(line.split('\t') [1])
					ident = int(float(line.split('\t') [2])*100)
					#king.append((id,strains,ident))
					name.append(id)
					y_data.append(strains)
					x_data.append(ident)
					#cat.append('NA')
			file1.close()
		else:
			file = open(out_path + 'summary/' + d , 'r')
			for line in file:
				#print(line)
				#try:
				if not line.startswith('gene_ID'):
					#id  = d.strip('_summary.tsv')
					id = line.split('\t') [1]
					strains = int(line.split('\t') [7])
					ident = int(float(line.split('\t') [8])*100)
					cats = line.split('\t') [6]
					name.append(id)
					name2.append(id)
					y_data.append(strains)
					x_data.append(ident)
					cat.append(cats)

					#king.append((id,strains,ident))
				#except Exception as e:
					#print(line + ' ' + str(e))
				#	pass
			file.close()
dict1['File ID'] = name
dict1['Number of Strains with Homolog to Gene'] = y_data
dict1['Paiwise Nuc. Identity (%) of Gene'] = x_data
items = int(len(Counter(name).keys())/2)







#dict1['cat'] = cat # vir categories
m = max(x_data)
#print(m)
df = pd.DataFrame(dict1)

dict12['File ID'] = name2
dict12['Phrog Category'] = cat
df2 = pd.DataFrame(dict12)

#########################
#make a specific color palette by File ID From : https://stackoverflow.com/questions/46173419/seaborn-change-color-according-to-hue-name
ls = []
for t in name:
	if t not in ls:
		ls.append(t)
#print(ls)
palette = dict(zip(ls, sns.color_palette(n_colors=len(ls))))
#print(palette)


#########################

#dotplot of # of individuals vs Pairwise identity
#graphic attribution: https://www.python-graph-gallery.com/43-use-categorical-variable-to-color-scatterplot-seaborn
#ax = sns.scatterplot( x="x", y="y", data=df, hue='name', legend=True, s=12)
#ax = sns.relplot(data=df, x="x", y="y", hue='name', legend=False, col = 'name')
sns.set_theme(style='darkgrid', rc={'figure.dpi': 147},		  
			  font_scale=0.7)
ax = sns.jointplot(data=df, x="Paiwise Nuc. Identity (%) of Gene", y="Number of Strains with Homolog to Gene", hue="File ID", palette=palette)
#ax.set(xlabel='% Average Pairwise Identity Across Gene', ylabel='# of Query Strains Gene Found In')
## Move the legend to an empty part of the plot
#plt.xlim(0,100)
#plt.ylim(1,max(y_data))
#sns.move_legend(ax, "lower center",bbox_to_anchor=(.5, 1), ncol=1, title=None, frameon=False,)
#plt.tight_layout()
leg = plt.legend(frameon=False,loc='lower center', ncol=items)
leg.set_in_layout(False)
plt.savefig(out_path + 'graphics/strains_vs_identity_dotplot.png', bbox_inches="tight")   # save the figure to file
plt.close()	# close the figure window


import textwrap
def wrap_labels(ax, width, break_long_words=False):
	labels = []
	for label in ax.get_yticklabels():
		text = label.get_text()
		labels.append(textwrap.fill(text, width=width,
					  break_long_words=break_long_words))
	ax.set_yticklabels(labels, rotation=0)
#https://seaborn.pydata.org/generated/seaborn.jointplot.html
sns.set_theme(style='darkgrid', rc={'figure.dpi': 147},			  
			  font_scale=0.7)
ax=sns.countplot(data=df2, y="Phrog Category", hue="File ID", palette=palette)
wrap_labels(ax, 30)
ax.figure
sns.move_legend(ax, "lower center",bbox_to_anchor=(.5, 1), ncol=items, title=None, frameon=False,)
plt.tight_layout()
plt.savefig(out_path + 'graphics/phrog_category_barchart.png', bbox_inches="tight")   # save the figure to file
plt.close()	# close the figure window


print('#########################')
print('VorB complete:')
print('Summary Tables Available in ' + out_path + 'summary/' + '\n')
print('Summary Graphics Available in ' + out_path + 'graphics/' + '\n')

print('Have a Great Day! :)')


#display_figures(ax,plot_data)
#boxplot based on code from: https://practicaldatascience.co.uk/data-science/how-to-visualise-data-using-boxplots-in-seaborn
## Boxplot Strain per gene
#ax = sns.boxplot(x='File ID', y='Number of Strains with Homolog to Gene', data=df)
## 
### Add jitter with the swarmplot function
#ax = sns.swarmplot(x='File ID', y='Number of Strains with Homolog to Gene', data=df, color="grey")
##ax.set(xlabel='', ylabel='Number of Strains Containing Gene')
##ax.set_xticklabels(ax.get_xticklabels(), fontsize=5)
##plt.ylim(0,max(y_data)+10)
##plt.xticks(rotation = 90)
#
#plt.tight_layout()
#plt.savefig(out_path + 'graphics/query_strains_per_gene_boxplot.png')   # save the figure to file
#plt.close()	# close the figure window

#

## Boxplot Avg PI per gene
#ax = sns.boxplot(x='File ID', y='Paiwise Nuc. Identity (%) of Gene', data=df)
### Add jitter with the swarmplot function
#ax = sns.swarmplot(x='File ID', y='Paiwise Nuc. Identity (%) of Gene', data=df, color="grey")
##ax.set(xlabel='', ylabel='Avg Nucleotide PI (%)')
##ax.set_xticklabels(ax.get_xticklabels(), fontsize=5)
##plt.ylim(0,120)
##plt.xticks(rotation = 90)
#
#plt.tight_layout()
##plt.subplots_adjust(bottom = 0.75)
#plt.savefig(out_path + 'graphics/avg_PI_across_gene_boxplot.png' )   # save the figure to file
#plt.close()	# close the figure window
#




























