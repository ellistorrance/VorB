#step4 is a summary step:
#number of total genes found in each seq for each virus 
#avg PI


#Reads: os.system('python step4.py ' + path + ' ' + header + ' ' + out_path + ' ' + query)
#Example: python step4.py ./out/ GCA_000963495.1_ASM96349v1_genomic ./VorB/ ./Pseudomonas_simiae/

import os 
import sys
import math

path = sys.argv[1]
header = sys.argv[2]
out_path = sys.argv[3]
query_in = sys.argv[4]

script_path = __file__.strip('step4.py')
#print(script_path)
#print('Initating gene annotation with PHROG database using HHsuite v3.3.0...')


try:
	os.mkdir(out_path + 'summary/')
except:
	pass

print('Summarizing Identity Characteristics...')

#############################################################
def nucleotide_PI(list_of_seqs):

	nucs = ['a','t','c','g','A','T','C','G']
	y=1
	pairwise = []
	for seqA in list_of_seqs [0:len(list_of_seqs)-1]:
		for seqB in list_of_seqs [y:len(list_of_seqs)]:
			x = 0
			total = 0
			same = 0
			while x < len(seqA):
				if seqA [x] in nucs and seqB [x] in nucs:
					total+=1
					if seqA [x] == seqB[x]:
						same+=1
				x+=1

			#print(str(same) + ' ' + str(total))
			try:
				pairwise.append(float(same/total))
			except:
				pairwise.append(float(0.0))
			 
		y+=1
	#deal with alignment pairs which might be all gaps: 
	if len(pairwise) != len(list_of_seqs):
		x = len(list_of_seqs) - len(pairwise)
		m = 0
		while m < x:
			pairwise.append(float(0.0))
			m+=1
	
	avg = float(sum(pairwise))/float(len(pairwise))


	#avg = np.mean(pairwise, dtype=np.float64)
	#med = np.median(pairwise)
	#stdev = np.std(pairwise, dtype=np.float64)
	#del(pairwise)
	#out_list = [avg,med,stdev]
	return round(avg,3)

################################################################

def GC(list_of_seqs):

	nucs = ['a','t','c','g','A','T','C','G']
	all_gc = []
	for seqA in list_of_seqs [0:len(list_of_seqs)]:
		total = 0
		gc = 0
		x = 0
		while x < len(seqA):
			if seqA [x] in nucs:
				total+=1
				if seqA [x] == 'g' or seqA [x] == 'c':
					gc+=1
			x+=1
		all_gc.append(float(gc/total))
	try:  
		avg = float(sum(all_gc))/float(len(all_gc))
		return round(avg,3)
	except:
		return 'NA'

################################################################
query = {}
query_pi = {}
q = os.listdir(out_path + 'alignments/query_set/')
for cog in q:
	if cog.endswith('.align'):
		query[cog] = {}
		query_pi[cog] = {}
		#if cog.endswith('.align'):
		file = open(out_path + 'alignments/query_set/' + cog,'r')
		x = 0 
		for line in file:
			if line.startswith('>'):
				x+=1
				name = line.strip('>').strip('\n')
				query[cog][name]= []
			else:
				query[cog][name].append(line.strip('\n'))
		file.close()
		tmp = []
		for name in query[cog]:
			tmp.append(''.join(query[cog][name]))
		query_pi[cog][x] = {}
		#print(cog + ' ' + str(nucleotide_PI(tmp)))
		query_pi[cog][x][str(nucleotide_PI(tmp))] = ''
		query_pi[cog][x][str(nucleotide_PI(tmp))] += str(GC(tmp))
#print(query_pi)
out = open(out_path + 'summary/' + 'query_marker_summary.tsv', 'w')
out.write('marker_ortholog_ID' + '\t' + 'number of query strains marker found in' + '\t' + 'average nucleotide PI' + '\t' + 'average GC content' + '\n')
for cog in query:
	for x in query_pi[cog]:
		for y in query_pi[cog][x]:
			out.write(str(cog).strip('.align') + '\t' + str(x) + '\t' + str(y)+ '\t' + str(query_pi[cog][x][y]) + '\n')
out.close()
del(query)
del(query_pi)

q1 = os.listdir(out_path + 'alignments/virus/')
virus = {}
virus_pi = {}
for folder in q1:
	virus[folder] = {}
	virus_pi[folder] = {}
	q = os.listdir(out_path + 'alignments/virus/' + folder)
	for cogs in q:
		if cogs.endswith('.align'):
			cog = cogs.strip('.align')
			virus[folder][cog] = {}
			virus_pi[folder][cog] = {}

			file = open(out_path + 'alignments/virus/' + folder + '/' + cogs,'r')
			num = 0
			for line in file:
				if line.startswith('>'):
					name = line.strip('>').strip('\n')
					virus[folder][cog][name] = []
					num+=1
				else:
					virus[folder][cog][name].append(line.strip('\n'))
			virus_pi[folder][cog][num] = {}
			file.close()
			tmp = []
			for name in virus[folder][cog]:
				tmp.append(''.join(virus[folder][cog][name]))
				#print(tmp)
			virus_pi[folder][cog][num][str(nucleotide_PI(tmp))] = ''
			virus_pi[folder][cog][num][str(nucleotide_PI(tmp))] += str(GC(tmp))
del(virus)
#print(virus)
#print(virus_pi)


#get all gene names for virus
v_dict = {}
gene_sum = open(path + header + '_summary/' + header + '_virus_genes.tsv', 'r')
for line in gene_sum:
	if not line.startswith('gene'):
		vir = line.replace('|','_')
		vir = vir.split('\t') [0]
		gene = vir.split('_') [0] + '_' + vir.split('_') [1] + '_' + vir.split('_') [2] + '_'+ vir.split('_') [3] + '_'+ vir.split('_') [4]
		vir = vir.split('_') [0] + '_' + vir.split('_') [1] + '_' + vir.split('_') [2] + '_'+ vir.split('_') [3]
		annot = line.strip('\n').split('\t') [18]
		if not vir in v_dict:
			v_dict[vir] = {}
			v_dict[vir][gene] = ''
			v_dict[vir][gene] += annot
		else:
			v_dict[vir][gene] = ''
			v_dict[vir][gene] += annot
gene_sum.close()
			
#print(v_dict)


##################
#initiate phrog stuff
os.system('python ' + script_path + 'step_phrog.py ' + path + ' ' + header + ' ' + out_path + ' ' + query_in)

#open phrog output files and get best hit from first line
ph = {}
phrogs = os.listdir(out_path + 'blast/phrog/')
for file in phrogs:
	if file.endswith('.faa.tsv'):
		with open(out_path + 'blast/phrog/' + file) as f:
			first_line = f.readline()
			id = (first_line.split() [1]).strip('phrog_')
			name = first_line.split() [0]
			ph[name] = {}
			ph[name][id] = ''



annot = open(script_path +'phrog_annot_v4.tsv', 'r')
for line in annot:
	n = line.split('\t') [0]
	for head in ph:
		if n in ph[head]:
			ph[head][n]+=line.strip('\n')
annot.close()

			






v= 1
for folder in v_dict:
	out = open(out_path + 'summary/' + folder + '_summary.tsv', 'w')
	out.write('gene_ID' + '\t' + 'plot_ID' +'\t' + 'geNomad_annotation' + '\t' + 'phrog_ID' + '\t' + 'phrog_color' + '\t' + 'phrog_annotation' + '\t' + 'phrog_category' + '\t' + 'number of query strains gene found in' + '\t' + 'average nucleotide PI' + '\t' + 'average GC content' + '\n')
	#print(folder)
	for cog in virus_pi[folder]:
		gene = str(cog).strip('.align')
		for num in virus_pi[folder][cog]:
			for x in virus_pi[folder][cog][num]:
				try:
					for text in ph[gene]:
						out.write(gene + '\t' + 'v' + str(v) +'\t' + str(v_dict[folder][gene])+ '\t' + 'phrog_' + str(ph[gene][text])  + '\t' + str(num) + '\t' +  str(x) + '\t' +  str(virus_pi[folder][cog][num][x]) + '\n')
				except:
					out.write(gene + '\t' + 'v' + str(v) +'\t' + str(v_dict[folder][gene])+ '\t' + 'NA' + '\t' + 'NA' + '\t' + 'NA' + '\t' + 'NA'  + '\t' + str(num) + '\t' +  str(x) + '\t' +  str(virus_pi[folder][cog][num][x]) + '\n')
	out.close()
	v+=1
#

print('Creating identity graphic...')
#print('python graphic1.py '+ path + ' ' + header + ' ' + out_path + ' ' + query_in)
#os.system('python graphic1.py '+ path + ' ' + header + ' ' + out_path + ' ' + query_in)



