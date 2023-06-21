#VorB step3 is aligning all genes found in the previous step (both query COGs and
# viral homologs to viruses detected with genomad) and building aligned gene concatenates for
#query seqs and viral gene sets for each virus found.


#This Script uses the tool MAFFT:
#Katoh and Toh (Bioinformatics 23:372-374, 2007) PartTree: an algorithm to build an approximate tree from a large number of unaligned sequences (describes the PartTree algorithm).
#Katoh, Kuma, Toh and Miyata (Nucleic Acids Res. 33:511-518, 2005) MAFFT version 5: improvement in accuracy of multiple sequence alignment (describes [ancestral versions of] the G-INS-i, L-INS-i and E-INS-i strategies)
#Katoh, Misawa, Kuma and Miyata (Nucleic Acids Res. 30:3059-3066, 2002) MAFFT: a novel method for rapid multiple sequence alignment based on fast Fourier transform (describes the FFT-NS-1, FFT-NS-2 and FFT-NS-i strategies)


#Reads: os.system('python step3.py ' + path + ' ' + header + ' ' + out_path + ' ' + query)
#Example: python step3.py ./out/ GCA_000963495.1_ASM96349v1_genomic ./VorB/ ./Pseudomonas_simiae/

import os 
import sys

path = sys.argv[1]
header = sys.argv[2]
out_path = sys.argv[3]
query = sys.argv[4]

##align query set
print('Aligning Marker Genes for Query Set...')
q = os.listdir(out_path + 'alignments/query_set/')
for item in q:
	if not item.endswith('.align'):
		st = ((os.path.splitext(out_path + 'alignments/query_set/' + item)[0]).split('/') [-1])
		os.system('mafft --quiet ' + out_path + 'alignments/query_set/' + item + ' > ' + out_path + 'alignments/query_set/' + st + '.align')

##align virus set
print('Aligning Viral Gene Homologs Found in Query Set...')
v = os.listdir(out_path + 'alignments/virus/')
for virus in v:
	p = os.listdir(out_path + 'alignments/virus/' + virus + '/')
	for prot in p:
		if not prot.endswith('.align'):
			st = ((os.path.splitext(out_path + 'alignments/virus/' + virus + '/' + prot)[0]).split('/') [-1])
			os.system('mafft --quiet ' + out_path + 'alignments/virus/' + virus + '/' + prot + ' > ' + out_path + 'alignments/virus/' + virus + '/' + st + '.align')


# build concatenates 
print('Building Multi-gene Concatenates...')
try:
	os.mkdir(out_path + 'concatenates/')
except:
	pass

total = 0
names = []
seqs = os.listdir(query)
for file in seqs:
	total += 1
	let = ((os.path.splitext(query + file)[0]).split('/') [-1])
	names.append(let)
#print(total)
#print(names)
#print(names)

### Concat for query COGS
concat = {}
for name in names:
	concat[name] = []

dict = {}
for cog in q:
	if cog.endswith('.align'):
		file = open(out_path + 'alignments/query_set/' + cog, 'r')
		dict[cog] = {}
		for line in file:
			if line.startswith('>'):
				name = line.strip('>').split('&') [0]
				dict[cog][name] = []
			else:
				if len(line) == 61 :
					dict[cog][name].append(line)
				if len(line) < 61 : #make end trailing lines of each gene equal in length
					line = line.strip('\n')
					N = 60 - len(line)
					K = '-'
					line = line.ljust(N + len(line), K) + '\n'
					dict[cog][name].append(line)
					#print(line)
		file.close()
#print(dict)
for cog in dict:
	x = 0
	for name in dict[cog]:
		x = len(dict[cog][name])
		concat[name].append(''.join(dict[cog][name]))
	for name in concat: #in instances where a COG was not found in a strain, this inputs an empty alignment sequence '-----' etc. for that COG for a given strain. 
		if name not in dict[cog]:
			empty = []
			while len(empty) < x:
				K = '-'
				empty.append(''.ljust(60,K) + '\n')
			concat[name].append(''.join(empty))

marker_out = open(out_path + 'concatenates/query_markers.align', 'w')
for name in concat:
	marker_out.write('>' + name + '\n' + ''.join(concat[name]))
marker_out.close()

	
#### Concat for virus gene homologs found in query strains
concat = {}
dict = {}
for virus in v:
	dict[virus] = {}
	concat[virus] = {}
	for name in names:
		concat[virus][name] = []
	q = os.listdir(out_path + 'alignments/virus/' + virus + '/')
	for cog in q:
		if cog.endswith('.align'):
			file = open(out_path + 'alignments/virus/' + virus + '/' + cog, 'r')
			dict[virus][cog] = {}
			for line in file:
				if line.startswith('>'):
					name = line.strip('>').split('&') [0]
					dict[virus][cog][name] = []
				else:
					if len(line) == 61 :
						dict[virus][cog][name].append(line)
					if len(line) < 61 : #make end trailing lines of each gene equal in length
						line = line.strip('\n')
						N = 60 - len(line)
						K = '-'
						line = line.ljust(N + len(line), K) + '\n'
						dict[virus][cog][name].append(line)
						#print(line)
			file.close()
#print(dict)
#print(dict)
#print(concat)
for virus in dict:
	for cog in dict[virus]:
		x = 0
		for name in dict[virus][cog]:
			x = len(dict[virus][cog][name])
			#print(x)
			concat[virus][name].append(''.join(dict[virus][cog][name]))
		for name in concat[virus]: #in instances where a COG was not found in a strain, this inputs an empty alignment sequence '-----' etc. for that COG for a gi
			if name not in dict[virus][cog]:
				#print(cog + ' ' + name)
				empty = []
				while len(empty) < x:
					K = '-'
					empty.append(''.ljust(60,K) + '\n')
				concat[virus][name].append(''.join(empty))
#print(concat)
for virus in concat:
	v_marker_out = open(out_path + 'concatenates/' + virus + '.align', 'w')
	for name in concat[virus]:
		v_marker_out.write('>' + name + '\n' + ''.join(concat[virus][name]))
	v_marker_out.close()


##### make tree - may need to be an optional parameter
#os.system('python step_tree.py ' + path + ' ' + header + ' ' + out_path + ' ' + query)
#
#print('python step4.py '+ path + ' ' + header + ' ' + out_path + ' ' + query)
#os.system('python step4.py '+ path + ' ' + header + ' ' + out_path + ' ' + query)
