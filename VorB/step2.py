# VorB Step2 is finding best blast hits to 40 bacterial and/or archaeal/bacterial COGs in 
# all query strains given by user and finding best blast hits to viral proteins found by 
# genomad in all strains given by user. the output is found in /VorB/alignments/ under 
# query for COG sequences and under virus for virus-like homologs found in query strains.
# feeds step3.py.

#This script uses the tool DIAMOND:
#Buchfink B, Reuter K, Drost HG, "Sensitive protein alignments at tree-of-life scale using DIAMOND", Nature Methods 18, 366–368 (2021). doi:10.1038/s41592-021-01101-x
# make the virus database. 
# blast query seqs against virus DB and COG DB


#Reads: os.system('python step2.py ' + path + ' ' + header + ' ' + out_path + ' ' + query)
#Example: python step2.py ./out/ GCA_000963495.1_ASM96349v1_genomic ./VorB/ ./Pseudomonas_simiae/

import os 
import sys
import math 

path = sys.argv[1]
header = sys.argv[2]
out_path = sys.argv[3]
query = sys.argv[4]

try:
	os.mkdir(out_path + 'db') # if restart = FALSE , remove
except:
	pass

try:
	os.mkdir(out_path + 'blast') # if restart = FALSE , remove
except:
	pass

try:
	os.mkdir(out_path + 'blast/virus/') # if restart = FALSE , remove
except:
	pass

try:
	os.mkdir(out_path + 'blast/query_set/') # if restart = FALSE , remove
except:
	pass

try:
	os.mkdir(out_path + 'alignments/') # if restart = FALSE , remove
except:
	pass

try:
	os.mkdir(out_path + 'alignments/query_set/') # if restart = FALSE , remove
except:
	pass

try:
	os.mkdir(out_path + 'alignments/virus/') # if restart = FALSE , remove
except:
	pass

## make the virus database. 
script_path = __file__.strip('step2.py')
os.system('diamond makedb --in ' + path + header + '_summary/' + header + '_virus_proteins.faa -d ' + out_path + 'db/' +  header + '_genomad --quiet')

###################
# blast viral faa seqs from genomad against phrog database. incorporate summary data into SOMETHING :
# might be useful for defining whether target is absolutely a phage and for further defining 
# subtype of bacteriocin/phage-like bacterial machinery if not a phage

try:
	os.mkdir(out_path + 'blast/phrog/') # if restart = FALSE , remove
except:
	pass

####################

#### blast virus DB against each seq query and blast COG db against each seq in query
#blasting identified vir proteins against nucleotide CDS from query set 
queries = os.listdir(query)
strains = []
for q in queries:
	cds = ((os.path.splitext(query + q)[0]).split('/') [-1])
	strains.append(cds)
	#### blast virus DB against each seq query.
	os.system('diamond blastx -d ' + out_path + 'db/' +  header + '_genomad  -q ' + query + q +  ' -o ' + out_path + 'blast/virus/' + cds + '_virus --quiet')
	#### blast bacteria DB against each seq query.
	os.system('diamond blastx -d '+ script_path + 'bacteria_archaea_COGS.dmnd  -q ' + query + q +  ' -o ' + out_path + 'blast/query_set/' + cds + '_Markers --quiet')





#print(strains)

print('Preparing query COG sequences for alignment...')
cogs = ['COG0359', 'COG0097', 'COG0149', 'COG0233', 'COG0018', 'COG0088', 'COG0094', 'COG0198', 'COG0305', 'COG0036', 'COG0099', 'COG0526', 'COG0126', 'COG0166', 'COG0571', 'COG0195', 'COG0537', 'COG0096', 'COG0533', 'COG0480', 'COG0006', 'COG0563', 'COG0013', 'COG0197', 'COG0092', 'COG1136', 'COG0495', 'COG0492', 'COG0484', 'COG0275', 'COG0202', 'COG0264', 'COG0536', 'COG0636', 'COG0242', 'COG0227', 'COG0522', 'COG0706', 'COG0250', 'COG0222', 'COG1132']
markers = os.listdir(out_path + 'blast/query_set/')
dict = {} #cog name and gene id 
seqs = {} # file name, gene id, sequence
for co in cogs:
	dict[co] = []
for b in markers:
	file = open(out_path + 'blast/query_set/' + b)
	temp = {}
	fn = b.strip('_Markers')
	seqs[fn] = {}
	for line in file:
		try:
			c = line.split() [1] # COG0526
			gene = line.strip('lcl|').split() [0] # AVQG01000001.1_cds_ERH61266.1_168
			bit = line.split() [10]
			if not c in temp:
				temp[c] = []
				temp[c].append([gene,bit])
			else:
				temp[c].append([gene,bit])
		except:
			pass
	#print(temp)
	file.close()
	for c in temp:
		#dict[c] = [] # cog, file_name&gene_name
		if len(temp[c])==1:
			name = b.strip('_Markers') + '&' + temp[c] [0][0]
			seqs[fn][temp[c] [0][0]] = ''
			#print(name) GCA_021606485.1_ASM2160648v1_cds_from_genomic&WKCL01000009.1_cds_MCF5340842.1_5605
			dict[c].append(name)
		#find the best bit score , small numbers keep reverting to zero in >/< equations and so had to use log in math function to avoid
		elif len(temp[c]) > 1:
			#print(temp[c]) #[['WKCH01000028.1', '1.54e-123'], ['WKCH01000004.1', '7.06e-05'], ['WKCH01000042.1', '1.33e-40']]
			x = 0
			char = ''
			for item in temp[c]:
				if float(item [1]) == 0.0:
					x = 100
					char = item[0]
					break
				else:
					if x == 0:
						x += math.log10(float(item[1])) 
					#x == math.log10(float(item[1]))
						char = item[0]
					elif x > math.log10(float(item[1])):
						char = item[0]
			name = b.strip('_Markers') + '&' + char
			seqs[fn][char] = ''
			dict[c].append(name)
#print(seqs)
#get sequences of found COGS for bacteria. move them to cog file for alignment. right now, these are untranslated. This is fine except if user is querying distantly related taxa. May want to add translation option.
#queries = os.listdir(query)
#for name in seqs:
#	if len(seqs[name]) > 1: # don't open files with 0 found COGS
#		#queries = os.listdir(query)
#		#print(queries)
#		#print(len(seqs[name]))

# do same as above for virus sets 
# folder number = to number of viruses identified. within each folder "COGs" will be viral genes instead
# rewrite because ran into problem where "real virus" has best hit in tailocin complex - artificially inflating the number of virus present:
# now, will cycle through column 1 first and when host gene mentioned >1 times: host gene will be "given" to virus with the best hit

print('Preparing virus homologs found in query sequences for alignment...')
dict2 = {} #vir gene name and gene id 
comp1 = {}
comp = {}
vir_full = []
summary = open(path + header + '_summary/'+ header + '_virus_genes.tsv', 'r')
for line in summary:
	if not line.startswith('gene'):
		vir_prot = line.split() [0] #CP005975.1|provirus_515927_531946_483
		vir_head = vir_prot.split('_') [0] + '_' + vir_prot.split('_') [1] + '_' + vir_prot.split('_') [2] #CP005975.1|provirus_515927_531946
		vir_full.append(vir_head)
		dict2[vir_prot] = []
		comp1[vir_prot] = []
		comp[vir_prot] = []
		#print(vir_prot)
		#print(vir_head)
summary.close()
seqs2 = {} # file name, gene id, sequence
markers = os.listdir(out_path + 'blast/virus/')
for b in markers:
	file = open(out_path + 'blast/virus/' + b)
	temp = {}
	temp1 = {}
	fn = b.strip('_virus') #file name # GCA_000465595.1_GS_De_Novo_Assembly_cds_from_genomic
	seqs2[fn] = {}
	for line in file:
		try:
			c = line.split() [1] # CP005975.1|provirus_4105422_4140707_3819 # virus hit 
			gene = line.strip('lcl|').split() [0] # JAEKCQ010000018.1_cds_MBJ2232224.1_2316 # host gene hit 
			bit = line.split() [10]
			#if not gene in temp1
			if not gene in temp1:
				temp1[gene] = []
				temp1[gene].append([c,bit])
			else:
				temp1[gene].append([c,bit])
			
			#IF NOT C IN TEMP
			if not c in temp:
				temp[c] = []
				temp[c].append([gene,bit])
			else:
				temp[c].append([gene,bit])
		except:
			pass
	file.close()
	for gene in temp1:
		if len(temp1[gene])==1:
			name = fn + '&' + gene
			#seqs2[fn][gene] = ''
			comp1[temp1[gene] [0][0]].append(name)

#		#find the best bit score , small numbers keep reverting to zero in >/< equations and so had to use log in math function to avoid
		elif len(temp1[gene]) > 1:
			#print(temp1[gene]) #[['WKCH01000028.1', '1.54e-123'], ['WKCH01000004.1', '7.06e-05'], ['WKCH01000042.1', '1.33e-40']]
			x = 0
			char = ''
			for item in temp1[gene]:
				if float(item [1]) == 0.0:
					x = 100
					char = gene
					break
				else:
					if x == 0:
						x += math.log10(float(item[1])) 
						char = gene
					elif x > math.log10(float(item[1])):
						char = gene
			name = fn + '&' + char
			#seqs2[fn][char] = ''
			comp1[temp1[char] [0][0]].append(name)
	for c in temp:
		#dict2[c] = [] # cog, file_name&gene_name
		if len(temp[c])==1:
			name = b.strip('_virus') + '&' + temp[c] [0][0]
			#seqs2[fn][temp[c] [0][0]] = ''
			#print(name) GCA_021606485.1_ASM2160648v1_cds_from_genomic&WKCL01000009.1_cds_MCF5340842.1_5605
			comp[c].append(name)
		#find the best bit score , small numbers keep reverting to zero in >/< equations and so had to use log in math function to avoid
		elif len(temp[c]) > 1:
			#print(temp[c]) #[['WKCH01000028.1', '1.54e-123'], ['WKCH01000004.1', '7.06e-05'], ['WKCH01000042.1', '1.33e-40']]
			x = 0
			char = ''
			for item in temp[c]:
				if float(item [1]) == 0.0:
					x = 100
					char = item[0]
					break
				else:
					if x == 0:
						x += math.log10(float(item[1])) 
					#x == math.log10(float(item[1]))
						char = item[0]
					elif x > math.log10(float(item[1])):
						char = item[0]
			name = b.strip('_virus') + '&' + char
			#seqs2[fn][char] = ''
			comp[c].append(name)

leave_out = []
for vir in comp:
	for gene in comp[vir]:
		if not gene in comp1[vir]:
			leave_out.append(gene)
	for gene in comp1[vir]:
		if not gene in comp[vir]:
			leave_out.append(gene)
for vir in comp:
	for gene in comp[vir]:
		for gene in comp1[vir]:
			if not gene in leave_out:
				fn = gene.split('&') [0]
				char = gene.split('&') [1]
				if not char in seqs2[fn]:
					seqs2[fn][char] = ''
				if not gene in dict2[vir]:
					dict2[vir].append(gene)

			
			


#print(comp) #'CP005975.1|provirus_515927_531946_483': ['GCA_022637535.1_ASM2263753v1_cds_from_genomic&CP093333.1_cds_UNK67599.1_1132', 'GCA_900111895.1_IMG-taxon_2663762772_annotated_assembly_cds_from_genomic&FOKB01000002.1_cds_SFB04354.1_804', 'GCA_000934565.1_ASM93456v1_cds_from_genomic&CP010896.1_cds_AJP50845.1_1136', 'GCA_016307615.1_ASM163076...
#print(comp1) #{'CP005975.1|provirus_515927_531946_483': ['GCA_022637535.1_ASM2263753v1_cds_from_genomic&CP093333.1_cds_UNK67599.1_1132', 'GCA_900111895.1_IMG-taxon_2663762772_annotated_assembly_cds_from_genomic&FOKB01000002.1_cds_SFB04354.1_804', 'GCA_000934565.1_ASM93456v1_cds_from_genomic&CP01
#print(seqs2) # {'GCA_022637535.1_ASM2263753v1_cds_from_genomic': {'CP093333.1_cds_UNK66981.1_497': '', 'CP093333.1_cds_UNK67575.1_1108': '', 'CP093333.1_cds_UNK67576.1_1109': '', 'CP093333.1_cds_UNK67577.1_1110': '', 'CP093333.1_cds_UNK67579.1_1112': '', 'CP093333.1_cds_UNK67580.1_1113': '', 'CP093333.1_cds_UNK67581.1_1114': '', 'CP093333.1_cds_UNK67582.1_1115': '', 'CP093333.1_cds_UNK67583.1_1116': ...
#print(dict2) # 'CP005975.1|provirus_4105422_4140707_3820': ['GCA_022637535.1_ASM2263753v1_cds_from_genomic&CP093333.1_cds_UNK68408.1_1972', 'GCA_900111895.1_IMG-taxon_2663762772_annotated_assembly_cds_from_genomic&FOKB01000001.1_cds_SFA83703.1_3933', 'GCA_900111895.1_IMG-taxon_2663762772_annotated_assembly_cds_from_genomic&FOKB01000013.1_cds_SFB57017.1_4884',...


#
#
for file in queries:
	name = ((os.path.splitext(query + file)[0]).split('/') [-1])
	tmp = {}
			#print(file)
	f = open(query + file, 'r')
	for line in f:
		if line.startswith('>'):
			head = line.strip('\n').strip('>').split() [0]
			tmp[head] = ''
		else:
			tmp[head] += line.strip('\n')
	f.close()
	for head in tmp:
		if head in seqs[name]:
			seqs[name][head] += tmp[head] #best COG homologs for each query strain
		if head in seqs2[name]:
			seqs2[name][head] += tmp[head] #best virus protein homologs for each query strain
#print(seqs2) # VIRUS HOMOLOGS {'GCA_022637535.1_ASM2263753v1_cds_from_genomic': {'CP093333.1_cds_UNK66981.1_497': 'ATGAACGATATTTCAGCTCCTGAGCAATACGATCTGCAAACCGCTGCCCTGAAGGTGCCGCCGCATTCCATCGAGGCCGAACAGGCTGTGCTCGGTGGTTTGATGCTGGACAACAACGCCTGGGAACGCGTGCTGGATCAAG...
#print(seqs) # COGS {'GCA_021605975.1_ASM2160597v1_cds_from_genomic': {'WKDM01000001.1_cds_MCF5184599.1_40': 'ATGATCCATATCCCCCGAGCGGAATACACCCGACGCCGCAAGGCGCTCATGGCGCAGATGGAACCCAACAGCATCGCGATCCTGCCGGCCGCTGCCGTGGCGA...
#print(dict) # COGS & genomes found in {'COG0359': ['GCA_021605975.1_ASM2160597v1_cds_from_genomic&WKDM01000010.1_cds_MCF5186723.1_398', 'GCA_021606505.1_ASM2160650v1_cds_from_genomic&WKCH01000003.1_cds_MCF5344800.1_2646', 'GCA_0143580...
for cog in dict:
	if not len(dict[cog]) == 0:
		out = open(out_path + 'alignments/query_set/' + cog + '.fna' , 'w')
		#out.write(str(dict[cog]) + '\n')
		#out.write(str(len(dict[cog])) + '\n')
		for name in dict[cog]:
			out.write('>' + name + '\n')
			thing = name.split('&') [0]
			char = name.split('&') [1]
			out.write(seqs[thing][char] + '\n')
		out.close()
#
#
#
#
#print(vir_full) # ['CP005975.1|provirus_515927_531946', 'CP005975.1|provirus_515927_531946', 'CP005975.1|provirus_515927_531946', 'CP005975.1|provirus_515927_531946', 'CP005975.1|provirus_515927_531946', 'CP005975.1|provirus_515927_531946', 'CP005975.1|provirus_515927_531946', 'CP005975.1|provirus_515927_531946', 'CP005975.1|provirus_515927_531946', 'CP005975.1|provirus_515927_531946', 'CP005975.1|provirus_515927_531946'
for v in vir_full:
	try:
		os.mkdir(out_path + 'alignments/virus/' + v.replace('|','_'))
	except:
		pass
	for prot in dict2:
		if prot.startswith(v):
			if not len(dict2[prot]) == 0: # some alignments have only one sequence? keep the way it is to make note or maybe lower threshhold to some percent of query seqs?
				out = open(out_path + 'alignments/virus/' + v.replace('|','_') + '/' + prot.replace('|','_') + '.fna' , 'w')
				for name in dict2[prot]:
					out.write('>' + name + '\n')
					thing = name.split('&') [0]
					char = name.split('&') [1]
					out.write(seqs2[thing][char] + '\n')
				out.close()
###print('python step3.py ' + path + ' ' + header + ' ' + out_path + ' ' + query)
#os.system('python step3.py ' + path + ' ' + header + ' ' + out_path + ' ' + query)