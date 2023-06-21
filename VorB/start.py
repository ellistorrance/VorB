#VorB
#environment location: /opt/anaconda3/envs/VorB
	#conda install python=3.7 (diamond needs 3.7)
	#python ver: /opt/anaconda3/envs/VorB/bin/python/
	#conda install -c bioconda diamond
	#additional reqs: numpy, matplotlib, seaborn, pandas (graphic1.py)
#VorB stands for "Virus or Bacteria?" and is a program to detect virus-like bacterial machinery
#uses outputs from geNomad on user-specified reference strain for detection in a broader bacterial species 

#overall - standalone version, version which includes genomad, potentially more "heavy" version with more bells and whistles (graphic generation, evolving database of known phage-like machinery for tentative annotation)
#call vorb 
#input = output folder of genomad

#step 1 
#read output of *genomic_virus_genes.tsv to determine if there is a literal mention of capsid - case insensitive 

#python step1.py -in ./out -query ./Pseudomonas_simiae
#out path is in genomad folder unless specified by user.
#python step1.py -in /Users/ellistorrance/Documents/Roux_Lab/Tailocin/Pipeline/VorB_Tests/Pseudomonas_aeruginosa/Pseudomonas_aeruginosa_100_GCA_000531435.1_PA38182_genomic -out /Users/ellistorrance/Documents/Roux_Lab/Tailocin/Pipeline/VorB_Tests/Pseudomonas_aeruginosa -query /Users/ellistorrance/Documents/Roux_Lab/Tailocin/Pipeline/VorB_Tests/Pseudomonas_aeruginosa/Genes

#python start.py -in /Users/ellistorrance/Documents/Roux_Lab/Tailocin/Pipeline/VorB_Tests/Vibrio_splendidus/Vibrio_splendidus_GCA_000152765.1_ASM15276v1_genomic -out /Users/ellistorrance/Documents/Roux_Lab/Tailocin/Pipeline/VorB_Tests/Vibrio_splendidus -query /Users/ellistorrance/Documents/Roux_Lab/Tailocin/Pipeline/VorB_Tests/Vibrio_splendidus/genes

# VorB requires:
# a reference strain to be analyzed by genomad
# nucleotide CDS of a population of related strains refered to as "query" population - file format = .fa or .fna

#possible option to translate prior to alignment for less related individuals? 


import os
import sys

from datetime import datetime
startTime = datetime.now()

inputs = sys.argv

#help output 
	#add colors :)

if '-h' in inputs:
	print("\nVorB Required Arguments:")
	print("-in <path to genomad output folder>")
	print("-query <path to related genomes> # putative virus or bacterial genes anotated by GeNomad will be queried against these genomes to look for evidence of orthology which would indicate phage domestication and/or phage-like bacterial machinery. ") 
	print("\nVorB Optional Arguments:")
	print("-out <path to VorB output> #Default VorB output is within specified geNomad output dir")
	print("-restart #if querying a previous GeNomad+VorB analysis on a new set of genomes you may want to use this option to avoid re-making DIAMOND database")


	#maybe: option for figure design


#parse input path 
path = ''
if "-in" in inputs:
	i = inputs.index("-in")
	path = inputs[i+1]
	if path[-1] == "/":
		pass
	else:
		path += "/"
elif not "-in" in inputs: #exit if no input given
	print("Error: User Must Specify Path to geNomad output with: '-in <path_to_genomad_output>'")
	quit()

#parse output path if specified
out_path = path # out_path == input path unless otherwised specified by user
if "-out" in inputs:
	i = inputs.index("-out")
	out_path = inputs[i+1]
	if out_path[-1] == "/":
		pass
	else:
		out_path += "/"

#parse output path if specified
script_path = __file__.strip('start.py')


##########################################
#CHECK STEP:
#Check user has met all requirements to run VorB

# Step 1:
# read virus summary from geNomad to determine if capsid has been annotated 
# capsid not yet known to be present in any incidence of phage-like bacterial machinery. 
# having a capsid indicates the virus is an active or degraded prophage.
# not having a capsid may indicate the virus has been degraded or is actually phage-like or phage-derived bacterial machinery. 
###########################################
print('\n#########################')
dict = {} # dictionary summarizing info needed for building database of gene fragments found to be viral by genomad 
header = '' # this is what the user chose to call genomad output files. 
genomad = os.listdir(path)
for folder in genomad:
	if folder.endswith('summary'):
		header += folder.strip('_summary')

		#if nothing was found in genomad: tell user and exit. 
		if os.stat(path + header + '_summary/' + header + '_virus_genes.tsv').st_size == 0:
			print("GeNomad did not detect any viral sequences for VorB to analyze:/n" + path + header + '_summary/' + header + '_virus_genes.tsv may be empty/n' + "Exiting Application.")
			quit()
		
		#look through file to find viruses with annotated capsid

		dict['virus'] = {}
		with open(path + header + '_summary/' + header + '_virus_genes.tsv', 'r') as f:
			for line in f:
				if 'capsid' in str.lower(line): # if capsid appears, get name of virus
					name = "_".join((line.split() [0]).split("_", 3)[:-1])
					if not name in dict['virus']:
						#virus.append(name)
						dict['virus'][name] = {}
						print(name + ' has an annotated capsid and is \x1B[3mnot\x1B[0m likely to be phage-like bacterial machinery.')
						#probably not phage-like bacterial machinery due to having capsid, but check for domestication/degredation. Don't exclude.				
		f.close()

		#look through file again to determine if there is any other virus(es) that was not found to have a capsid. 

		dict['unknown'] = {}
		with open(path + header + '_summary/' + header + '_virus_genes.tsv', 'r') as f:
			for line in f:
				if not "_".join((line.split() [0]).split("_", 3)[:-1]) in dict['virus']:
					if not "_".join((line.split() [0]).split("_", 3)[:-1]) in dict['unknown'] and not "_".join((line.split() [0]).split("_", 3)[:-1]) == '':
						name = "_".join((line.split() [0]).split("_", 3)[:-1])
						dict['unknown'][name] = {}
						print(name + ' may be virus-like bacterial machinery or a degraded prophage.')
		f.close()

		#if nothing was found in genomad, exit. Redundant. Necessary?
		if dict['virus'] == 0 and dict['unknown'] == 0:
			print("GeNomad did not detect any viral sequences for VorB to analyze:/n" + path + header + '_summary/' + header + '_virus_genes.tsv may be empty/n' + "Exiting Application.")
			quit()
print('#########################')
query = ''
if "-query" in inputs:
	i = inputs.index("-query")
	query = inputs[i+1]
	if query[-1] == "/":
		pass
	else:
		query += "/"
elif not "-query" in inputs: #exit if no input given
	print("Error: User Must Specify Path Dataset to Query With: '-query <path to folder containing CDS files>'")
	quit()

#check step
cds = os.listdir(query)
for file in cds:
	f = open(query + file, 'r')
	for line in f:
		if not line.startswith('>'):
			if 'e' in line or 'E' in line:
				print(line)
				print('Issue Found: all files in ' + query + ' must be predicted nucleotide CDS in fasta format. Files containing amino acid sequences were found.')
				print('Exiting VorB')
				quit()
	f.close()

#move forward to create diamond database if genomad found virus-like thingoids. 
if not dict['virus'] == 0 and not dict['unknown'] == 0:
	if not "-restart" in inputs:
		print('\n Continuing to Step2: Creating Diamond Database...')
		#print('python step2.py ' + path + ' ' + header + ' ' + out_path + ' ' + query)
		os.system('python '+ script_path +'step2.py ' + path + ' ' + header + ' ' + out_path + ' ' + query)
		os.system('python '+ script_path +'step3.py ' + path + ' ' + header + ' ' + out_path + ' ' + query)
		os.system('python '+ script_path +'step4.py '+ path + ' ' + header + ' ' + out_path + ' ' + query)
		os.system('python '+ script_path +'graphic1.py '+ path + ' ' + header + ' ' + out_path + ' ' + query)
		
	####add restart argument -see step2 for list of folders to remove

if "-restart" in inputs:
		print('\nRestart Indicated: Removing all VorB documents from specified path...')
		os.system('rm -r ' + out_path + 'db')
		os.system('rm -r ' + out_path + 'blast')
		os.system('rm -r ' + out_path + 'alignments')
		os.system('rm -r ' + out_path + 'concatenates')
		os.system('rm -r ' + out_path + 'graphics')
		os.system('rm -r ' + out_path + 'summary')

		print('#########################')
		print('\nContinuing to Step2: Creating Diamond Database...')
		#print('python step2.py ' + path + ' ' + header + ' ' + out_path + ' ' + query)
		os.system('python '+ script_path +'step2.py ' + path + ' ' + header + ' ' + out_path + ' ' + query)
		os.system('python '+ script_path +'step3.py ' + path + ' ' + header + ' ' + out_path + ' ' + query)
		os.system('python '+ script_path +'step4.py '+ path + ' ' + header + ' ' + out_path + ' ' + query)
		os.system('python '+ script_path +'graphic1.py '+ path + ' ' + header + ' ' + out_path + ' ' + query)




