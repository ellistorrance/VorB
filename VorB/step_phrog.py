
##phrog with Hhsuite 3.3.0
##comes after step3 and before step4
## requires AA seqs in Genomad 
## output in blast/phrog
#
##Example: python step_phrog.py ./out/ GCA_000963495.1_ASM96349v1_genomic ./VorB/ ./Pseudomonas_simiae/
#
import os 
import sys

path = sys.argv[1]
header = sys.argv[2]
out_path = sys.argv[3]
query = sys.argv[4]

script_path = __file__.strip('step_phrog.py')
#print(script_path)
print('Initating gene annotation with PHROG database using HHsuite v3.3.0...')

try:
	os.mkdir(out_path + 'blast/phrog/')
except:
	pass

try:
	os.mkdir(out_path + 'blast/phrog/genes')
except:
	pass
#################
##################
input = path + header + '_summary/' + header + '_virus_proteins.faa'
file = open(input)
dict = {}
for line in file:
	if line.startswith('>'):
		f = line.split() [0]
		f = f.replace('|','_')
		dict[f] = ''
	else:
		dict[f] += line
file.close()
for f in dict:
	file = open(out_path + 'blast/phrog/genes/' + f.strip('>') + '.faa', 'w')
	file.write(f + '\n' + dict[f])
	file.close()

genes = os.listdir(out_path + 'blast/phrog/genes/')
for gene in genes:
	if gene.endswith('.faa'):
		os.system('hhblits -i ' + out_path + 'blast/phrog/genes/' + gene + ' -d '+ script_path + 'Phrog_HHsuite/phrogs_hhsuite_db/phrogs -n 1 -o ' + out_path + 'blast/phrog/genes/' + gene + '.hmm_out -v 0'  + ' -blasttab '+ out_path + 'blast/phrog/' + gene + '.tsv')




#### MMSEQS was trying to take 8+ hours...ain't nobody got time for that
#phrog with MMseqs2 v14.7e284
#comes after step3 and before step4
# requires AA seqs in Genomad 
# output in blast/phrog
#
#Example: python step_phrog.py ./out/ GCA_000963495.1_ASM96349v1_genomic ./VorB/ ./Pseudomonas_simiae/
#
#import os 
#import sys
##
#path = sys.argv[1]
#header = sys.argv[2]
#out_path = sys.argv[3]
#query = sys.argv[4]
##
#print('Initating gene annotation with PHROG database using MMseqs2 v14.7e284...')
##
#try:
#	os.mkdir(out_path + 'blast/phrog/')
#except:
#	pass
##
#
#input = path + header + '_summary/' + header + '_virus_proteins.faa'
#
## make mmseqsdb 
#
#
#os.system('mmseqs createdb ' + input + ' ' + out_path + 'db/' + header + '_mmseqs_db')
#os.system('mmseqs search ' + './phrogs_mmseqs_db/phrogs_profile_db ' +  out_path + 'db/' + header + '_mmseqs_db ' + out_path + 'blast/phrog/' + header + '_mmseqs_out  ./tmp -s 7')
#os.system('mmseqs createtsv ' + './phrogs_mmseqs_db/phrogs_profile_db ' +  out_path + 'db/' + header + '_mmseqs_db '+ out_path + 'blast/phrog/' + header + '_mmseqs_out ' + out_path + 'blast/phrog/' + header + '_mmseqs_output.tsv ')
#
#