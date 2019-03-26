import glob
import sys 

if len(sys.argv) != 2:
	print('usage: python add_geneinfo.py file')
	sys.exit()

genes = open('../data/genes.txt')

mapping = {}

for line in genes:
	line = line.split('\t')
	mapping[line[0]] = {'gene_id': line[0], 'gene_symbol': line[3], 'gene_type': line[2], 'go_term_name': line[5], 'go_term_definition': line[6]}


file = sys.argv[1]
output = open(file[:-4] + '.geneinfo.tsv', 'w')
file = open(file)

firstline = file.readline().strip()
output.write(firstline + '\tgene_type\tgo_term_name\tgo_term_definition\n')

for line in file:
	output.write(line.strip())
	line = line.strip().split('\t')
	gene = line[firstline.split('\t').index('ensembl_gene_id')]
	if gene in mapping:
		if len(mapping[gene]['go_term_name']) == 0:
			output.write('\t' + mapping[gene]['gene_type'] + '\t---\t---\n')
		else:
			output.write('\t' + mapping[gene]['gene_type'] + '\t' + mapping[gene]['go_term_name'] + '\t' + mapping[gene]['go_term_definition'] + '\n')
	else:
		output.write('\t---\t---\t---\n' )

output.close()

