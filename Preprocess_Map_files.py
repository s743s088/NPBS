import sys

'''
  This script produces a text file of mapped read counts from a sam file. 
	Reads are only counted if both pais map to the same scaffold and if both
	reads are sequenced sucessfully (neither pair is composed of all N's).

	To use:
	python Preprocess_Map_files.py infile.sam

	out put:
	.txt file that is a list of the transcripts mapped to by each raw read

	out file format:
	transcript # read 1 (forward and reverse mapped to)
	transcript # read 2 (forward and reverse mapped to)

'''

sam_files = sys.argv[1:]
gene_dic = []

print "reading in sam header"

bad_seq = "N" * 100

for file in sam_files:
	print "processing file "+file
	count = 0
	with open(file,'r') as infile:
		AT1 = infile.readline()
		line1 = infile.readline()
		line2 = infile.readline()
		while line2:
			line1 = infile.readline()
			line2 = infile.readline()
			if 'SN:' in line1 and 'SN:' in line2 or '@PG' in line2:
				continue
			else:
				elements1 = line1.split('\t')
				elements2 = line2.split('\t')
				count +=2
				if count%1000 == 0:
					print "on line "+str(count)
				if bad_seq in elements1 or bad_seq in elements2:
					continue
				else:
					try:
						if elements1[0] == elements2[0] and elements1[2] == elements2[2] != '*':
							try:
								gene = elements1[2].split('_')
								print gene
								gene_dic.append(gene[3])
							except KeyError:
								pass
					except IndexError:
						pass
	with open('NAME_OF_OUTFFILE.txt', 'w') as outfile:
		print "writing out file"
		print gene_dic
		for j in gene_dic:
			outfile.write('%s\n' % j)
		outfile.flush()
		outfile.close()
