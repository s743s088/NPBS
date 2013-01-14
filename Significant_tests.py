## Differential gene expression non-parametric bootstrap program

## Import Libraries ##
import gc
import os


## Define Functions ##

def get_Counts():
  
## OUTFILES IN CMD1 AND CMD3 SHOULD OBVIOUSLY BE FROM DIFFERENT LIBRARIES TO WISH TO COMPARE
	
	cmd1 = os.popen('wc -l OUTPUT_FILE_FROM_Preprocess_Map_files.txt').readline().split(' ')
	numRR_1 = float(cmd1[1])
	
	cmd2 = os.popen('wc -l SAM_HEADER_FILE.txt').readline().split(' ')
	numGenes = int(cmd2[3])
	
	cmd3 = os.popen('wc -l OUTPUT_FILE_FROM_Preprocess_Map_files.txt').readline().split(' ')
	numRR_2 = float(cmd3[1])
	
	return numRR_1, numGenes, numRR_2
	

def normalization():
	
	f = open('SAM_Header_FILE.txt','r')
	
	GeneLengths = []
	
	for idx,row in enumerate(f):
		a = row.replace('\n','').replace('LN:','').split('\t')
		
		GeneLengths.append(int(a[2]))
	
	return GeneLengths
	

def init_Gene_Count_list(sim):
	
	2_BS_Counts = [0]*numGenes
	
	for i in xrange(len(2_BS_Counts)):
		2_BS_Counts[i] = []
	
	gc.disable()
	
	BSfile = open('BOOTSTRAP_FILES_CORRESPONDING_TO_CMD3%g.csv' % sim, 'r')
	
	L1 = []
	
	for idx, row in enumerate(BSfile):
	
		L1 = row.replace('\n','').split(',')
		
		for i in xrange(len(L1)-1):
			
			2_BS_Counts[idx].append(int(L1[i])/numRR_2/GeneLengths[idx])
		
	
	1_BS_Counts = [0]*numGenes
	
	for i in xrange(len(1_BS_Counts)):
		1_BS_Counts[i] = []
	
	gc.disable()
	
	BSfile = open('BOOTSTRAP_FILES_CORRESPONDING_TO_CMD1%g.csv' % sim, 'r')
	
	L2 = []
	
	for idx, row in enumerate(BSfile):
	
		L2 = row.replace('\n','').split(',')
		
		for i in xrange(len(L2)-1):
			
			1_BS_Counts[idx].append(int(L2[i])/numRR_1/GeneLengths[idx])
	
	gc.enable()
	
	return 2_BS_Counts, 1_BS_Counts

def init_sig_gene_storage():
	
	UP = [0]*numSIM
	DOWN = [0]*numSIM
	NS = [0]*numSIM
	
	for i in xrange(numSIM):
		UP[i] = []
		DOWN[i] = []
		NS[i] = []
	
	return UP, DOWN, NS

def sig_test_function(sim):
	
	gc.disable()
	
	print 'Testing for significant differential expression between 2zooid and gonozooid of each gene in simulation %g of %g.' % (sim+1, numSIM)
	
	for i in xrange(numGenes):
		2_BS_Counts[i].sort()
		
		highCrit1 = 2_BS_Counts[i][int(numBS*.975)]
		lowCrit1 = 2_BS_Counts[i][int(numBS*.025)]
		
		1_BS_Counts[i].sort()
		
		highCrit2 = 1_BS_Counts[i][int(numBS*.975)]
		lowCrit2 = 1_BS_Counts[i][int(numBS*.025)]
		
		if lowCrit2 > highCrit1:
			
			UP[sim].append(i+1)
			
		elif highCrit2 < lowCrit1:
			
			DOWN[sim].append(i+1)
			
		else:
			
			NS[sim].append(i+1)
	
	gc.enable()
	
	return 

def sum_over_all_sim():
	
	gc.disable()
	
	print 'Summarizing differentially expressed genes over all simulation runs.'
	
	upCOUNT = []
	downCOUNT = []
	nsCOUNT = []
	
	upTOT = []
	downTOT = []
	nsTOT = []
	
	for i in xrange(3):
		
		if i == 0:
			list1 = UP
			list2 = upTOT
			list3 = upCOUNT
		elif i == 1:
			list1 = DOWN
			list2 = downTOT
			list3 = downCOUNT
		else:
			list1 = NS
			list2 = nsTOT
			list3 = nsCOUNT
		
		for j in xrange(len(list1)):
			
			for k in xrange(len(list1[j])):
				
				list2.append(list1[j][k])
		
		list2.sort()
		
		a = 0
		c = 0
		
		for m in xrange(len(list2)):
			
			if list2[m] != list2[m-1]:
				list3.append([0,list2[m]])
				while a < len(list2):
					if list2[m] == list2[a]:
						
						a += 1
						list3[c][0] += 1
						
					else:
						c +=1
						
						break
	
	gc.enable()
	
	return upCOUNT, downCOUNT, nsCOUNT

def make_output_files():
	
	rawfile = open('TRANSCRIPTOME_FASTA_FILE', 'r')
	
	edfile = open('TRANSCRIPTOME_FASTA_FILE_EDITED', 'w')
	
	for idx, row in enumerate(rawfile, start=1):
		g = row.replace('\n', ',')
		edfile.write(g)
	
	edfile = open('TRANSCRIPTOME_FASTA_FILE_EDITED', 'r')
	x = 0
	
	for idx, row in enumerate(edfile, start=1):
		
		g = row.split(',')
	
	print 'Creating files with list of differentially expressed genes.'
	
	for k in xrange(3):
		
		if k == 0:
			f = open('UP_regulated_gene.csv', 'w')
			x = upCOUNT
		elif k == 1:
			f = open('DOWN_regulated_gene.csv', 'w')
			x = downCOUNT
		else:
			f = open('NS_regulated_gene.csv', 'w')
			x = nsCOUNT
		
		f.write('# sig, Gene ID, Gene Seq., \n')
		
		for i in xrange(len(x)):
			
			f.write('%g, TRANSCRIPT_PREFIX_%g, %s, \n' % (x[i][0], x[i][1], g[((x[i][1]-1)*2)+1]))
			
		f.flush()
		f.close()
	
	return

def make_fasta_files():
	#Creates fasta files of genes found significant in a specified number of runs
	#These files can be the input into an1 program that runs a blast search on all the genes
	
	print 'Creating .fasta files of genes found significantly up and down regulated in %g of %g independent runs.' % (minSig, numSIM)
	
	for k in xrange(3):
		
		if k == 0:
			f = open('UP_regulated_gene.csv', 'r')
			g = open('UP_reg.fasta', 'w')
		elif k == 1:
			f = open('DOWN_regulated_gene.csv', 'r')
			g = open('DOWN_reg.fasta', 'w')
		else:
			f = open('NS_regulated_gene.csv', 'r')
			g = open('NS_reg.fasta', 'w')
		
		for idx, row in enumerate(f):
			
			l = row.replace(' ', '').replace('\n',',').split(',')
			
			if idx>0:
				if k <= 1:
					if int(l[0]) >= minSig:
						g.write('>%s\n%s\n' % (l[1], l[2]))
				else:
					g.write('>%s\n%s\n' % (l[1], l[2]))
	
	f.close()
	g.close()
	return

## MAIN ##

## Define Globals
numSIM = 10
minSig = 10 #Number of times a gene must be significant in all simulations to be reported in the fasta files
numBS = 100

(numRR_1, numGenes, numRR_2) = get_Counts()

GeneLengths = normalization()

(UP, DOWN, NS) = init_sig_gene_storage()

for sim in xrange(numSIM):
	
	(2_BS_Counts, 1_BS_Counts) = init_Gene_Count_list(sim)
	
	sig_test_function(sim)

(upCOUNT, downCOUNT, nsCOUNT) = sum_over_all_sim()

make_output_files()

make_fasta_files()

print 'Program complete.'

