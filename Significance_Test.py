## This code will look through BS files and assess significance
## using both the average log2 fold change in expression value 
## from the BS files and the confidence intervals they create. 

## Import Libraries ##
import gc
import os
from math import log as log
from numpy import mean as mn

## Define Functions ##

def get_Counts():
	
	cmd1 = os.popen('wc -l OUTPUT_FROM_Preprocess_Map_file.py_FOR_FIRST_SAMPLE').readline().split(' ')
	numRR_Sample1 = float(cmd1)
	
	cmd2 = os.popen('wc -l FILE_THAT_CONTAINS_ONLY_SAM_HEADER').readline().split(' ')
	numGenes = int(cmd2)
	
	cmd3 = os.popen('wc -l OUTPUT_FROM_Preprocess_Map_file.py_FOR_SECOND_SAMPLE').readline().split(' ')
	numRR_Sample2 = float(cmd3)
	
	return numRR_Sample1, numGenes, numRR_Sample2
	

def normalization():
	
	f = open('FILE_THAT_CONTAINS_ONLY_SAM_HEADER','r')
	
	GeneLengths = []
	
	for idx,row in enumerate(f):
		a = row.replace('\n','').replace('LN:','').split('\t')
		
		GeneLengths.append(int(a[2]))
	
	return GeneLengths
	

def init_Gene_Count_list(sim):
	
	Sample2_BS_Counts = [0]*numGenes
	
	for i in xrange(len(Sample2_BS_Counts)):
		Sample2_BS_Counts[i] = []
	
	gc.disable()
	
	BSfile = open('Sample2_BS_FILE_PREFIX%g.csv' % sim, 'r')
	
	L1 = []
	
	for idx, row in enumerate(BSfile):
	
		L1 = row.replace('\n','').split(',')
		
		for i in xrange(len(L1)-1):
			
			Sample2_BS_Counts[idx].append(int(L1[i])/numRR_Sample2/GeneLengths[idx])
		
	
	Sample1_BS_Counts = [0]*numGenes
	
	for i in xrange(len(Sample1_BS_Counts)):
		Sample1_BS_Counts[i] = []
	
	gc.disable()
	
	BSfile = open('Sample1_BS_FILE_PREFIX%g.csv' % sim, 'r')
	
	L2 = []
	
	for idx, row in enumerate(BSfile):
	
		L2 = row.replace('\n','').split(',')
		
		for i in xrange(len(L2)-1):
			
			Sample1_BS_Counts[idx].append(int(L2[i])/numRR_Sample1/GeneLengths[idx])
	
	gc.enable()
	
	return Sample2_BS_Counts, Sample1_BS_Counts

def Get_Log_Change():##THIS SECTION WILL REQUIRE SOME EDITING DEPENDING ON HOW MANY SAMPLE TYPES YOU HAVE IN THE 'all_counts.csv' file
	f = open('csv_SUMMARY_FILE','r')
	
	log_FOLD_change=[]
	
	for idx, row in enumerate(f):
		
		if idx > 0:
			
			index = [int(row.split(',')[0].split('_')[-1])-1,0,0,0]
			
			Sample1 = float(row.split(',')[5])
			Sample2 = float(row.split(',')[7])
			
			if Sample1 != 0.0 and Sample2 != 0.0:
				if abs(log(Sample1/Sample2,2)) >=foldChange:
					log_FOLD_change.append(index)
	
	return log_FOLD_change

def sig_test_function(sim):
	
	gc.disable()
	
	print 'Testing for significant differential expression of each transcript between samples in simulation %g of %g.' % (sim+1, numSIM)
	
	for i in xrange(len(log_FOLD_change)):
		Sample2_BS_Counts[log_FOLD_change[i][0]].sort()
		
		highCrit1 = Sample2_BS_Counts[log_FOLD_change[i][0]][-1]
		lowCrit1 = Sample2_BS_Counts[log_FOLD_change[i][0]][0]
		average1 = mn(Sample2_BS_Counts[log_FOLD_change[i][0]])
		
		Sample1_BS_Counts[log_FOLD_change[i][0]].sort()
		
		highCrit2 = Sample1_BS_Counts[log_FOLD_change[i][0]][-1]
		lowCrit2 = Sample1_BS_Counts[log_FOLD_change[i][0]][0]
		average2 = mn(Sample1_BS_Counts[log_FOLD_change[i][0]])
		
		if (lowCrit2 > highCrit1) and highCrit1 != 0.0:
			
			if log(average2/average1,2) >= foldChange:
				log_FOLD_change[i][1]+=1
			else:
				log_FOLD_change[i][3]+=1
		
		elif (highCrit2 < lowCrit1) and highCrit2 != 0.0:
			
			if log(average2/average1,2) <= -foldChange:
				log_FOLD_change[i][2]+=1
			else:
				log_FOLD_change[i][3]+=1
		else:
				log_FOLD_change[i][3]+=1
	
	gc.enable()
	
	return 

def make_output_files():
	
	rawfile = open('OUT_ASSEMBLY_FILE_PREFIX_FROM_Format_Assembly.py_.fasta', 'r')
	
	edfile = open('OUT_ASSEMBLY_FILE_FROM_REFIX_Format_Assembly.py_ed.fasta', 'w')
	
	for idx, row in enumerate(rawfile, start=1):
		g = row.replace('\n', ',')
		edfile.write(g)
	
	edfile = open('OUT_ASSEMBLY_FILE_FROM_REFIX_Format_Assembly.py_ed.fasta', 'r')
	x = 0
	
	for idx, row in enumerate(edfile, start=1):
		
		g = row.split(',')
	
	print 'Creating files with list of differentially expressed genes.'
	
	f = open('RESULTS.csv', 'w')

	f.write('Gene ID, #UP, #DOWN, #NS, Gene Seq., \n')
	
	for i in xrange(len(log_FOLD_change)):
		
		f.write('NEW_ID_PREFIX_FROM_Format_Assembly.py_FILE_%g, %g, %g, %g, %s, \n' % (log_FOLD_change[i][0]+1, log_FOLD_change[i][1], log_FOLD_change[i][2], log_FOLD_change[i][3], g[((log_FOLD_change[i][0])*2)+1]))
		
	f.flush()
	f.close()
	
	return

def make_fasta_files():
	#Creates fasta files of genes found significant in a specified number of runs
	#These files can be the input into another program that runs a blast search on all the genes
	
	print 'Creating .fasta files of genes found significantly up and down regulated in %g of %g independent runs.' % (minSig, numSIM)
	
	f = open('RESULTS.csv', 'r')
	
	g = open('UP_reg_RESULTS.fasta', 'w')
	
	h = open('DOWN_reg_RESULTS.fasta', 'w')
		
	for idx, row in enumerate(f):
		
		l = row.replace(' ', '').replace('\n',',').split(',')
		
		if idx>0:
			if int(l[1]) >= minSig:
				g.write('>%s\n%s\n' % (l[0], l[4]))
			elif int(l[2]) >= minSig:
				h.write('>%s\n%s\n' % (l[0], l[4]))
	
	f.close()
	g.close()
	return

## MAIN ##

## Define Globals
numSIM = 10
minSig = 10 #Number of times a gene must be significant in all simulations to be reported in the fasta files
numBS = 100
foldChange = 2.0

(numRR_Sample1, numGenes, numRR_Sample2) = get_Counts()

GeneLengths = normalization()

log_FOLD_change = Get_Log_Change()

for sim in xrange(numSIM):
	
	(Sample2_BS_Counts, Sample1_BS_Counts) = init_Gene_Count_list(sim)
	
	sig_test_function(sim)

make_output_files()

make_fasta_files()

print 'Program complete.'

