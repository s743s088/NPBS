# This program expects the input file to be a single string: 
#  'Transcript 1 name sequence Transcript 2 name sequence ...etc..'
#This is easiliest done in a text editor or in bash. Open your assembly and remove all line breaks.

#This script will format your assembly so it can be used with the bootstrapping code.
#Output with be a new fasta file of the assembly with sequentially number transcript ids and a spreadsheet that keeps track of the original name.
#As with all my code, you might have to adjust it a bit to work accurately with your assembly.

#Use this new assembly when mapping with bowtie2

f = open('INPUT FILE','r')

list=[]

for idx, row in enumerate(f):
	
	list=row.replace(' ','').split('>')

list = list[1:]

list2 = []

for i in range(len(list)):
	
	list2.append(list[i].replace(']', '] ').split(' '))

g = open('OUTPUT FILE WITH NEW IDS','w')
h = open('SPREADSHEET OF NEW AND ORIGINAL IDS','w')

h.write('New Id, Original Id,\n')

for i in range(len(list2)):
	
	g.write('>NEW_ID_PREFIX_%g\n' % (i+1))
	g.write('%s\n' % list2[i][1])
	
	h.write('NEW_ID_PREFIX_%g, %s,\n' % ((i+1), list2[i][0]))

