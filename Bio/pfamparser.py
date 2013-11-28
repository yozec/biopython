#!/usr/bin/env python

import sys
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

if len(sys.argv)<3:
	print "\n\t**Pfamparser takes 2 arguments: pfamseq.move i swisspfam**\n\t**You can download it from ftp://ftp.sanger.ac.uk/pub/databases/Pfam/releases/ **\n\t**After You choosen Pfam version, you sholud download pfamseq.gz and swisspfam.gz **\n"
	exit(1)

f_seq = open(sys.argv[1])
lines = f_seq.readlines()
f_seq.close()

proteins_hash={}


for i in xrange(len(lines)):
	
	line=lines[i].split()
	
	if line[0][0]=='>':
		
		sequence=''
		
		for j in xrange(i+1, len(lines)):
			if lines[j][0]!='>':
				sequence=sequence+lines[j].rstrip()
			else:

				protein=SeqRecord(sequence)
				protein.id=line[0][1:]
				protein.description=lines[i][15:].rstrip()
				proteins_hash[protein.id]=protein
				break

	

f = open(sys.argv[2])

lines = f.readlines()
f.close()


for i in xrange(0, len(lines)):

	line=[]
	
	line=lines[i].split()
	
	if len(line)>0:
		
		if line[0][0]=='>' and line[0][1:] in proteins_hash:
			protein_id=line[0][1:]
			coordinate=[]

			for j in xrange(i+1, len(lines)):
							
				if lines[j]!='\n':

					line_domain=lines[j].split()

					ile_domain=int(line_domain[1])
					
					for y in xrange(ile_domain):		
						coordinate_one_domain=line_domain[-1-y].split('-')
						
						for z in xrange(2):
							coordinate_one_domain[z]=int(coordinate_one_domain[z])
							coordinate.append(coordinate_one_domain[z])

						domain=SeqRecord(proteins_hash[protein_id].seq[coordinate_one_domain[0]-1:coordinate_one_domain[1]])
						domain.id=line_domain[0]
						domain.features=coordinate_one_domain
						proteins_hash[protein_id].annotations[domain.seq]=domain

				else:
					i=j+1
					break



#print proteins_hash['P74763'].protein_sequence
#print proteins_hash['11S3_HELAN'].annotation['MASKATLLLAFTLLFATCIARHQQRQQQQNQCQLQNIEALEPIEVIQAEAGVTEIWDAYDQQFQCAWSILFDTGFNLVAFSCLPTSTPLFWPSSREGVILPGCRRTYEYSQEQQFSGEGGRRGGGEGTFRTVIRKLENLKEGDVVAIPTGTAHWLHNDGNTELVVVFLDTQNHENQLDENQRRFFLAGNPQAQAQSQQQQQRQPRQQSPQRQRQRQRQGQGQNAGNIFNGFTPELIAQSFNVDQETAQKLQGQNDQRGHIVNVGQDLQIVRPPQDRRSPRQQQEQATSPRQQQEQQQGRRGGWSNGVEETICSMKFKVNIDNPSQADFVNPQAGSIANLNSFKFPILEHLRLSVERGELRPNAIQSPHWTINAHNLLYVTEGALRVQIVDNQGNSVFDNELREGQVVVIPQNFAVIKRANEQGSRWVSFKTNDNAMIANLAGRVSASAASPLTLWANRYQLSREEAQQLKFS'].domain.id
