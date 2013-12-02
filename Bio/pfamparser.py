#!/usr/bin/env python

# Copyright 2013 by xyz.  All rights reserved.
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.

import sys
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord


'''
Pfamparser needs 2 file: pfamseq.move (FASTA proteins file) and swisspfam (file with protein domains).
You can download this files from ftp://ftp.sanger.ac.uk/pub/databases/Pfam/releases/.

Script has 4 methods:

FASTAparser takes file pfamseq.move and return a dictionary: {proteinId:[proteinObject]}. The value is a list.

PFAMparser takes file swisspfam and dictionary with proteins from FASTAparse, 
create domain object, which sequence is from slicing protein sequence and add it to the proteinObject list. 
Proteins and domains are SeqReceord objects. It returns a dictionary: {proteinId:[proteinObject,[domainObject,...,domainObject].
The coordinates domains in protein sequence are in domain.feature. 

Methods proteinsViewer view single protein record and domainsViewer view all of the protein domains.


Exemple:

If you want to print record of ACON_BRAJA protein:

###print proteins_viewer(proteins_dict['ACON_BRAJA'])

and you should get:

ID: ACON_BRAJA
Name: <unknown name>
Description: 920 ACONITATE HYDRATASE (EC 4.2.1.3) (CITRATE HYDRO-LYASE) (ACONITASE).
Number of features: 0
'MTSLDSFKCKKTLKVGAKTYVYYSLPTAEK[...]ICQSFERIHRSNLVGMGVLPLTFEEGTSWSSLGLKGDEKVTLRGLVGDLKPRQKLTAEIVSGDGSLQRVSLLCRIDTLDELDYYRNGGILHYVLRKLAA'

###print domains_viewer(proteins_dict['120K_RICRI'])


ID: Pfam-B_112
Name: <unknown name>
Description: <unknown description>
Number of features: 2
'IVNATTLYAGISTLNNNQGTVTLSGGVPNTP[...]TSIETTLTLANGNIGHIVILEGAQVNTTTTGTTTIKVQDNANANFSGTQTYT'

ID: Pfam-B_112
Name: <unknown name>
Description: <unknown description>
Number of features: 2
'GTVINGKVNQTALVGGALAAGTITLDGSAT[...]DANVGSFVFNAGGTNIVSGTVGGQQGNKFNTVALENGTTVKFLGNATFNGNTTIAANSTL'

ID: Pfam-B_19744
Name: <unknown name>
Description: <unknown description>
Number of features: 2
'MVIQSANATGQVNFRHIVDVGADGTTAFK[...]VAVTNNITAIEASGAGVVQLSGTHAAELRLGNAGSIFKLAD'

'''

def FASTAparser(pfamseq):

	proteinsDict={}

	for i in xrange(len(pfamseq)):
	
		line=pfamseq[i].split()
	
		if line[0][0]=='>':
			proteinId=line[0][1:]
			proteinAnnotation=pfamseq[i][15:].rstrip()
			sequence=''
		
			for j in xrange(i+1, len(pfamseq)):
				if pfamseq[j][0]!='>':
					sequence=sequence+pfamseq[j].rstrip()
				else:
	
					protein=SeqRecord(sequence)
					protein.id=proteinId
					protein.description=proteinAnnotation
					proteinsDict[proteinId]=[protein]
					break
	return proteinsDict
	

f = open(sys.argv[2])

lines = f.readlines()
f.close()

def PFAMparser(swisspfam, protDict):

	for i in xrange(0, len(swisspfam)):

		line=[]
	
		line=swisspfam[i].split()
	
		if len(line)>0:
		
			if line[0][0]=='>' and line[0][1:] in protDict:
				proteinId=line[0][1:]
				thisProteinSequence=protDict[proteinId][0].seq
				coordinate=[]
				domainsList=[]

				for j in xrange(i+1, len(swisspfam)):
							
					if swisspfam[j]!='\n':

						domainLine=swisspfam[j].split()

						numberOfDomains=int(domainLine[1])
					
						for y in xrange(numberOfDomains):		
							coordinatesOfSingleDomain=domainLine[-1-y].split('-')
						
							for z in xrange(2):
								coordinatesOfSingleDomain[z]=int(coordinatesOfSingleDomain[z])
								coordinate.append(coordinatesOfSingleDomain[z])

							domain=SeqRecord(thisProteinSequence[coordinatesOfSingleDomain[0]-1:coordinatesOfSingleDomain[1]])
							domain.id=domainLine[0]
							domain.features=coordinatesOfSingleDomain
							domainsList.append(domain)

					else:
						protDict[proteinId].append(domainsList)
						i=j+1
						break
	return protDict


def domainsViewer(domains):
	for i in domains[1]:
		print "\n",i

def proteinsViewer(protein):
	print protein[0]



#This is default settings - files pfamparser.py, pfamseq.move and swisspfam should be in the same direcetion.
fSeq = open('pfamseq.move')
proteins = fSeq.readlines()
fSeq.close()


fSeq = open('swisspfam')
domains = fSeq.readlines()
fSeq.close()
proteinsWithDomains=PFAMparser(domains, FASTAparser(proteins))


