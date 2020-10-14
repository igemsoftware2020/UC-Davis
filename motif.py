#!/usr/bin/env python3
import random
import math



def read_memetxt(memetxt):
	'''Opens text file output of meme and extracts the sites discovered by 
	meme'''
	sites = []
	motif_stats = []
	with open(memetxt) as mt:
		while True:
			line = mt.readline()
			if line == '': break
			if line.startswith('MOTIF'):
				l = line.split()
				motifid = l[2]
				nsites = l[8]
				#skip to sequence information
				while True:
					line = mt.readline()
					if line.startswith('Sequence name'):
						line = mt.readline()
						break
				#parse through sequence information
				while True:
					line = mt.readline()
					if line.startswith('---'): break
					l = line.split()
					seq = l[0]
					if l[1] == '+' or l[1] == '-':
						strand = l[1]
						beg = int(l[2])
						pval = float(l[3])	
					else:
						strand = '+'
						beg = int(l[1])
						pval = float(l[2])
					sites.append((motifid,nsites,strand,seq,beg,pval))
				#could potentially read out info about the motifs
				while True:
					line = mt.readline()
					if line.startswith('letter-probability matrix'):
						l = line.split()
						wid = int(l[5])
						e_val = float(l[9])
						break						
		return sites				
				


def memepwm(memeouttxt): #'meme_out/meme.txt'
	'''Opens text file output of meme and extracts the position weight matrix
	 and it's corresponding information for every motif discovered by meme'''
	motifs = []
	motif_stats = []
	nummot = 1
	nt = ['A','C','G','T']
	with open(memeouttxt) as mt:
		#iterates through the whole file
		while True:
			line = mt.readline()
			if line == '': break
			#finds the beginning of the pwms
			while True:
				if line == '': break
				line = mt.readline()
				if line.startswith('letter-probability matrix'):
					m = []
					l = line.split()
					wid = int(l[5])
					nsites = int(l[7])
					e_val = float(l[9])
					for i in range(wid):
						m.append({})
						#iterates through the pwms
					for i in range(len(m)):
						line = mt.readline()
						l = line.split()
						for n in range(len(nt)):
							m[i][nt[n]] = float(l[n])
					motifs.append(m)
					bits,p_bits = score_motifbit(m)	
					motif_stats.append((f'MEME-{nummot}',wid,nsites,e_val,bits\
					,p_bits))	
					nummot += 1
	return motifs, motif_stats

#are variable names too specific?
def read_testmotif(motif_file):
	'''Opens the fasta file generated  '''
	jpositions = []
	numjsites = 0
	with open(motif_file) as tm:
		for line in tm.readlines():
			if line.startswith('>'):
				positions = [] 
				line = line.split()
				numjsites += (int((len(line[1:(len(line))]))/2))
				seq  = line[0].strip('>')
				for i in range(1,len(line)):
					if i == 1: 
						if line[1] == '[]':
							line[1] = line[1].strip("]")
						positions.append(line[1].strip("['"))
					elif i == len(line)-1:
						positions.append(line[i].strip("]'"))
					else:
						positions.append(line[i].strip("',"))
				jpositions.append((seq,positions))
	return jpositions, numjsites

	


def read_JASPAR(jasparfile):
	'''Opens a jaspar file containing the binding motif of interest and creates 
	a position weight matrix representing the binding motif'''
	motif = []
	linenum = 0
	nt = ['A', 'C', 'G', 'T']
	with open(jasparfile) as jf:
		for line in jf.readlines():
			if line.startswith('>'): continue
			line = line.split()	
			for i in range(2,len(line)-1):
				if len(motif) < (len(line)-3): motif.append({})
				motif[i-2][nt[linenum]] = int(line[i])
			linenum+= 1 
	for i in range(len(motif)):
		total = 0 
		for nt in motif[i]: total += motif[i][nt]
		for nt in motif[i]: motif[i][nt] /= total
	return motif 

def generate_site(motif):
	'''Generates a single DNA sequence from the position weight matrix '''
	motifseq = []
	total = 0
	for i in range(0,len(motif)):
		pA = motif[i]['A']
		pC = motif[i]['C']
		pG = motif[i]['G']
		pT = motif[i]['T']
		r = random.random()
		if r < pA:               motifseq += 'A'
		elif r < (pA + pC):      motifseq += 'C'
		elif r < (pA + pC + pG): motifseq += 'G' 
		else:                    motifseq += 'T'
	return motifseq
	
def new_pwm(sites):
	'''Accepts an array of binding site sequences and produces a position 
	weight matrix from them '''
	newmot =[]
	for i in range(len(sites[0])):
		A = 0
		C = 0
		G = 0
		T = 0
		newmot.append({})
		for j in range(len(sites)):
			s = sites[j][i]
			if s == 'A': A += 1
			elif s == 'C': C += 1
			elif s == 'G': G += 1
			elif s == 'T': T += 1
		newmot[i]['A'] = A/len(sites)
		newmot[i]['C'] = C/len(sites)
		newmot[i]['G'] = G/len(sites)
		newmot[i]['T'] = T/len(sites)
	return newmot
	
#meme.txt files do not find false negatives 

def pos_accuracy(mpos,jpos,mw,jw):
	'''Accepts and compares the starting positions of the meme motifs and the
	 embedded jaspar motifs to determine how far apart they are and how much
	  overlap there is   '''
	fp = 0
	fl = 0
	posdis = 0
	if jpos == '':
		fp += 1
		fl += 1
		posdis = 100
		overlap = 0
		overlap_p = 0
	elif jpos != '':
		mpos = int(mpos)
		jpos = int(jpos)
		mend = mpos + mw
		jend = jpos + jw
		if mpos >= jpos and mpos <= jend:
			posdis = mpos-jpos
			if mend >= jend:
				overlap = jend - mpos
			else:
				overlap = mend - mpos
			overlap_p = overlap/jw
			if overlap_p < 0.2:
				fp += 1
				fl += 1		
		elif mpos <= jpos and jpos <= mend:
			posdis = jpos-mpos
			if mend <= jend:
				overlap = mend - jpos
			else:
				overlap = jend - jpos
			overlap_p = overlap/jw
			if overlap_p < 0.2:
				fp += 1
				fl += 1
		else:
			fl += 1	
			posdis = 98
			overlap = 0
			overlap_p = 0
	return fp, posdis, fl, overlap, overlap_p


def global_motcompare(motif1,motif2,background): #background info needs to be a dictionary
	'''Finds global similarity score between the meme motif and the jaspar 
	motif by finding the best matching region between the two motifs and 
	filling unmatched positions with background nucleotide distribution 
	frequency  '''
	distances = []
	score = 0 
	if len(motif1)==0 or len(motif2)== 0:
		score = 0
		distances.append(score)
		max_score = 1
	elif len(motif1) == len(motif2):
		max_score = 2*len(motif1)
		for i in range(0, len(motif1)):
			d = 0
			for nt in motif1[i]:
				d += abs(motif1[i][nt] - motif2[i][nt])
			score += 2-d
		distances.append(score)		
	else:
		if len(motif1) > len(motif2):
			big = motif1
			short = motif2
		elif len(motif2) > len(motif1):
			big = motif2
			short = motif1
		max_score = 2*len(big)
		for i in range(0,len(big)-len(short)+1):
			minmotif = []
			for j in range(0,len(big)):
				minmotif.append(background)
			for k in range(0,len(short)):
				minmotif[i+k] = short[k]
			score = 0
			for l in range(len(minmotif)):
				d = 0
				for nt in minmotif[l]:
					#manhattan similarity: allows for highest score
					d += abs(minmotif[l][nt]-big[l][nt])
				score += 2-d
			distances.append(score)
	bestfit = 0
	fitindex = 0
	for i in range(0,len(distances)):
		if i == 0: 
			bestfit = distances[i]
			fitindex = i
		else:
			if distances[i] >= bestfit:
				bestfit = distances[i]
				fitindex = i
	return bestfit, bestfit/max_score
		
def local_motcompare(motif1, motif2):
	'''Finds local similarity score between the meme motif and the jaspar 
	motif by finding the best matching region between the two motifs.'''
	distances = []
	if len(motif1) == 0 or len(motif2)==0:
		score = 0
		distances.append(score)
		max_score = 1
	elif len(motif1) == len(motif2):
		max_score = 2*len(motif1)
		score = 0
		for i in range(len(motif1)):
			d = 0
			for nt in motif1[i]:
				d += abs(motif1[i][nt]-motif2[i][nt])
			score += 2-d
		distances.append(score)
	else:
		if len(motif1) > len(motif2):
			big = motif1
			short = motif2
		elif len(motif2) > len(motif1):
			big = motif2
			short = motif1
		max_score = 2*len(big)
		for i in range(0,len(big)-len(short)+1):
			score = 0
			window = big[i:i+len(short)]
			for j in range(len(window)):
				d = 0
				for nt in window[j]:
					d += abs(window[j][nt]-short[j][nt])
				score += 2-d
			distances.append(score)
	bestfit = 0
	fitindex = 0
	for i in range(len(distances)):
		if i == 0: 
			bestfit = distances[i]
			fitindex = i
		else: 
			if distances[i] > bestfit:
				bestfit = distances[i]
				fitindex = i
	return bestfit, bestfit/max_score
	
def score_motifbit(motif):
	'''Finds the informational content of the motif by summing up the information content of each individual position  '''
	score = 0
	max_score = 2*len(motif)
	for i in range(len(motif)):
		entropy = 0
		pseudos = 0
		for nt in motif[i]:
			prob = motif[i][nt]
			if prob == 0:
				continue
			entropy -= prob* math.log2(prob)
		bits = 2 - entropy - pseudos*.1
		score += bits
	if max_score > 0:
		p = score/max_score
	else:
		p = ''
	return score, p