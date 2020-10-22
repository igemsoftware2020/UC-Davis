import random
import math
import os
import statistics

NT = ['A', 'C', 'G', 'T']

IUPAC = {
	'A': {'A':1, 'C':0, 'G':0, 'T':0},
	'C': {'A':0, 'C':1, 'G':0, 'T':0},
	'G': {'A':0, 'C':0, 'G':1, 'T':0},
	'T': {'A':0, 'C':0, 'G':0, 'T':1},
	'R': {'A':0.5, 'C':0.0, 'G':0.5, 'T':0.0},
	'Y': {'A':0.0, 'C':0.5, 'G':0.0, 'T':0.5},
	'M': {'A':0.5, 'C':0.5, 'G':0.0, 'T':0.0},
	'K': {'A':0.0, 'C':0.0, 'G':0.5, 'T':0.5},
	'W': {'A':0.5, 'C':0.0, 'G':0.0, 'T':0.5},
	'S': {'A':0.0, 'C':0.5, 'G':0.5, 'T':0.0},
	'B': {'A':0.01, 'C':0.33, 'G':0.33, 'T':0.33},
	'D': {'A':0.33, 'C':0.01, 'G':0.33, 'T':0.33},
	'H': {'A':0.33, 'C':0.33, 'G':0.01, 'T':0.33},
	'V': {'A':0.33, 'C':0.33, 'G':0.33, 'T':0.01},
	'N': {'A':0.25, 'C':0.25, 'G':0.25, 'T':0.25},
}

def manhattan(d1, d2):
	'''Calculates the Manhattan distance between two nucleotide distribution 
	frequencies each represented as a python list'''
	d = 0
	for nt in NT:
		d += abs(d1[nt] - d2[nt])
	return d

def motif2seq(m):
	'''Converts a position weight matrix (represented as an array of 
	dictionaries) into concensus sequence '''
	seq = []
	for i in range(len(m)):
		min = 1e9
		best = None
		for s in IUPAC:
			d = manhattan(m[i], IUPAC[s])
			if d < min:
				min = d
				best = s
		seq.append(best)
	return ''.join(seq)
	
def anti_motif(m):
	'''Returns the reverse complement of a motif formatted as a position weight
	 matrix'''
	a = []
	for i in range(len(m)):
		d = {}
		d['A'] = m[i]['T']
		d['C'] = m[i]['G']
		d['G'] = m[i]['C']
		d['T'] = m[i]['A']
		a.append(d)
	reversed(a)
	return a

def motif_similarity(m1, m2):
	'''Finds and scores all possible overlaps of two motifs (position weight
	matrix format) using an adapted Manhattan distance score and returns the
	highest score and its percentage of the maximum possible score  '''
	m = []
	maxpos = None
	if len(m1) <= len(m2):
		maxpos = 2*len(m1)
	else:
		maxpos = 2*len(m2)
	
	for i in range(len(m1)):
		m.append([])
		for j in range(len(m2)): m[i].append([])
	
	# init first row and column
	for i in range(len(m1)): m[i][0] = 2 - manhattan(m1[i], m2[0])
	for j in range(len(m2)): m[0][j] = 2 - manhattan(m1[0], m2[j])

	# compute diagonals and find best score
	max_score = 0	
	for i in range(1, len(m1)):
		for j in range(1, len(m2)):
			m[i][j] = m[i-1][j-1] + 2 - manhattan(m1[i], m2[j])
			if m[i][j] > max_score:
				max_score = m[i][j]
	score_p = max_score / maxpos
	return max_score, score_p
	
def motif_score(m1,m2):
	'''Accepts two position weight matrices (motifs) and determines if they
	are reverse complements before scoring and returning the highest score
	between the two '''
	s1 = motif_similarity(m1, m2)
	s2 = motif_similarity(m1, anti_motif(m2))
	if s1 >= s2:
		return s1
	else:
		return s2

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



	
def score_motifbit(motif):
	'''Finds the informational content of the motif by summing up the 
	information content of each individual position  '''
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


def run_meme(memepath,promoterfile,m,o,nummotifs,maxw,minw):
	'''Will compile selected meme parameters into a command to run downloaded
	MEME software, and then will extract important information from meme output
	file. '''
	meme = f'{memepath} {promoterfile} -dna -markov_order {o} -mod {m}\
	-nmotifs {nummotifs} -maxw {maxw} -minw {minw} -revcomp 2>/dev/null'
	os.system(meme)
	meme_info = read_memetxt('meme_out/meme.txt')
	motifs, motif_info = memepwm('meme_out/meme.txt')
	return motifs, meme_info, motif_info
	



def performance_bg(motif,motifs):
	'''Finds global similarity score between the given motif and the motif
	found by meme '''
	scores = []
	for i in range(len(motifs)):
		memepwm = motifs[i]
		score = motif_score(motif,memepwm)
		scores.append(score)
	return scores
	


