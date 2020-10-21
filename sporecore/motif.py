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
	d = 0
	for nt in NT:
		d += abs(d1[nt] - d2[nt])
	return d

def motif2seq(m):
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
	m = []
	maxpos = None
	#print(len(m1))
	#print(len(m2))
	if len(m1) <= len(m2):
		#print('true')
		maxpos = 2*len(m1)
	else:
		#print('false')
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


"""
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
"""
	
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
	
#memetools

'''
def convert_argtovar(minprom, maxprom, promstep, minseq,maxseq,seqstep):
	Accepts command line parameter arguments and converts them to ranges of
	parameters to be tested in MEME
	promoter = []
	num_seq = []
	for i in range(minprom,maxprom+1,promstep):
		if i != 0:
			promoter.append(i)
	for i in range(minseq,maxseq+1,seqstep):
		if i != 0:
			num_seq.append(i)
	return promoter, num_seq
	
	

#fix generate promter to create files for condenseddata
def generate_promoter_mo(dnafile, jasparfile, p, n, o, freq,i): 
	Creates a fasta file of promoter sequences with embedded binding motifs.
	Background sequence is generated using markov models of DNA sequence from a
	fasta file input
	tmpfile = f'/tmp/testmotif{os.getpid()}_{p}_{n}_{i}.fa' 
	#f'python3 markov_Mbed.py --jasparfile {jasparfile} \
	#--dnafile {dnafile} --markov_order {o} --numseq {n} --seqlen {p} --freq \
	#{freq} --bothstrands > {tmpfile} '
	cmd = f'python3 markov_Mbed.py --jasparfile {jasparfile} \
	--dnafile {dnafile} --markov_order {o} --numseq {n} --seqlen {p} --freq \
	{freq} --bothstrands > {tmpfile} '
	os.system(cmd)
	return tmpfile

def generate_promoter_bg(jasparfile, p, n, freq, background,i,condenseddata):
	Creates a fasta file of promoter sequences with embedded binding motifs.
	Background sequence is generated using a nucleotide distribution frequency
	input from the command line.
	tmpfile = f'/tmp/testmotif{os.getpid()}_{p}_{n}_{i}.fa' 
	cmd = f'python3 Mbed.py --jasparfile {jasparfile} \
	--numseq {n} --seqlen {p} --freq {freq} --PA {background["A"]} --PC \
	{background["C"]} --PG {background["G"]} --PT {background["T"]} \
	--bothstrands > {tmpfile} '
	os.system(cmd)
	return tmpfile
'''
	

def calcbg_frompromfile(promoterfile):
	'''Calculates a nucleotide distribution frequency for the background 
	DNA sequence from the promoters in the generated FASTA file'''
	A = 0
	C = 0
	G = 0
	T = 0
	with open(promoterfile) as pf:
		for line in pf.readlines():
			if not line.startswith('>'): 
				line = line.strip()
				#print('line', line)
				for i in range(len(line)):
					if line[i] == 'a': A += 1 
					if line[i] == 'c': C += 1 
					if line[i] == 'g': G += 1 
					if line[i] == 't': T += 1 
	total = A+C+G+T
	if total > 0:
		A = A/total
		C = C/total
		G = G/total
		T = T/total 
		background = {'A':A, 'C':C, 'G': G, 'T':T}
		return background
	else:
		raise ValueError('Promoter file empty, check file parameters')	
	
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
	
#can combine these now	
#performance hands back one number 
def performance_bg(motif,motifs,background):
	'''Finds global similarity score between the given motif and the motif
	found by meme '''
	scores = []
	for i in range(len(motifs)):
		memepwm = motifs[i]
		score = motif_score(motif,memepwm)
		scores.append(score)
	return scores
	
def performance_mo(motif,motifs):
	'''Finds local similarity score between the given motif and the motif found
	 by meme.'''
	scores = []
	for i in range(len(motifs)):
		memepwm = motifs[i]
		score = local_motcompare(motif,memepwm)
		scores.append(score)
	return scores

'''
def get_memedata(promoter_file, meme_info, j_info, distance_scores, motif_info\
,jpwm,p,n,m,o,iteration):
	#Assembles data extracted from the meme output into a single array that
	 #can be displayed or further analyzed
	result = []
	fp = 0
	fn = []
	for i in range(len(meme_info)):
		seq = meme_info[i][3]
		motif = meme_info[i][0]
		nsites = meme_info[i][1]
		m_strand = meme_info[i][2]
		m_pos = meme_info[i][4]
		p_val = meme_info[i][5]
		for j in range(len(j_info)):
			if seq == j_info[j][0]:
				j_pos = j_info[j][1][0]
				if j_pos != '':
					j_strand = j_info[j][1][1]
				else:
					j_strand = ''
				break
		for j in range(len(motif_info)):
			if motif == motif_info[j][0]:
				score = distance_scores[j][0] 
				p_score = distance_scores[j][1]
				meme_wid = motif_info[j][1]
				meme_eval = motif_info[j][3]
				break	
		fpos,posdis,fl,overlap,overlap_p = pos_accuracy(m_pos, j_pos,\
		meme_wid,len(jpwm))
		result.append((promoter_file, seq, motif,nsites,m_pos, m_strand,\
		meme_wid,j_pos,j_strand,len(jpwm),p_val,score,p_score,meme_eval,fpos\
		,0,posdis,fl,overlap,overlap_p,p,n,m,o,iteration))
		fn.append(seq)
	return result,fn



#need to fix false negs everything 
def find_falsenegs(fn,j_info,p,n,m,o,promoter_file,jpwm,iteration):
	Determines which promoters with embedded motifs from the fasta files
	were not detected by meme. 
	#print('j_info',j_info)
	#print('fn',fn)
	false_negs = []
	for i in range(len(j_info)):
		#print('j_info[i][0]',j_info[i][0])
		#print(fn[i])
		if j_info[i][0] not in fn:
			if j_info[i][1][0] != '':
				#print('is it ""?',j_info[i][1][0])
				false_negs.append(( promoter_file,j_info[i][0],'','','','','',\
				j_info[i][1][0],j_info[i][1][1],len(jpwm),'','','','',0,1,99,1\
				,'','',p,n,m,o,iteration))
	return false_negs	



def find_condensedstats(result,fn,motif_info,j_info,promoter_file,bits,p_bits,\
p,n,m,o,false_neg,numjsites,iteration):
	Creates condensed meme results focused on the performance of each motif
	 found in meme 
	condensed_stats = []
	for i in range(len(motif_info)):
		fl_ct = 0
		fp_ct = 0
		for j in range(len(result)):
			if result[j][2] == motif_info[i][0]:
				nsites = int(result[j][3])
				fl_ct += (result[j][17])
				fp_ct += (result[j][14])
				score = (result[j][11])
				score_p = (result[j][12])
				evalue = (result[j][13])
		condensed_stats.append((promoter_file,motif_info[i][0],nsites,bits,p_bits,\
		evalue,score,score_p,fl_ct/nsites,(nsites-fl_ct)/nsites,\
		fp_ct/nsites,len(false_neg)/numjsites,p,n,m,o,iteration))
	return condensed_stats	


#make sure this is good, make uncondensed have file name
def present_info(promoter_file,results, bits,p_bits,scores,fn,\
	j_info,motif_info,p,n,m,o,condenseddata):
	Determines what presentation of data should look like based on the input 
	parameters 
	present = []
	false_neg = find_falsenegs(fn, j_info,p,n,m,o,promoter_file)
	if condenseddata:
		present = (find_condensedstats(results,fn,motif_info,j_info,\
		promoter_file,bits,p_bits,p,n,m,o))		
	else:
		 for i in range(len(results)):
		 	present.append((results[i]))
		 for i in range(len(false_neg)):
		 	present.append(false_neg[i])
	return present
	


def avg_condensedstats(final,motif_info,p,n,m,o,r):
	Averages performance of condensed motif information for a more accurate
	 understanding of motif performance 
	output = []
	if len(final) > 1:
		for i in range(len(motif_info)):
			mot = motif_info[i][0]
			avg_nsites = []
			avg_score = []
			avg_eval = []
			avg_score = []
			avg_pscore = []
			avg_frate = []
			avg_srate = []
			avg_fprate = []
			avg_fnrate = []
			for j in range(len(final)):
				for k in range(len(final[j])):
					if  mot == final[j][k][1]:
						avg_nsites.append(final[j][k][2])
						avg_eval.append(final[j][k][5])
						avg_score.append(final[j][k][6])
						avg_pscore.append(final[j][k][7])
						avg_frate.append(final[j][k][8])
						avg_srate.append(final[j][k][9])
						avg_fprate.append(final[j][k][10])
						avg_fnrate.append(final[j][k][11])
			output.append((final[j][k][0],mot,statistics.mean(avg_nsites),\
			final[j][k][2],final[j][k][3], statistics.mean(avg_eval),\
			statistics.mean(avg_score), statistics.mean(avg_pscore),\
			statistics.mean(avg_frate), statistics.mean(avg_srate),\
			statistics.mean(avg_fprate),statistics.mean(avg_fnrate),p,n,m,o,\
			r+1,statistics.stdev(avg_eval),statistics.stdev(avg_score)\
			,statistics.stdev(avg_pscore),statistics.stdev(avg_frate),\
			statistics.stdev(avg_srate),statistics.stdev(avg_fprate)\
			,statistics.stdev(avg_fnrate)))
		return output
	else:
		for i in range(len(final)):
			for j in range(len(final[i])):
				output.append(final[i][j])
		return output
'''