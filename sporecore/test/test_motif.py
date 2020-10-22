import unittest
import sporecore.motif as sm
import math
import os
import shutil

#NT = ['A', 'C', 'G', 'T']

class TestMotif(unittest.TestCase):
	
	def test_checkmeme(self):
		assert(shutil.which('meme')!= None)
		
	def test_manhattan(self):
		d1 = [{'A':0.5, 'C':0.0, 'G':0.5, 'T':0.0}]
		d2 = [{'A':0.0, 'C':0.5, 'G':0.0, 'T':0.5}]
		score = sm.motif_score(d1,d2)
		assert(score[0]==0)
		
	def test_antimotif(self):
		motif = sm.read_JASPAR("data/testjaspar.1.jaspar")
		anti = sm.anti_motif(motif)
		assert(anti[0]['G']==0.1)
		
	def test_motifsimilarity(self):
		m1 = sm.read_JASPAR("data/testjaspar.1.jaspar")
		m2 = sm.read_JASPAR("data/testjaspar.1.jaspar")
		sim = sm.motif_similarity(m1, m2)
		assert(sim[0]==18.0)
		
	def test_motifscore(self):
		m1 = sm.read_JASPAR("data/testjaspar.1.jaspar")
		m2 = sm.read_JASPAR("data/testjaspar2.1.jaspar")
		score = sm.motif_score(m1,m2)
		assert(score[0] == 16.0)	
	
	def test_readmemetxt(self):
		sites = sm.read_memetxt("data/meme.txt")
		assert(sites[0][4]==99)
		assert(len(sites)==63)
	
	def test_memepwm(self):
		motifs, stats = sm.memepwm("data/meme.txt")
		assert(math.isclose(motifs[0][0]["C"], 0.304348))
		assert(stats[0][1]==12)
	
	def test_readtestmotif(self):
		jpositions, numjsites = sm.read_testmotif("data/fasta_output.fa")
		assert(jpositions[1][1][0]== '64')
		assert(numjsites==10)
		
	def test_readjaspar(self):
		motif = sm.read_JASPAR("data/testjaspar.1.jaspar")
		assert(motif[0]["G"]== 0.9)
		
	def test_motif2seq(self):
		m = sm.read_JASPAR("data/testjaspar.1.jaspar")
		seq = sm.motif2seq(m)
		assert(seq == 'GAGGTMAAT')
		
		
	#is this legal?
	def test_generatesite(self):
		motif = sm.read_JASPAR("data/testjaspar.1.jaspar")
		motifseq = sm.generate_site(motif)
		assert(len(motifseq)==9)
		
	'''
	def test_newpwm(self):
		sites = sm.read_memetxt("data/meme.txt")
		print(sites)
		newmot = sm.new_pwm(sites)
		print(newmot)
	'''
	
	def test_posaccuracy(self):
		sites = sm.read_memetxt("data/meme.txt")
		jpositions, numjsites = sm.read_testmotif("data/fasta_output.fa")
		motifs, stats = sm.memepwm("data/meme.txt")
		jpwm = sm.read_JASPAR("data/testjaspar.1.jaspar")
		memepwm = motifs[1]
		jpos = int(jpositions[1][1][0])
		mpos = sites[0][4]
		fp, posdis, fl, overlap, overlap_p = sm.pos_accuracy(mpos,jpos,\
		len(memepwm),len(jpwm))
		assert(fp == 0)
		assert(posdis == 98)
		assert(fl == 1)
		assert(overlap == 0)
		assert(overlap_p == 0)
		
	'''
	def test_globalmotcompare(self):
		#bg = {'A':0.25,'C':0.25,'G':0.25,'T':0.25}
		motif1 = sm.read_JASPAR("data/testjaspar.1.jaspar")
		motifs, stats = sm.memepwm("data/meme.txt")
		motif2 = motifs[1]
		score, pscore = sm.global_motcompare(motif1,motif2)
		print(score)
		print(pscore)
		#assert(score == 9.233334000000001)
		#assert(pscore == 0.38472225000000004)
	'''
	
	'''
	def test_localmotcompare(self):
		motif1 = sm.read_JASPAR("data/testjaspar.1.jaspar")
		motifs, stats = sm.memepwm("data/meme.txt")
		motif2 = motifs[1]
		score, scorep = sm.local_motcompare(motif1,motif2)
		assert(score == 6.800002)
		assert(scorep == 0.28333341666666667)
	'''
		
	def test_scoremotifbit(self):
		motif = sm.read_JASPAR("data/testjaspar.1.jaspar")
		score, p = sm.score_motifbit(motif)
		assert(score == 12.152431191966443)
		assert(p == 0.6751350662203579)

	"""
	def test_calcbg(self):
		promfile = "data/fasta_output.fa"
		bg = sm.calcbg_frompromfile(promfile)
		assert(bg["G"] == 0.2098901098901099)
	"""
	
	def test_runmeme(self):
		promoterfile = "data/fasta_output.fa"
		memepath = 'meme'
		o = 0
		m = 'zoops'
		nummotifs = 2
		cmd = f'{memepath} {promoterfile} -dna -markov_order {o} -mod {m}\
		-nmotifs {nummotifs} -revcomp'
		os.system(cmd)
		meme_info = sm.read_memetxt('meme_out/meme.txt')
		motifs, motif_info = sm.memepwm('meme_out/meme.txt')
		assert(meme_info[2][3]=='seq-5')
		assert(motifs[0][1]["A"]==0.9)
		assert(motif_info[1][1]==11)
		
	def test_performancebg(self):
		promfile = "data/fasta_output.fa"
		#bg = sm.calcbg_frompromfile(promfile)
		motif = sm.read_JASPAR("data/testjaspar.1.jaspar")
		motifs, stats = sm.memepwm("data/meme.txt")
		perf = sm.performance_bg(motif,motifs)
		#assert(perf[1][1]==0.38119660897435903)
		
	'''
	def test_performancemo(self):
		promfile = "data/fasta_output.fa"
		motif = sm.read_JASPAR("data/testjaspar.1.jaspar")
		motifs, stats = sm.memepwm("data/meme.txt")
		perf = sm.performance_mo(motif,motifs)
		assert(perf[1][0]==6.800002)
	'''		
		
if __name__ == '__main__':
	unittest.main(verbosity=2)
 
