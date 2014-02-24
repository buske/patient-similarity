from node import *
from phenotype import *
from hpo_data import *
from math import *
import pickle as pickle
import sys
import re

f = open(sys.argv[1], 'r', encoding='utf-8')
line = ['']
name_to_dis = {}
freq_mods = {'very rare':0.01, 'rare':0.05, 'occasional': 0.075, 'frequent':0.33, 'typical': 0.5, 'variable':0.5, 'common':0.75, 'hallmark':0.9, 'obligate': 1.0}
default_freq = 0.25
eps = 1e-9
percent = re.compile(r"\s*([\d.]*)\s*%\s*")
frac = re.compile(r"\s*(\d*)\s*(of|/)\s*(\d*)\s*")

for l in f:
	line = str.split(l.rstrip(),'\t')
	if len(line) > 2:
		line[2] = line[2].rstrip()
	if line[0] == 'OMIM':
		if line[2] in name_to_dis:
			curr_dis = name_to_dis[line[2]]
		else:
			curr_dis = Disease(line[2])
			name_to_dis[line[2]] = curr_dis
			curr_dis.codes = []
		try:
			trait = code_to_term[line[4]]
		except KeyError:
			pass
		else:
			curr_dis.traits.append(trait)
			curr_dis.codes.append(line[4])
			freq = line[8].lower()
			m1=percent.match(freq)
			m2=frac.match(freq)
			if m1:
				curr_dis.freqs.append(float(m1.group(1))/100)
			elif m2:
				curr_dis.freqs.append(float(m2.group(1))/float(m2.group(3)))
			elif freq in freq_mods:
				curr_dis.freqs.append(freq_mods[freq])
			else:
				if freq:
					print('Couldnt parse: {!r}'.format(freq), file=sys.stderr)
				curr_dis.freqs.append(None)

for dis in name_to_dis.values():
	freqs = {}
	for t in dis.traits:
		freqs[t] = []
	for tup in zip(dis.traits, dis.freqs):
		if tup[1] is not None:
			freqs[tup[0]].append(tup[1])
	dis.traits = []
	dis.freqs = []
	for t in freqs:
		dis.traits.append(t)
		dis.freqs.append(float(sum(freqs[t]))/len(freqs[t]) if len(freqs[t]) > 0 else 0.25)
			
	for i in range(len(dis.traits)):
		dis.traits[i].raw_freq += dis.freqs[i]


assert abs(code_to_term['HP:0000032'].raw_freq - 0.25) < 1e-6
assert abs(code_to_term['HP:0000137'].raw_freq - 0.075) < 1e-6
assert abs(code_to_term['HP:0000177'].raw_freq - 0.05) < 1e-6

term_codes = []
for t in terms:
	t.raw_freq += 0.001
	term_codes.append(t.code)
term_codes.sort()

def sd(f):
	def fp(s):
		return set.intersection(*[f(x) for x in s])
	return fp

#returns the common descendants of a set of phenotypes
sd_get_descendants = sd(Node.get_descendants)

def info(s):
	"""The information of a set of traits is the total probability mass in their common descendants"""
	prob = sum([x.freq for x in sd_get_descendants(s)])
	prob = min(max(prob,eps), 1 - eps)
	return -log(prob)

def compute_probs():
	"""This is where we compute all the probabilities"""
	
	#t.raw_freq is the total number of annotations for a term t weighted by their
	#probabilities in phenotype_annotation.tab, after some cleanup.
	tot = sum(t.raw_freq for t in terms)

	#Here we normalize the frequencies so that sum_t t.freq=1, giving a proper distribution over terms
	for t in terms:
		t.freq = t.raw_freq / tot
		t.freq = min(max(t.freq, eps), 1-eps)
  

	print('Total frequency mass:', tot)
	for code in term_codes[:5]:
		t = code_to_term[code]
		print(code)
		print(t.raw_freq)
		print(t.freq)

	#The t.prob is the total probability mass in the descendants of t
	for t in terms:
		t.prob = sum([x.freq for x in t.get_descendants()])
		t.prob = min(max(t.prob,eps), 1-eps)
		t.inf = -log(t.prob)

	#t.cond_inf is a bound on log P(t|t.parents). It is the information content of t
	#minus the information content of the most informative ancestor.

	#t.cond_inf2 is a sligthly better bound. It is the information content of t
	#minus log sum_{c a common child of t.parents} c.freq
	for t in terms:
		t.cond_inf = t.inf - info(t.get_ancestors().difference(set([t])) | set([hpo_root]))
		t.cond_inf2 = t.inf - max([p.inf for p in t.parents] + [0])

compute_probs()

diseases = [name_to_dis[name] for name in sorted(name_to_dis)]
for i in range(len(diseases)):
	diseases[i].id = i
print("Diseases:", len(name_to_dis))

sys.setrecursionlimit(100000)
data = {'code_to_term': code_to_term, 'hpo_root': hpo_root, 'name_to_dis': name_to_dis, 'terms': terms, 'diseases': diseases}
out = open('all_data', 'wb')
pickle.dump(data, out)
out.close()
