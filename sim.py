from all_data import *
from math import *

def sf(f):
	def fp(s):
		return set.union(*[f(x) for x in s])
	return fp

#returns the union of all the descendants/ancestors of a set
sf_get_descendants = sf(Node.get_descendants)
sf_get_ancestors = sf(Node.get_ancestors)

def cost(traits):
	return sum([t.inf for t in traits])

def sim(p1,p2,x):
	"""Returns the similarity between two sets of traits, p1 and p2
	If x is false, it uses the bound P(t|t.parents) > exp(max_{p in t.parents} t.inf - p.inf)
	If x is true, it uses the bound P(t|t.parents) > t.freq / (sum_{c a common child of t.parents} t.freq)"""

	#intersection is the set of all common ancestors of p1 and p2
	intersection = sf_get_ancestors(p1.traits) & sf_get_ancestors(p2.traits)
	for t in intersection:
		t.p1count = 0
		t.p2count = 0

	print('Patient 1 terms and IC')
	for t in p1.traits:
		print('  {} {:.6f}'.format(t.code, t.inf))
	print('Patient 2 terms and IC')
	for t in p2.traits:
		print('  {} {:.6f}'.format(t.code, t.inf))


	#for each t in intersection, t.p1count is the number of traits in p1 which are a descendant of t
	for t in p1.traits:
		for y in (t.get_ancestors() & intersection):
			y.p1count += 1
	for t in p2.traits:
		for y in (t.get_ancestors() & intersection):
			y.p2count += 1
	print('Found {} common ancestors'.format(len(intersection)))
	scores = []
	for t in intersection:
		t.mincount = min(t.p1count, t.p2count)
		scores.append((t.code, t.mincount, t.cond_inf2))
	scores.sort()
	for L in scores:
		print('{} {} {:.6f}'.format(*L))

	#costp1 is sum of the information contents of all of p1's traits
	costp1 = cost(p1.traits)
	costp2 = cost(p2.traits)
	print('Patient 1 cost: {:.6f}'.format(costp1))
	print('Patient 2 cost: {:.6f}'.format(costp2))
	#t.cond_inf and t.cond_inf2 are pre-computed in parse_diseases.py
	if x:
		scores = [t.cond_inf * t.mincount for t in intersection]
	else:
		scores = [t.cond_inf2 * t.mincount for t in intersection]

	shared_info = sum(scores)
	print('Shared cost: {:.6f}'.format(shared_info))
	names = [t.name for t in intersection]
	return 2*shared_info/(costp1 + costp2), zip(scores, names)
