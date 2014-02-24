from node import *
import pickle as pickle
import sys
import re
import queue

f = open(sys.argv[1], 'r', encoding='utf-8')

code_to_term = {}
line = ('','','')
curr_term = None

for l in f:
	line = str.partition(l.rstrip(),': ')
	if line[0] == '[Term]':
		curr_term = Term()
	elif line[0] == 'id':
		curr_term.code = line[2]
		code_to_term[curr_term.code] = curr_term
	elif line[0] == 'name':
		curr_term.name = line[2]
	elif line[0] == 'is_a':
		curr_term.parents.append(str.partition(line[2]," ")[0])
	elif line[0] == 'alt_id':
		code_to_term[line[2]] = curr_term
		
f.close()

for code in code_to_term:
	curr_term = code_to_term[code]
	if len(curr_term.parents)>0 and isinstance(curr_term.parents[0],str):
		for i in range(len(curr_term.parents)):
			curr_term.parents[i] = code_to_term[curr_term.parents[i]]
			curr_term.parents[i].children.append(curr_term)

# Changed by OJB to drop onset and inheritance terms
hpo_root = code_to_term['HP:0000001']
#hpo_root = code_to_term['HP:0000118']
#hpo_root.parents = []

hpo_root.map(Node.set_proper_children)
for node in hpo_root:
	node.clean()

#terms = sorted([t for t in set(code_to_term.values()) if t.adam() is hpo_root])
terms = []
for t in set(code_to_term.values()):
	t.count = 0

q = queue.Queue()
q.put(hpo_root)

while not q.empty():
	t = q.get()
	for c in t.children:
		c.count += 1
		if c.count == len(c.parents):
			q.put(c)
	terms.append(t)


for t in set(code_to_term.values()):
        del t.count

#to_drop = set()
#for code, t in code_to_term.items():
#	if t.count == 0 and t is not hpo_root:
#		to_drop.add(code)

#print('Dropping {} terms not under root'.format(len(to_drop)), file=sys.stderr)
#for code in to_drop:
#	del code_to_term[code]

#for t in set(code_to_term.values()):
#	del t.count

#for t in terms:
#	print(t.name)

print('Found {} terms'.format(len(terms)))

sys.setrecursionlimit(10000)
data = {'code_to_term': code_to_term, 'hpo_root': hpo_root, 'terms': terms}
out = open('hpo_data', 'wb')
pickle.dump(data, out)
out.close()
