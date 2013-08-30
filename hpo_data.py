import cPickle as pickle
from node import *

f = open('hpo_data', 'rb')
data = pickle.load(f)
code_to_term = data['code_to_term']
hpo_root = data['hpo_root']
terms = data['terms']
f.close()
