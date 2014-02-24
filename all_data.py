import pickle as pickle
from node import *
import phenotype
from math import *
import numpy as np
import queue

f = open('all_data', 'rb')
data = pickle.load(f)
code_to_term = data['code_to_term']
hpo_root = data['hpo_root']
name_to_dis = data['name_to_dis']
terms = data['terms']
diseases = data['diseases']

f.close()
