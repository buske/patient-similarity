from all_data import *
import numpy as np
import random


#Generates some random patients
#num_patients is a list containing the number of patients to generate for each disease
#min_dis_traits and max_dis_traits specify the minimum and maximum number of annotations for candidate diseases
#The differences between the similarity algorithms show up better with well annotated patients

num_diseases = 300
num_patients = [501] + [10]*(num_diseases-1)
min_dis_traits = 40
max_dis_traits = 80

randwalk_patients = []
diseases = random.sample([d for d in diseases if min_dis_traits <= len(d.traits) <= max_dis_traits], num_diseases)
for d,n in zip(diseases,num_patients):
	for i in xrange(n):
		randwalk_patients.append(d.gen_rand(steps = 3))

print 'Average trait number:', np.mean([len(p.traits) for p in randwalk_patients])

f = open('test.txt', 'w')
for p in randwalk_patients:
	f.write(p.print_traits() + '\n')
f.close()
