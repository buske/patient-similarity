from node import *
import numpy as np

class Phenotype:
	def __init__(self, traits = []):
		self.traits = traits

	def get_all_traits(self):
		if not hasattr(self, 'all_traits'):
			self.all_traits = set.union(*[t.get_ancestors() for t in self.traits])
		return self.all_traits

class Disease(Phenotype):
	def __init__(self, name = "", traits = []):
		self.name = name
		self.traits = []
		self.freqs = []

	def gen_rand(self, steps = 3):
		new_traits = []
		while len(new_traits) < 5:
			for i in range(len(self.traits)):
				f = self.freqs[i]
				if random.random() < f:
					new_traits.append(self.traits[i].rand_walk(steps))
		return Patient(new_traits, disease = self)

	def gen_prototypical(self):
		raise NotImplementedError()

	def show(self):
		print(self.name)
		for t in sorted([t.name for t in self.traits]):
			print('\t' + t)

class Patient(Phenotype):
	def __init__(self, traits = [], disease = None):
		self.traits = traits
		self.disease = disease
		for t in self.get_all_traits():
			t.patients.add(self)

	def print_compact(self):
		return self.disease.name + '\t' + ''.join([t.code + ' ' for t in self.traits])

	def print_traits(self):
		return ''.join([t.code + ' ' for t in self.traits])
	
	def show(self):
		print(self.disease.name)
		for t in sorted([t.name for t in self.traits]):
			print('\t' + t)
