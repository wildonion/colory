

import matplotlib.pyplot as plt
import plotly.graph_objects as go
import numpy as np
import random
plt.style.use('ggplot')


__all__ = ['genetic_process', 'population'] # export only these classes



class gene:
	def __init__(self, allele):
		self.allele = allele # value of gene


class chromosome:
	def __init__(self, genes):
		self.genes_objects = np.array([gene(g) for g in genes]) # list of all genes objects
		self.genes = genes # list of genes value

	def fitness(self):
		pass

	def __getitem__(self, locus):
		return self.genes_objects[locus]

	def __len__(self):
		return len(self.genes)


class population:
	def __init__(self, amount=200, colors=3, chromosomes=None):
		self.amount = amount
		self.colors = colors
		self.pop = [] # list of all chromosomes (solutions)
		if not chromosomes: self.__init_pop()
		else: self.pop = [chromosome(c) for c in chromosomes]
	
	def __init_pop(self):
		for i in range(self.amount):
			c = None # build each chromosome with at least self.color genes
			self.pop.append(chromosome(c))

	def fitness_score(self):
		scores = [] # scores
		for chromosome in self.pop:
			scores.append(chromosome.fitness())
		scores, population = np.array(scores), np.array([c.genes for c in self.pop]) # list of all chromosomes' scores, population of all chromosomes with their genes
		indices = np.argsort(scores) # return the indices of sorted scores in ascending order - used in rank selection
		descending_scores = scores[indices][::-1] # sorted scores in descending order
		descending_population_of_scores = population[indices, :][::-1] # sorted population of chromosomes scores in descending order
		return list(descending_scores), list(descending_population_of_scores) # return descending order of population of none object genes chromosome and scores

		
	def __len__(self):
		return len(self.pop)

	def __getitem__(self, idx):
		return self.pop[idx]


class genetic_process:
	def __init__(self, generations, population, parents, selection_method, crossover_method, mutation_method, mutation_rate, crossover_rate):
		self.generations = generations
		self.population = population
		self.parents = parents
		self.mutation_rate = mutation_rate
		self.crossover_rate = crossover_rate
		self.selection_method = selection_method
		self.crossover_method = crossover_method
		self.mutation_method = mutation_method
		self.population_after_fitness = []
		self.parents_population = []
		self.population_after_crossover = []
		self.best_chromosomes = []
		self.best_scores = []


	def run(self):
		for i in range(self.generations):
			print(f"🧬 Generation --- {i+1}")
			scores, self.population_after_fitness = self.population.fitness_score(self.model, self.data)
			print(f"\t▶  Best Score for Two Chromosomes --- {scores[:2]}\n")
			# =================== GA Operators ===================
			self.__selection() # select best fitness as parents
			self.__crossover() # parents mating pool
			self.__mutation() # mutating genes
			# ====================================================
			self.best_chromosomes.append(self.population_after_fitness[0])
			self.best_scores.append(scores[0])


	def __crossover(self):
		offspring = self.parents_population
		if self.crossover_method == "single_point":
			raise NotImplementedError # TODO
		elif self.crossover_method == "two_point":
			raise NotImplementedError # TODO
		elif self.crossover_method == "multi_point":
			raise NotImplementedError
		else:
			raise NotImplementedError


	def __mutation(self):
		offspring_after_mutation = []
		if self.mutation_method == "flipping":
			raise NotImplementedError
		elif self.mutation_method == "reversing":
			raise NotImplementedError # TODO
		elif self.mutation_rate == "interchanging":
			raise NotImplementedError # TODO
		else:
			raise NotImplementedError

	def __selection(self):
		population_next_generation = []
		if self.selection_method == "roulette_wheel":
			fitness_population = sum(self.population.fitness_score(self.model, self.data)[0]) # sum of all scores (fitnesses)
			individual_expected_values = [c.fitness(self.model, self.data)/fitness_population for c in self.population] # all chromosomes prob (exprected values)
			cum_prob = [sum(individual_expected_values[:i+1]) for i in range(len(individual_expected_values))] # cumulative sum of chromosomes exprected values (prob)
			for i in range(self.parents):
				r = random.random()
				for j, chromosome in enumerate(self.population):
					if cum_prob[j] >= r:
						population_next_generation.append(self.population[j].genes)
			self.parents_population = population_next_generation
		elif self.selection_method == "rank":
			raise NotImplementedError
		elif self.selection_method == "tournament":
			raise NotImplementedError # TODO
		else:
			raise NotImplementedError

	def plot(self, lib="plotly"):
		if lib == "matplotlib":			
			plt.plot(self.best_scores)
			plt.xlabel("Generation")
			plt.ylabel("Best Fitness")
			plt.savefig("fitness_generation.png")
		elif lib == "plotly":
			fig = go.Figure()
			fig.add_trace(go.Scatter(self.best_scores, mode='lines', name='Fitness Generation'))
			fig.update_xaxes(title_text='Generation')
			fig.update_yaxes(title_text='Best Fitness')
			fig.show()

	def save(self):
		print(f"▶ Saving Best chromosomes and best scores")
		np.save(f"utils/best_chromo_in_{self.generation}_generations.npy", self.best_chromosomes)
		np.save(f"utils/best_scores_in_{self.generation}_generations.npy", self.best_scores)