


import networkx as nx
import plotly.graph_objects as go
import numpy as np
import random


__all__ = ['genetic_process', 'population'] # export only these classes



class gene:
	def __init__(self, allele):
		self.allele = allele # value of gene


class chromosome:
	def __init__(self, genes):
		self.gene_objects = np.array([gene(g) for g in genes]) # list of all gene objects
		self.genes = np.array(genes) # list of all gene values

	def fitness(self, adj_mat, colors):
		alpha, beta              = 0.8, 0.2
		total_edges              = int(adj_mat[adj_mat==1].sum()/2) # sum up thoes place in adjacency matrix that their value is 1, because there is an edge between that row and that column of the matrix
		invalid_genes            = self.__invalid_genes_objective(adj_mat) # total number of invalid genes, it might be equal to zero which means that there is no invalid coloring
		minimum_number_of_colors = self.__minimum_genes_objective(colors) # this is the m′, the minimum number of colors used in this chromosome, it might be equal to the number of total input colors 
		m_prime_normalized       = minimum_number_of_colors - 1 / colors.shape[0] # normalizing the range of m′
		invalid_normalized       = invalid_genes / total_edges # normalizing the range of invalid genes
		total_fitness            = (alpha*invalid_normalized) + (beta*m_prime_normalized) # total fitness of this genotype (chromosome) 
		return total_fitness


	def __invalid_genes_objective(self, adj_mat):
		invalid = 0
		for i in range(len(self.gene_objects)):
			if self.gene_objects[i].allele == self.gene_objects[(i+1)%len(self.gene_objects)].allele:
				if adj_mat[i][(i+1)%len(self.gene_objects)] == 1: # there is an edge between this node and the next node
					invalid += 1
		return invalid


	def __minimum_genes_objective(self, colors):
		_, unique_indices = np.unique(self.genes, return_index=True)
		unique_colors = self.genes[unique_indices]
		m_prime = unique_colors.shape[0] if unique_colors.shape[0] < colors.shape[0] else colors.shape[0]
		return m_prime # minimum value is 1 and maximum value is colors.shape[0]

	def __getitem__(self, locus):
		return self.gene_objects[locus]

	def __len__(self):
		return self.genes.shape[0]


class population:
	def __init__(self, amount=200, colors=3, adj_mat=None, chromosomes=None):
		self.amount  = amount
		self.colors  = colors
		self.adj_mat = adj_mat
		if not chromosomes: self.__init_pop()
		else: self.pop = np.array([chromosome(c) for c in chromosomes])
	
	def __init_pop(self): # permutation encoding
		self.pop = [chromosome([np.where(self.colors == np.random.choice(self.colors, 1))[0] for _ in range(self.adj_mat.shape[0])]) for _ in range(self.amount)] # list of all chromosomes (solutions) - build each chromosome genes with a random inex of colors list

	def fitness_scores(self):
		fitness_scores = [chromosome.fitness(self.adj_mat, self.colors) for chromosome in self.pop] # all chromosomes fitness inside the generated population

	def __len__(self):
		return self.pop.shape[0]

	def __getitem__(self, idx):
		return self.pop[idx]


class genetic_process:
	def __init__(self, generations, population, parents, selection_method, crossover_method, mutation_method, replacement_method, mutation_rate, crossover_rate):
		self.generations = generations
		self.population = population
		self.parents = parents
		self.mutation_rate = mutation_rate
		self.crossover_rate = crossover_rate
		self.selection_method = selection_method
		self.crossover_method = crossover_method
		self.mutation_method = mutation_method
		self.replacement_method = replacement_method
		self.population_after_fitness = []
		self.population_after_selection = []
		self.population_after_crossover = []
		self.population_after_mutation = []
		self.best_chromosomes = []
		self.best_fitness_scores = []


	def run(self):
		for i in range(self.generations): # TODO - stopping criterion: if fitness == 0
			print(f"🧬 Generation --- {i+1}")
			scores, self.population_after_fitness = self.population.fitness_scores(self.model, self.data)
			print(f"\t▶  Best Score --- {scores[:2]}\n")
			# =================== GA Operators ===================
			self.__selection() # select best chromosomes as parents
			self.__crossover() # parents mating pool
			self.__mutation() # mutating genes
			self.__replacement() # replacing old population
			# ====================================================
			self.best_chromosomes.append(self.population_after_fitness[0])
			self.best_fitness_scores.append(scores[0])

	def __selection(self):
		population_after_selection = []
		if self.selection_method == "roulette_wheel":
			fitness_population = sum(self.population.fitness_scores(self.model, self.data)[0]) # sum of all scores (fitnesses)
			individual_expected_values = [c.fitness(self.model, self.data)/fitness_population for c in self.population] # all chromosomes prob (exprected values)
			cum_prob = [sum(individual_expected_values[:i+1]) for i in range(len(individual_expected_values))] # cumulative sum of chromosomes exprected values (prob)
			for i in range(self.parents):
				r = random.random()
				for j, chromosome in enumerate(self.population):
					if cum_prob[j] >= r:
						population_after_selection.append(self.population[j].genes)
			self.population_after_selection = population_after_selection # parents population
		elif self.selection_method == "rank":
			self.population_after_selection = population_after_selection # parents population
			raise NotImplementedError # TODO
		elif self.selection_method == "tournament":
			self.population_after_selection = population_after_selection # parents population
			raise NotImplementedError # TODO
		else:
			raise NotImplementedError

	def __crossover(self):
		# http://ijcsit.com/docs/Volume%205/vol5issue06/ijcsit2014050673.pdf
		offspring = []
		if "point" in self.crossover_method:
			point = self.crossover_method.split("_")[0] # point to crossover
			self.population_after_crossover = offspring
			raise NotImplementedError # TODO
		if self.crossover_method == "uniform":
			self.population_after_crossover = offspring
			raise NotImplementedError # TODO
		elif self.crossover_method == "pmx":
			self.population_after_crossover = offspring
			raise NotImplementedError # TODO
		elif self.crossover_method == "ox":
			self.population_after_crossover = offspring
			raise NotImplementedError # TODO
		else:
			raise NotImplementedError


	def __mutation(self):
		# http://ijcsit.com/docs/Volume%205/vol5issue03/ijcsit20140503404.pdf
		offspring_after_mutation = []
		if self.mutation_method == "swap":
			self.population_after_mutation = offspring_after_mutation
			raise NotImplementedError # TODO
		elif self.mutation_method == "creep":
			self.population_after_mutation = offspring_after_mutation
			raise NotImplementedError # TODO
		elif self.mutation_rate == "interchanging":
			self.population_after_mutation = offspring_after_mutation
			raise NotImplementedError # TODO
		elif self.mutation_rate == "reversing":
			self.population_after_mutation = offspring_after_mutation
			raise NotImplementedError # TODO
		else:
			raise NotImplementedError

	def __replacement(self):
		population_next_generation = []
		if self.mutation_method == "generational_elitism":
			self.population = population_next_generation
			raise NotImplementedError # TODO
		elif self.mutation_method == "generational_gap":
			self.population = population_next_generation
			raise NotImplementedError # TODO
		elif self.mutation_rate == "steady_state":
			self.population = population_next_generation
			raise NotImplementedError # TODO
		else:
			raise NotImplementedError

	def plot(self):
		fig = go.Figure()
		fig.add_trace(go.Scatter(self.best_fitness_scores, mode='lines', name='Fitness Generation'))
		fig.update_xaxes(title_text='Generation')
		fig.update_yaxes(title_text='Best Fitness')
		fig.show()


	def draw(self):
		# https://plotly.com/python/v3/3d-network-graph/
		# https://plotly.com/python/network-graphs/
		pass

	def save(self):
		print(f"▶ Saving Best chromosomes and best scores")
		np.save(f"utils/best_chromo_in_{self.generation}_generations.npy", self.best_chromosomes)
		np.save(f"utils/best_fitness_scores_in_{self.generation}_generations.npy", self.best_fitness_scores)