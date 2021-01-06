


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
		self.genes = np.array(genes) # list of all gene values
		self.gene_objects = np.array([gene(g) for g in self.genes]) # list of all gene objects

	def fitness(self, adj_mat, colors):
		alpha, beta              = 0.8, 0.2
		total_edges              = int(adj_mat[adj_mat==1].sum()/2) # sum up thoes place in adjacency matrix that their value is 1, because there is an edge between that row and that column of the matrix
		invalid_genes            = self.__invalid_genes_objective(adj_mat) # total number of invalid genes, it might be equal to zero which means that there is no invalid coloring
		minimum_number_of_colors = self.__minimum_genes_objective(colors) # this is the m′, the minimum number of colors used in this chromosome, it might be equal to the number of total input colors 
		m_prime_normalized       = (minimum_number_of_colors - 1) / colors.shape[0] # normalizing the range of m′
		invalid_normalized       = invalid_genes / total_edges # normalizing the range of invalid genes
		total_fitness            = (alpha*invalid_normalized) + (beta*m_prime_normalized) # total fitness of this genotype (chromosome) 
		return total_fitness


	def __invalid_genes_objective(self, adj_mat):
		invalid = 0
		for i in range(self.gene_objects.shape[0]):
			next_index = (i+1)%self.gene_objects.shape[0]
			if next_index != 0 and self.gene_objects[i].allele == self.gene_objects[next_index].allele:
				if adj_mat[i][next_index] == 1: # there is an edge between this node and the next node
					invalid += 1
		return invalid


	def __minimum_genes_objective(self, colors):
		_, unique_indices = np.unique(self.genes, return_index=True)
		unique_colors = self.genes[unique_indices] # all unique colors inside a genotype (choromosome)
		m_prime = unique_colors.shape[0] if unique_colors.shape[0] < colors.shape[0] else colors.shape[0]
		return m_prime # minimum value is 1 and maximum value is colors.shape[0]

	def __getitem__(self, locus):
		return self.gene_objects[locus]

	def __len__(self):
		return self.genes.shape[0]


class population:
	def __init__(self, colors=None, amount=None, adj_mat=None, chromosomes=None):
		self.amount  = amount
		self.colors  = colors
		self.adj_mat = adj_mat
		if not chromosomes: self.__init_pop()
		else: self.pop = np.array([chromosome(c) for c in chromosomes])
	
	def __init_pop(self): # permutation encoding
		self.pop = [chromosome([np.where(self.colors == np.random.choice(self.colors, 1))[0][0] for _ in range(self.adj_mat.shape[0])]) for _ in range(self.amount)] # build each chromosome genes with a random index of colors list

	def fitness_scores(self):
		fitness_scores = [chromosome.fitness(self.adj_mat, self.colors) for chromosome in self.pop] # all chromosomes fitness inside the generated population
		fitness_scores, genes_population, chromosome_objects_population = np.array(fitness_scores), np.array([c.genes for c in self.pop]), self.pop
		indices = np.argsort(fitness_scores)
		ascending_fitness_scores = fitness_scores[indices]
		ascending_genes_population_based_on_fitness_scores = genes_population[indices]
		return ascending_fitness_scores, ascending_genes_population_based_on_fitness_scores, chromosome_objects_population

	def __len__(self):
		return self.pop.shape[0]

	def __getitem__(self, idx):
		return self.pop[idx]

	def get_adj_mat(self):
		return self.adj_mat

	def get_colors(self):
		return self.colors


class genetic_process(population):
	def __init__(self, generations, population, parents, selection_method, crossover_method, mutation_method, replacement_method, mutation_rate, crossover_rate):
		self.adj_mat = population.get_adj_mat()
		self.colors = population.get_colors()
		self.generations = generations
		self.population = population
		self.parents = parents
		self.mutation_rate = mutation_rate
		self.crossover_rate = crossover_rate
		self.selection_method = selection_method
		self.crossover_method = crossover_method
		self.mutation_method = mutation_method
		self.replacement_method = replacement_method
		self.genes_population_after_fitness = []
		self.population_after_selection = []
		self.population_after_crossover = []
		self.population_after_mutation = []
		self.best_chromosomes = []
		self.best_fitness_scores = []


	def run(self):
		for i in range(self.generations):
			print(f"🧬 Generation --- {i+1}")
			fitness_scores, self.genes_population_after_fitness, self.chromosome_objects_population = self.population.fitness_scores()
			print(f"\t▶  Best Fitness Scores of Two First Chromosomes --- {fitness_scores[:2]}\n") # minimum fitness scores are the best ones
			# =================== GA Operators ===================
			self.__selection() # select best chromosomes as parents
			self.__crossover() # parents mating pool
			# self.__mutation() # mutating genes
			# self.__replacement() # replacing old population
			# ====================================================
			self.best_chromosomes.append(self.genes_population_after_fitness[0])
			self.best_fitness_scores.append(fitness_scores[0])

	def __selection(self):
		population_after_selection = []
		if self.selection_method == "roulette_wheel": # =====================================================================================
			fitness_population = sum([c.fitness(self.adj_mat, self.colors) for c in self.population]) # sum of all scores (fitnesses)
			individual_expected_values = [(1/c.fitness(self.adj_mat, self.colors))/fitness_population for c in self.population] # all chromosomes prob (exprected values) - we scaled up every chromosome fitness cause big roulette of the wheel belongs to minimum fitnesses 
			cum_prob = [sum(individual_expected_values[:i+1]) for i in range(len(individual_expected_values))] # cumulative sum of chromosomes exprected values (prob)
			for i in range(self.parents):
				r = random.random()
				if cum_prob[i] >= r: # add this chromosome only if its cum_prob is greater than the generated random number
					population_after_selection.append(self.population[i].genes)
			self.population_after_selection = np.array(population_after_selection) # parents population
		
		elif self.selection_method == "rank": # =====================================================================================
			for p in range(self.parents): # the first self.parents chromosomes are the best ones
				population_after_selection.append(self.genes_population_after_fitness[p])
			self.population_after_selection = np.array(population_after_selection) # parents population
		
		elif self.selection_method == "tournament": # =====================================================================================
			k = int(np.log2(len(self.chromosome_objects_population)))
			population_after_tournament = []
			for _ in range(len(self.chromosome_objects_population)):
				tournament_population = random.sample(self.chromosome_objects_population, k)
				tournament_population_fitness_scores = np.array([c.fitness(self.adj_mat, self.colors) for c in tournament_population])
				indices = np.argsort(tournament_population_fitness_scores)
				sorted_tournament_population_based_on_fitness_scores = [tournament_population[idx] for idx in indices]
				population_after_tournament.append(sorted_tournament_population_based_on_fitness_scores[0])
			population_after_tournament_fitness_scores = np.array([c.fitness(self.adj_mat, self.colors) for c in population_after_tournament])
			indices = np.argsort(population_after_tournament_fitness_scores)
			sorted_population_after_tournament = [population_after_tournament[idx] for idx in indices]
			for p in range(self.parents):
				population_after_selection.append(sorted_population_after_tournament[p].genes)
			self.population_after_selection = np.array(population_after_selection) # parents population
		else:
			raise NotImplementedError

	def __crossover(self):
		offspring = []
		if "point" in self.crossover_method:
			points = int(self.crossover_method.split("_")[0]) # point to crossover
			if points >= self.adj_mat.shape[0]-2:
				raise ValueError("\t❌ the point is too large")
			else:
				point_indices = random.sample(range(self.adj_mat.shape[0]), points) # getting random crossover lines
				# point_indices = [0, 1, 3]
				# p1 = [0 | 1 | 0 2 | 1]
				# p2 = [1 | 0 | 1 1 | 0]
				# o1 = [0 0 0 1 1]
				# o2 = [1 1 1 2 0]
				for curr_idx in range(self.parents):
					next_index = (curr_idx+1)%self.parents if next_index != 0 else: break
					first_parent  = self.population_after_selection[curr_idx]
					second_parent = self.population_after_selection[next_index]
					first_child   = np.zeros(self.adj_mat.shape[0])
					second_child  = np.zeros(self.adj_mat.shape[0])
					for line_idx in np.sort(point_indices):
						pass 
					offspring.append(first_child)
					offspring.append(second_child)
				self.population_after_crossover = np.array(offspring)
				

			self.population_after_crossover = np.array(offspring)
		elif self.crossover_method == "uniform": # =====================================================================================
			self.population_after_crossover = np.array(offspring)
			raise NotImplementedError # TODO
		elif self.crossover_method == "pmx": # =====================================================================================
			self.population_after_crossover = np.array(offspring)
			raise NotImplementedError # TODO
		elif self.crossover_method == "ox": # =====================================================================================
			self.population_after_crossover = np.array(offspring)
			raise NotImplementedError # TODO
		else:
			raise NotImplementedError


	def __mutation(self):
		offspring_after_mutation = []
		if self.mutation_method == "swap": # =====================================================================================
			self.population_after_mutation = np.array(offspring_after_mutation)
			raise NotImplementedError # TODO
		elif self.mutation_method == "creep": # =====================================================================================
			self.population_after_mutation = np.array(offspring_after_mutation)
			raise NotImplementedError # TODO
		elif self.mutation_rate == "interchanging": # =====================================================================================
			self.population_after_mutation = np.array(offspring_after_mutation)
			raise NotImplementedError # TODO
		elif self.mutation_rate == "reversing": # =====================================================================================
			self.population_after_mutation = np.array(offspring_after_mutation)
			raise NotImplementedError # TODO
		else:
			raise NotImplementedError

	def __replacement(self):
		population_next_generation = []
		if self.mutation_method == "generational_elitism":
			self.population = np.array(population_next_generation)
			raise NotImplementedError # TODO
		elif self.mutation_method == "generational_gap":
			self.population = np.array(population_next_generation)
			raise NotImplementedError # TODO
		elif self.mutation_rate == "steady_state":
			self.population = np.array(population_next_generation)
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
		print(self.best_chromosomes[0]) # the first one is the best one because its sorted in ascending order and is minimum

	def save(self):
		print(f"▶ Saving Best chromosomes and best scores")
		np.save(f"utils/best_chromo_in_{self.generations}_generations.npy", self.best_chromosomes)
		np.save(f"utils/best_fitness_scores_in_{self.generations}_generations.npy", self.best_fitness_scores)
