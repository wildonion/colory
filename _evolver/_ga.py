

import matplotlib.pyplot as plt
import networkx as nx
import plotly.graph_objects as go
import numpy as np
import random


__all__ = ['genetic_process', 'population'] # export only these classes



class gene:
	"""
		🔴 This is the gene class, an allele of a gene.


			 💡 Attributes:
					allele:
						type:
							public 
						value:
							- value of a gene
	"""
	def __init__(self, allele):
		self.allele = allele


# ≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣

class chromosome:
	"""
		🔴 This is the chromosome class, the heart of the genetic algorithm.
		⚠️ A valid chromosome (graph) is a chromosome that none of its two neighboring genes (nodes) have same colors.
		⚠️ A minimum valid chromosome (graph) is chromosome that is valid and its genes (nodes) colored with minimum number of colors. 


			 💡 Methods:
					fitness:
						type: 
							public
						params:
							- numpy array of adjacency matrix
							- numpy array of colors
						return: 
							- total fitness of the chromosome
						jobs:
							- calculate the total fitness of the chromosome 
					-------------------------------------------------------
					__invalid_genes_objective:
						type: 
							private
						params: 
							- numpy array of adjacency matrix
						return:
							- number of invalid genes inside the chromosome
						jobs:
							- calculate the number of invalid genes
					-------------------------------------------------------
					__minimum_genes_objective:
						type: 
							private 
						params:
							- numpy array of colors
						return:
							- minimum number of colors inside the chromosome
						jobs:
							- calculate the minimum number of colors
					-------------------------------------------------------
					__getitem__:
						type:
							magic
						params:
							- locus
						return:
							- a gene value by its locus
						jobs:
							- return a gene value by its locus 
					-------------------------------------------------------
					__len__:
						type:
							magic
						return:
							- length of the chromosome
						jobs:
							- return the length of the chromosome 

			 💡 Attributes:
					genes:
						type:
							public
						value:
							- numpy array of all gene values
					-------------------------------------------------------
					gene_objects:
						type:
							public
						value:
							- numpy array of all gene objects
					-------------------------------------------------------
					is_valid:
						type:
							public
						value:
							- is this a valid chromosome
					-------------------------------------------------------
					has_min_colors:
						type:
							public
						value:
							- does this chromosome have minimum number of colors

	"""
	def __init__(self, genes):
		self.genes = np.array(genes)
		self.gene_objects = np.array([gene(g) for g in self.genes])
		self.is_valid = False
		self.has_min_colors = False

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
			for j in range(self.gene_objects.shape[0]):
				if self.gene_objects[i].allele == self.gene_objects[j].allele:
					if adj_mat[i][j] == 1: # there is an edge between this node and the next node
						invalid += 1
		if invalid == 0: self.is_valid = True
		return invalid


	def __minimum_genes_objective(self, colors):
		_, unique_indices = np.unique(self.genes, return_index=True)
		unique_colors = self.genes[unique_indices] # all unique colors inside a genotype (choromosome)
		if unique_colors.shape[0] < colors.shape[0]: 
			m_prime = unique_colors.shape[0]
			self.has_min_colors = True
		else: 
			m_prime = colors.shape[0]
		return m_prime # minimum value is 1 and maximum value is colors.shape[0]

	def __getitem__(self, locus):
		return self.gene_objects[locus]

	def __len__(self):
		return self.genes.shape[0]


# ≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣

class population:
	"""
		🔴 This is the population class, a population of chromosomes.
		⚠️ We built our chromosomes as a list of random values in range 0 to the number of colors - 1.
		⚠️ If the chromosomes attribute is empty means that we are in first generation.
		⚠️ The first fitness scores in ascending_fitness_scores is the best chromosome.


			 💡 Methods:
					__init_pop:
						type: 
							private
						jobs:
							- generate chromosomes population
					-------------------------------------------------------
					fitness_scores:
						type:
							public
						return:
							- all chromosomes fitness scores in ascending order
							- the population of all chromosomes in ascending order
						jobs:
							- calculate each fitness score of all chromosomes
					-------------------------------------------------------
					get_adj_mat:
						type: 
							public
						return:
							- numpy array of adjacency matrix
						jobs:
							- return the numpy array of adjacency matrix
					-------------------------------------------------------
					get_colors:
						type: 
							public
						return:
							- numpy array of colors
						jobs:
							- return the numpy array of colors
					-------------------------------------------------------
					__call__:
						type: 
							magic
						return:
							- the population of all chromosomes
						jobs:
							- return the population of all chromosomes 
					-------------------------------------------------------
					__getitem__:
						type:
							magic
						params:
							- locus
						return:
							- a gene value by its locus
						jobs:
							- get the idx-th chromosome of the population 
					-------------------------------------------------------
					__len__:
						type:
							magic
						return:
							- length of the population
						jobs:
							- return the length of the population 

			 💡 Attributes:
					amount:
						type:
							public
						value:
							- number of chromosomes inside the population
					-------------------------------------------------------
					colors:
						type:
							public
						value:
							- numpy array of colors
					-------------------------------------------------------
					adj_mat:
						type:
							public
						value:
							- numpy array of adjacency matrix
					-------------------------------------------------------
					pop:
						type:
							public
						value:
							- numpy array of population of all chromosomes

	"""
	def __init__(self, amount=None, colors=None, adj_mat=None, chromosomes=[]):
		self.amount  = amount
		self.colors  = colors
		self.adj_mat = adj_mat
		if chromosomes == []: self.__init_pop()
		else: self.pop = [chromosome(c) for c in chromosomes]
	
	def __init_pop(self):
		self.pop = [chromosome([np.where(self.colors == np.random.choice(self.colors, 1))[0][0] for _ in range(self.adj_mat.shape[0])]) for _ in range(self.amount)] # build each chromosome genes with a random index of colors list

	def fitness_scores(self):
		fitness_scores = [chromosome.fitness(self.adj_mat, self.colors) for chromosome in self.pop] # all chromosomes fitness inside the generated population
		fitness_scores, genes_population = np.array(fitness_scores), np.array([c.genes for c in self.pop])
		indices = np.argsort(fitness_scores)
		ascending_fitness_scores = fitness_scores[indices]
		ascending_genes_population_based_on_fitness_scores = genes_population[indices]
		return ascending_fitness_scores, ascending_genes_population_based_on_fitness_scores

	def __len__(self):
		return self.pop.shape[0]

	def __getitem__(self, idx):
		return self.pop[idx]

	def __call__(self):
		return self.pop

	def get_adj_mat(self):
		return self.adj_mat

	def get_colors(self):
		return self.colors


# ≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣

class genetic_process():
	"""
		🔴 This is the genetic_process class, it controls the process of our genetic algorithm.
		⚠️ We built a numpy array called total_generations_valid_chromosomes of all generations that contains valid chromosomes.
		⚠️ We built a numpy array called total_generations_minimum_colors_valid_chromosomes of all generations that contains valid chromosomes with minimum number of colors.
		⚠️ We can't color the graph at all if there exist invalid chromosomes, means total_generations_valid_chromosomes.size == 0. 
		⚠️ We have to color the graph only using valid chromosomes if total_generations_minimum_colors_valid_chromosomes.size == 0.
		⚠️ We can color the graph using valid chromosomes with minimum number of colors if we have valid chromosomes or total_generations_valid_chromosomes.size != 0.
		⚠️ We Selected k chromosomes with uniform probability in our tournament method.


			 💡 Methods:
					run:
						type: 
							public
						jobs:
							- run our GA operators in every generation
							- build a numpy array of all generations that contains valid chromosomes
							- build a numpy array of all generations that contains valid chromosomes with minimum number of colors
					-------------------------------------------------------
					__selection:
						type:
							private
						jobs:
							- select best chromosomes based on defined method and the number of parents
					-------------------------------------------------------
					__single_point_crossover:
						type: 
							private
						params:
							- first parent
							- second parent
							- point to cross
						return:
							- two generated offspring 
						jobs:
							- do a single point crossover on selected parents
					-------------------------------------------------------
					__multi_point_crossover:
						type: 
							private
						params:
							- first parent
							- second parent
							- points to cross
						return:
							- two generated offspring 
						jobs:
							- do a recursive multi point crossover on selected parents using __single_point_crossover method
					-------------------------------------------------------
					__crossover:
						type:
							private
						jobs:
							- generate offsprings based on defined method
					-------------------------------------------------------
					__mutation:
						type:
							private
						jobs:
							- mutate offsprings based on defined method
					-------------------------------------------------------
					__replacement:
						type:
							private
						jobs:
							- replace the old population with new population based on defined method
					-------------------------------------------------------
					plot:
						type:
							public
						jobs:
							- plot the numpy array of best fitnesse scores collected in every generation using plotly
					-------------------------------------------------------
					draw:
						type:
							public
						jobs:
							- select a random chromosome either from valid chromosomes or valid chromosomes with minimum number of colors
							- color the graph using selected random chromosome and networkx
					-------------------------------------------------------
					save:
						type:
							public
						jobs:
							- save the numpy arrays of best fitnesse scores, valid chromosomes and valid chromosomes with minimum number of colors collected from every generation

			 💡 Attributes:
					colors:
						type:
							public
						value:
							- numpy array of colors
					-------------------------------------------------------
					adj_mat:
						type:
							public
						value:
							- numpy array of adjacency matrix
					-------------------------------------------------------
					graph:
						type:
							public
						value:
							- networkx graph from adjacency matrix
					-------------------------------------------------------
					generations:
						type:
							public
						value:
							- number of generations
					-------------------------------------------------------
					population:
						type:
							public
						value:
							- numpy array of population of all chromosome objects
					-------------------------------------------------------
					parents:
						type:
							public
						value:
							- number of parents to mate for breeding offspring
					-------------------------------------------------------
					mutation_rate:
						type:
							public
						value:
							- mutation ratio
					-------------------------------------------------------
					crossover_rate:
						type:
							public
						value:
							- crossover ratio
					-------------------------------------------------------
					alpha_rate:
						type:
							public
						value:
							- alpha ratio for replacement method
					-------------------------------------------------------
					genes_population_after_fitness:
						type:
							public
						value:
							- sorted chromosomes population with theit gene values
					-------------------------------------------------------
					population_after_selection:
						type:
							public
						value:
							- chromosomes population after selection process
					-------------------------------------------------------
					population_after_crossover:
						type:
							public
						value:
							- chromosomes population after crossover process
					-------------------------------------------------------
					population_after_mutation:
						type:
							public
						value:
							- chromosomes population after mutation process
					-------------------------------------------------------
					total_generations_valid_chromosomes:
						type:
							public
						value:
							- all none empty valid chromosomes collected from every generation
					-------------------------------------------------------
					total_generations_minimum_colors_valid_chromosomes:
						type:
							public
						value:
							- all none empty valid chromosomes with minimum number of colors collected from every generation
					-------------------------------------------------------
					best_fitness_scores:
						type:
							public
						value:
							- all best fitness scores collected from every generation


	"""
	def __init__(self, generations, population, parents, selection_method, crossover_method, mutation_method, alpha_rate, mutation_rate, crossover_rate):
		self.adj_mat 						  					= population.get_adj_mat()
		self.colors 						  					= population.get_colors()
		self.graph                                              = nx.from_numpy_matrix(self.adj_mat)
		self.generations 					  					= generations
		self.population 					  					= population
		self.parents 						  					= parents
		self.mutation_rate 					  					= mutation_rate
		self.crossover_rate 				  					= crossover_rate
		self.selection_method 				  					= selection_method
		self.crossover_method 				  					= crossover_method
		self.mutation_method 				  					= mutation_method
		self.alpha_rate 					  					= alpha_rate
		self.genes_population_after_fitness   					= None
		self.population_after_selection       					= None
		self.population_after_crossover       					= None
		self.population_after_mutation        					= None
		self.total_generations_valid_chromosomes                = []
		self.total_generations_minimum_colors_valid_chromosomes = []
		self.best_fitness_scores                                = []


	# ≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣
	
	def run(self):
		for g in range(self.generations):
			print("\n≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣\n")
			print(f"🧬 Generation --- {g+1}")
			fitness_scores, self.genes_population_after_fitness = self.population.fitness_scores()
			self.best_fitness_scores.append(fitness_scores[0])
			valid_chromo_in_this_generation     = np.array([chromosome.genes for chromosome in self.population() if chromosome.is_valid == True])
			min_valid_chromo_in_this_generation = np.array([chromosome.genes for chromosome in self.population() if chromosome.has_min_colors == True and chromosome.is_valid == True]) # a list of all best chromosomes in the world!
			if valid_chromo_in_this_generation.size != 0:
				self.total_generations_valid_chromosomes.append(valid_chromo_in_this_generation) # total valid chromosomes in one generation
			if min_valid_chromo_in_this_generation.size != 0:
				self.total_generations_minimum_colors_valid_chromosomes.append(min_valid_chromo_in_this_generation) # total valid chromosomes which have minimum colors in one generation
			print(f"\t▶  Population Shape --- {self.genes_population_after_fitness.shape}\n")
			print(f"\t▶  Best Fitness Scores of Two First Sorted Chromosomes --- {fitness_scores[:2]}\n") # minimum fitness scores are the best ones
			print(f"\t▶  Total Valid Chromosomes in this Generation --- {valid_chromo_in_this_generation.shape}\n")
			print(f"\t▶  Total Valid Chromosomes with Minimum Colors in this Generation --- {min_valid_chromo_in_this_generation.shape}\n")
			# =================== GA Operators ===================
			self.__selection() # select best chromosomes as parents
			self.__crossover() # parents mating pool
			self.__mutation() # mutating genes
			self.__replacement() # replacing old population
			# ====================================================
		self.total_generations_valid_chromosomes                = np.array(self.total_generations_valid_chromosomes, dtype=object)
		self.total_generations_minimum_colors_valid_chromosomes = np.array(self.total_generations_minimum_colors_valid_chromosomes, dtype=object)

	# ≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣
	# 									SELECTION OPERATORS

	def __selection(self):
		print(f"\t▶  Generating Parents Population using {self.selection_method} Method\n")
		population_after_selection = []
		# ==========================================================================
		#								ROULETTE WHEEL
		
		if self.parents <= self.genes_population_after_fitness.shape[0]: # check that the number of parents is smaller than the total population
			if self.selection_method == "roulette_wheel" or self.selection_method == "FPS":
				fitness_population = sum([c.fitness(self.adj_mat, self.colors) for c in self.population]) # sum of all scores (fitnesses)
				individual_expected_values = [(1/c.fitness(self.adj_mat, self.colors))/fitness_population for c in self.population] # all chromosomes prob (exprected values) - we scaled up every chromosome fitness cause big roulette of the wheel belongs to minimum fitnesses 
				cum_prob = [sum(individual_expected_values[:i+1]) for i in range(len(individual_expected_values))] # cumulative sum of chromosomes exprected values (prob)
				for i in range(self.parents):
					r = random.random()
					if cum_prob[i] >= r: # add this chromosome only if its cum_prob is greater than the generated random number
						population_after_selection.append(self.population[i].genes)
				self.population_after_selection = np.array(population_after_selection) # parents population
				
			# ==========================================================================
			#								  RANK
			
			elif self.selection_method == "rank":
				for p in range(self.parents): # the first self.parents chromosomes are the best ones
					population_after_selection.append(self.genes_population_after_fitness[p])
				self.population_after_selection = np.array(population_after_selection) # parents population
			
			# ==========================================================================
			#								TOURNAMENT
			
			elif self.selection_method == "tournament":
				k = int(np.log2(len(self.population())))
				population_after_tournament = []
				for _ in range(len(self.population())):
					tournament_population = random.sample(self.population(), k) # none repetitive k indices with uniform probability
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
			print(f"\t▶  Population Shape After Selection --- {self.population_after_selection.shape}\n")
		else:
			raise ValueError("\t❌ number of parents is greater than the total population")
	
	# ≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣

	def __single_point_crossover(self, p1, p2, point):
		first_child = np.append(p1[:point], p2[point:])
		second_child = np.append(p2[:point], p1[point:])
		return first_child, second_child


	def __multi_point_crossover(self, p1, p2, points):
		for p in points:
			p1, p2 = self.__single_point_crossover(p1, p2, p)
		return p1, p2
	
	# ≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣
	# 									CROSSOVER OPERATORS
	
	def __crossover(self):
		parents_indices = random.sample(range(self.population_after_selection.shape[0]), 2) # none repetitive indices
		parents         = self.population_after_selection[parents_indices]
		print(f"\t▶  Selected Parents to Mate --- {list(parents)}\n")
		offspring = []
		# ==========================================================================
		#								T-POINT
		
		if "point" in self.crossover_method:
			points = int(self.crossover_method.split("-")[0]) # point to crossover
			if points > self.adj_mat.shape[0]-2:
				raise ValueError("\t❌ the point is too large")
			else:
				do_cross = random.random() <= self.crossover_rate
				point_indices   = random.sample(range(self.adj_mat.shape[0]), points) # getting random crossover lines - none repetitive indices
				print(f"\t\t▶  Point Indices --- {point_indices}\n")
				if len(point_indices) == 1 and do_cross:
					first_child, second_child = self.__single_point_crossover(parents[0], parents[1], point_indices[0])
				elif len(point_indices) >= 2 and do_cross:
					first_child, second_child = self.__multi_point_crossover(parents[0], parents[1], point_indices)
				else:
					first_child, second_child = parents[1], parents[0]  
				offspring.append(first_child)
				offspring.append(second_child)
				_msg = "Generated" if do_cross else "No Need to Generate" 
				print(f"\t\t▶  {_msg} Offspring using t-point Crossover --- {offspring}\n")
				self.population_after_crossover = np.vstack((self.population_after_selection, np.array(offspring)))
		# ==========================================================================
		#								UNIFORM
		
		elif self.crossover_method == "uniform":
			do_cross 		  = random.random() <= self.crossover_rate
			swap_probability  = 0.5
			first_child       = np.zeros(self.adj_mat.shape[0], dtype=int)
			second_child      = np.zeros(self.adj_mat.shape[0], dtype=int)
			chromosome_length = self.adj_mat.shape[0] 
			if do_cross:
				for g in range(chromosome_length):
					u = random.random()
					if u <= swap_probability:
						first_child[g]  = parents[1][g]
						second_child[g] = parents[0][g]
					else:
						first_child[g]  = parents[0][g]
						second_child[g] = parents[1][g]
			else:
				first_child, second_child = parents[1], parents[0]
			offspring.append(first_child)
			offspring.append(second_child)
			_msg = "Generated" if do_cross else "No Need to Generate" 
			print(f"\t\t▶  {_msg} Offspring using uniform Crossover --- {offspring}\n")
			self.population_after_crossover = np.vstack((self.population_after_selection, np.array(offspring)))
		# ==========================================================================
		else:
			raise NotImplementedError
		# ==========================================================================
		print(f"\t▶  Population Shape After Crossing Over --- {self.population_after_crossover.shape}\n")

	# ≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣
	# 									MUTATION OPERATORS
	
	def __mutation(self):
		offspring_after_mutation = []
		for p in range(self.population_after_crossover.shape[0]):
			do_mutate                  = random.random() <= self.mutation_rate
			offspring_before_mutation  = self.population_after_crossover[p]
			print(f"\t\t▶  Offspring Before {self.mutation_method} --- {offspring_before_mutation}")
			# ==========================================================================
			# 									SWAP
			
			if self.mutation_method == "swap":
				locuses               = random.sample(range(self.adj_mat.shape[0]), 2) # none repetitive indices
				print(f"\t\t▶  First Locus --- {locuses[0]}")
				print(f"\t\t▶  Second Locus --- {locuses[1]}")
				first_gene_allele     = self.population_after_crossover[p][locuses[0]]
				print(f"\t\t▶  First Allele --- {first_gene_allele}")
				second_gene_allele    = self.population_after_crossover[p][locuses[1]]
				print(f"\t\t▶  Second Allele --- {second_gene_allele}\n")
				if do_mutate:
					self.population_after_crossover[p][locuses[0]] = second_gene_allele
					self.population_after_crossover[p][locuses[1]] = first_gene_allele
				else:	
					pass
			# ==========================================================================
			# 									CREEP
			
			elif self.mutation_method == "creep":
				locus                  = random.sample(range(self.adj_mat.shape[0]), 1)[0] # none repetitive indices
				print(f"\t\t▶  Selected Locus --- {locus}")
				selected_gene_allele   = self.population_after_crossover[p][locus]
				print(f"\t\t▶  Selected Gene Allele --- {selected_gene_allele}\n")
				if do_mutate:
					self.population_after_crossover[p][locus] = random.sample(range(self.colors.shape[0]), 1)[0] # none repetitive indices
				else: pass
			# ==========================================================================
			# 								   INVERSION
			
			elif self.mutation_method == "inversion":
				inversion_indices = np.sort(random.sample(range(self.population_after_crossover[p].shape[0]), 2)) # none repetitive indices
				print(f"\t\t▶  Inversion Indices --- {inversion_indices}")
				selected_genes_allele    = self.population_after_crossover[p][inversion_indices[0]:inversion_indices[1]]
				print(f"\t\t▶  Selected Genes --- {selected_genes_allele}")
				inverted_genes_allele    = selected_genes_allele[::-1]
				print(f"\t\t▶  Inverted Genes --- {inverted_genes_allele}\n")
				if do_mutate:
					self.population_after_crossover[p][inversion_indices[0]:inversion_indices[1]] = inverted_genes_allele
				else: pass
			# ==========================================================================
			else:
				raise NotImplementedError
			# ==========================================================================
			_msg = "Mutated" if do_mutate else "No Need to Mutate"
			print(f"\t\t▶  {_msg} Offspring using {self.mutation_method} Mutation --- {self.population_after_crossover[p]}")
			print("\t\t——————————————————————————————\n")
			offspring_after_mutation.append(self.population_after_crossover[p])
		self.population_after_mutation = np.array(offspring_after_mutation)
		print(f"\t▶  Population Shape After Mutation --- {self.population_after_mutation.shape}\n")

	# ≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣
	# 									REPLACEMENT OPERATORS

	def __replacement(self):
		print(f"\t▶  Replacing with {self.alpha_rate*100} ratio\n")
		# ==========================================================================
		#								   STEADY STATE
		
		if self.alpha_rate == 0.01:
			print(f"\t▶  Steady State Method is Selected")
			index                             = random.sample(range(self.population_after_mutation.shape[0]), 1)[0]
			print(f"\t\t▶  Random Selected Index --- {index}")
			worst_chromosome                  = self.genes_population_after_fitness[-1] # the last one is the worst cause it's sorted
			print(f"\t\t▶  Worst Chromosome from Old Generation --- {worst_chromosome}")
			random_mutated_chromosome         = self.population_after_mutation[index]
			print(f"\t\t▶  Random Selected Mutated Chromosome --- {random_mutated_chromosome}")
			worst_chromosome_fitness          = chromosome(worst_chromosome).fitness(self.adj_mat, self.colors)
			print(f"\t\t▶  Worst Chromosome Fitness --- {worst_chromosome_fitness}")
			random_mutated_chromosome_fitness = chromosome(random_mutated_chromosome).fitness(self.adj_mat, self.colors)
			print(f"\t\t▶  Random Selected Chromosome Fitness --- {random_mutated_chromosome_fitness}")
			if random_mutated_chromosome_fitness <= worst_chromosome_fitness:
				self.genes_population_after_fitness[-1] = random_mutated_chromosome
				print(f"\t\t▶  Replaced Worst Chromosome with Random Selected Chromosome --- {self.genes_population_after_fitness[-1]}\n")
			else:
				print(f"\t\t▶  No Need to Use Steady State Replacement Method\n")
			np.random.shuffle(self.genes_population_after_fitness)
			self.population = population(colors=self.colors, adj_mat=self.adj_mat, chromosomes=self.genes_population_after_fitness)
		# ==========================================================================
		#									GENERATIONAL
		
		elif self.alpha_rate == 1:
			print(f"\t▶  Generational+Elitism Method is Selected\n")
			old_population_best_chromosome = self.genes_population_after_fitness[0]
			population_next_generation = np.vstack((self.population_after_mutation, old_population_best_chromosome))
			self.population = population(colors=self.colors, adj_mat=self.adj_mat, chromosomes=population_next_generation)
		# ==========================================================================
		#								  GENERATIONAL GAP
		
		else:
			print(f"\t▶  Generational Gap Method is Selected")
			number_of_random_indices = random.sample(range(self.population_after_mutation.shape[0]), int(self.alpha_rate * self.population_after_mutation.shape[0]))
			print(f"\t\t▶  Selected {self.alpha_rate*100}% of Mutated Population Chromosomes --- {list(self.population_after_mutation[number_of_random_indices])}")
			print(f"\t\t▶  Selected {self.alpha_rate*100}% of Old Population Chromosomes --- {list(self.genes_population_after_fitness[number_of_random_indices])}")
			self.genes_population_after_fitness[number_of_random_indices] = self.population_after_mutation[number_of_random_indices]
			print(f"\t\t▶  Replaced {self.alpha_rate*100}% of Old Population Chromosomes with {self.alpha_rate*100}% of Mutated Population Chromosomes --- {list(self.genes_population_after_fitness[number_of_random_indices])}\n")
			np.random.shuffle(self.genes_population_after_fitness)
			self.population = population(colors=self.colors, adj_mat=self.adj_mat, chromosomes=self.genes_population_after_fitness)
		# ==========================================================================
		print("\t▶  Population Shape After Replacement --- ({}, {})\n".format(len(self.population()), len(self.population()[0]) ))

	# ≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣
	# 										  PLOT FITNESS GENERATIONS

	def plot(self):
		fig = go.Figure()
		fig.add_trace(go.Scatter(x=np.arange(1, self.generations+1), y=self.best_fitness_scores, mode='lines', name='Fitness Generation'))
		fig.update_xaxes(title_text='Generation')
		fig.update_yaxes(title_text='Best Fitness')
		fig.show()

	# ≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣
	# 		COLOR THE GRAPH USING VALID CHROMOSOMES OR CHROMOSOMES WITH MINIMUM NUMBER OF COLORS 

	def draw(self):
		if self.total_generations_valid_chromosomes.size == 0: # we can't color the graph cause there exist invalid chromosomes
			print(f"\t\t❌  It's not Possible to Color the graph with the Given Number of Colors. Detected Invalid Chromosomes\n")
		if self.total_generations_minimum_colors_valid_chromosomes.size == 0 and self.total_generations_valid_chromosomes.size != 0: # we have to plot only if, first all chromosomes are valid and second there are minimum colors inside their genes
			print(f"\t\t❌  It's not Possible to Color the Graph with Less than {self.colors.shape[0]} Colors")
			print(f"\t\t✅  It's Possible to Color the Graph with {self.colors.shape[0]} Colors using Random Selected Valid Chromosome")
			generation_index                     = random.sample(range(self.total_generations_valid_chromosomes.shape[0]), 1)[0]
			selected_generation_chromosome_index = random.sample(range(self.total_generations_valid_chromosomes[generation_index].shape[0]), 1)[0]
			selected_valid_chromosome            = self.total_generations_valid_chromosomes[generation_index][selected_generation_chromosome_index]
			print(f"\t\t🔮  Random Selected Valid Chromosome --- {selected_valid_chromosome}")
			colors = self.colors[selected_valid_chromosome]
			nx.draw(self.graph, with_labels=True, node_color=colors, node_size=700, alpha=0.9)
			plt.savefig(f"utils/results/latest-test-case/colored_graph_using_valid_chromosomes_after_{self.generations}_generations.png")
			print(f"\t\t📸  Colored Graph Saved at utils/results/latest-test-case\n")

		if self.total_generations_minimum_colors_valid_chromosomes.size != 0: # we can color the graph with minimum colors cause there are some valid chromosomes
			print(f"\t\t✅  It's Possible to Color the Graph with Less than {self.colors.shape[0]} Colors using Random Selected Valid Chromosome with Minimum Colors")
			generation_index    = random.sample(range(self.total_generations_minimum_colors_valid_chromosomes.shape[0]), 1)[0]
			selected_generation = self.total_generations_minimum_colors_valid_chromosomes[generation_index]
			selected_generation_nb_colors_for_each_chromosome = np.array([c[np.unique(c, return_index=True)[1]].shape[0] for c in selected_generation])
			minimum_nb_colors   = np.min(selected_generation_nb_colors_for_each_chromosome)
			nb_color_index      = random.sample(range(selected_generation_nb_colors_for_each_chromosome.shape[0]), 1)
			selected_nb_color   = selected_generation_nb_colors_for_each_chromosome[nb_color_index]
			while selected_nb_color != minimum_nb_colors:
				nb_color_index    = random.sample(range(selected_generation_nb_colors_for_each_chromosome.shape[0]), 1)
				selected_nb_color = selected_generation_nb_colors_for_each_chromosome[nb_color_index]
			selected_minimum_colors_valid_chromosome = selected_generation[nb_color_index][0]
			print(f"\t\t🔮  Random Selected Valid Chromosome with Minimum Colors --- {selected_minimum_colors_valid_chromosome}")
			minimum_colors = self.colors[selected_minimum_colors_valid_chromosome]
			nx.draw(self.graph, with_labels=True, node_color=minimum_colors, node_size=700, alpha=0.9)
			plt.savefig(f"utils/results/latest-test-case/colored_graph_using_minimum_colors_after_{self.generations}_generations.png")
			print(f"\t\t📸  Colored Graph Saved at utils/results/latest-test-case\n")

	# ≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣
	
	def save(self):
		print("\n≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣\n")
		print(f"⏳ Saving Valid Chromosomes, Chromosomes with Minimum Colors and Best Fitness Scores")
		np.save(f"utils/results/latest-test-case/valid_chromosomes_in_{self.generations}_generations.npy", self.total_generations_valid_chromosomes)
		np.save(f"utils/results/latest-test-case/minimum_colors_valid_chromosomes_in_{self.generations}_generations.npy", self.total_generations_minimum_colors_valid_chromosomes)
		np.save(f"utils/results/latest-test-case/best_fitness_scores_in_{self.generations}_generations.npy", self.best_fitness_scores)
		print()