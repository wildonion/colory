




# coding: utf-8

'''
	Codded By : 
 █     █░ ██▓ ██▓    ▓█████▄  ▒█████   ███▄    █  ██▓ ▒█████   ███▄    █ 
▓█░ █ ░█░▓██▒▓██▒    ▒██▀ ██▌▒██▒  ██▒ ██ ▀█   █ ▓██▒▒██▒  ██▒ ██ ▀█   █ 
▒█░ █ ░█ ▒██▒▒██░    ░██   █▌▒██░  ██▒▓██  ▀█ ██▒▒██▒▒██░  ██▒▓██  ▀█ ██▒
░█░ █ ░█ ░██░▒██░    ░▓█▄   ▌▒██   ██░▓██▒  ▐▌██▒░██░▒██   ██░▓██▒  ▐▌██▒
░░██▒██▓ ░██░░██████▒░▒████▓ ░ ████▓▒░▒██░   ▓██░░██░░ ████▓▒░▒██░   ▓██
 -------------------------------------------------------------------------------------------------
| Graph Coloring With Minimum Number of Colors Using Genetic Algorithm
|-------------------------------------------------------------------------------------------------
|
| USAGE EXAMPLE : 
|			python colory.py --adj-mat utils/matrices/tree/adj_mat.txt --colors o r b g \
|							 --chromosomes 50 --generations 60 --parents 35 \
|							 --selection-method tournament --crossover-method 3-point \
|							 --mutation-method creep --alpha-rate 20 \
|							 --mutation-rate 0.20 --crossover-rate 0.80
|			
|
|
|
|
|
|
'''


import os, sys
import numpy as np
import argparse
from _evolver import population, genetic_process





# ------------ argument options
# -------------------------------
parser = argparse.ArgumentParser(description='【Graph Coloring using Genetic Algorithm】')
parser.add_argument('--adj-mat', help='Path to adjacency matrix', type=argparse.FileType('r', encoding='UTF-8'), required=True)
parser.add_argument('--colors', action='store', nargs="+", help='Name of colors separated by space.', required=True)
parser.add_argument('--chromosomes', action='store', type=int, help='The number of total chromosomes in a population.', required=True)
parser.add_argument('--generations', action='store', type=int, help='The number of generations.', required=True)
parser.add_argument('--parents', action='store', type=int, help='The number of parents to mate for breeding offspring.', required=True)
parser.add_argument('--selection-method', action='store', type=str, help='Selection method for crossover operation (roulette_wheel, tournament or rank).', required=True)
parser.add_argument('--crossover-method', action='store', type=str, help='Crossover method to generate offspring (n_point[where n is an integer] or uniform).', required=True)
parser.add_argument('--mutation-method', action='store', type=str, help='Mutation method to mutate offspring (swap, creep or inversion).', required=True)
parser.add_argument('--alpha-rate', action='store', type=float, help='Alpha rate for replacing the old population.', required=True)
parser.add_argument('--mutation-rate', action='store', type=float, help='Mutation rate (between 0.01 and 0.05 based on 20 <= chromosomes <= 30). You can use 1/chromosomes to setup the ratio.', required=True)
parser.add_argument('--crossover-rate', action='store', type=float, help='Crossover rate (between 0.75 and 0.95 based on 20 <= chromosomes <= 30).', required=True)
args = parser.parse_args()






# ------------ building adjacency matrix
# -----------------------------------------------
ADJ_MAT                               = np.array([list(map(lambda x : int(x), list(filter(lambda x: x != '', \
													[x for x in row.replace('\n', '').split(" ")])))) for row in args.adj_mat.readlines()])
COLORS                                = np.array(args.colors)





# ------------ it's not possible to paint the graph using available colors
# --------------------------------------------------------------------------------------
if COLORS.shape[0] == len(ADJ_MAT) or COLORS.shape[0] > len(ADJ_MAT):
	print("\t❌ no need to use GA, you can paint each node with a specific color\n")
	sys.exit(1)





# ------------ we can paint the graph using available colors
# ----------------------------------------------------------------
elif COLORS.shape[0] < len(ADJ_MAT):
	



	# ------------ initialize the population using defined arguments
	# --------------------------------------------------------------------
	pop = population(args.chromosomes, COLORS, ADJ_MAT)





	# ------------ testing design patterns
	# -----------------------------------------
	print(f"\n\t🧬 third gene of third chromosome with length {len(pop[2])} is ::: {pop[2].gene_objects[2].allele}")

	





	# ------------ running a genetic process to solve our problem
	# --------------------------------------------------------------------
	print(f"\n‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗finding minimum valid colors for each node through a genetic process‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗\n")
	app = genetic_process(generations=args.generations, population=pop, parents=args.parents, selection_method=args.selection_method, 
						  crossover_method=args.crossover_method, mutation_method=args.mutation_method, alpha_rate=args.alpha_rate,
						  mutation_rate=args.mutation_rate, crossover_rate=args.crossover_rate)
	app.run() # run the process
	app.save() # save valid chromosomes, minimum valid chromosomes and best fitness scores
	# app.plot() # plot fitness score in each generation after finishing the process 
	app.draw() # draw the colored graph with best chromosome
	best_fitness_scores 			 				   = app.best_fitness_scores # all best fitness scores collected from every generation
	total_generations_valid_chromosomes 			   = app.total_generations_valid_chromosomes # all valid chromosomes collected from every generation
	total_generations_minimum_colors_valid_chromosomes = app.total_generations_minimum_colors_valid_chromosomes # all valid chromosomes with minimum colors collected from every generation





	'''
	≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣
		 genetic process public attributes

	app.adj_mat
	app.colors
	app.graph    
	app.generations
	app.population    
	app.parents
	app.mutation_rate
	app.crossover_rate
	app.selection_method
	app.crossover_method
	app.mutation_method
	app.alpha_rate
	app.genes_population_after_fitness
	app.population_after_selection
	app.population_after_crossover
	app.population_after_mutation
	app.total_generations_valid_chromosomes
	app.total_generations_minimum_colors_valid_chromosomes
	app.best_fitness_scores

	≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣≣
	'''






# ------------ logging statistical infos of the genetic process  
# ------------------------------------------------------------------
print('\t▶ Average Fitness Scores --- ', np.mean(best_fitness_scores))
print('\t▶ Median Fitness Scores --- ', np.median(best_fitness_scores))
if total_generations_valid_chromosomes.size > 0: # if there are some valid chromosomes 
	print(f"\t▶ Total Valid Chromosomes in Generation 1 --- {total_generations_valid_chromosomes[0].shape}")
if total_generations_minimum_colors_valid_chromosomes.size > 0: # if there are some valid chromosomes with minimum colors
	print(f"\t▶ Total Valid Chromosomes with Minimum Colors in Generation 1 --- {total_generations_minimum_colors_valid_chromosomes[0].shape}\n")
