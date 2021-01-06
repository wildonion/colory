




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
| USAGE : 
|			python colory.py --adj-mat utils/adj_mat.txt --colors o r b \
|							 --chromosomes 50 --generations 15 --parents 10 \
|							 --selection-method tournament --crossover-method 2_point \
|							 --mutation-method creep --replacement-method generational_elitism \
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
parser.add_argument('--crossover-method', action='store', type=str, help='Crossover method to generate offspring (n_point[where n is an integer], uniform, pmx or ox).', required=True)
parser.add_argument('--mutation-method', action='store', type=str, help='Mutation method to mutate offspring (swap, creep, reversing or interchanging).', required=True)
parser.add_argument('--replacement-method', action='store', type=str, help='Replacement method to replace the old population (generational_elitism, generational_gap, or steady_state).', required=True)
parser.add_argument('--mutation-rate', action='store', type=float, help='Mutation rate (between 0.01 and 0.05 based on 20 <= chromosomes <= 30). You can use 1/chromosomes to setup the ratio.', required=True)
parser.add_argument('--crossover-rate', action='store', type=float, help='Crossover rate (between 0.75 and 0.95 based on 20 <= chromosomes <= 30).', required=True)
args = parser.parse_args()






# ------------ defining the paths and building adjacency matrix
# -------------------------------------------------------------------------
BEST_CHROMOSOMES_PATH = f"best_chromo_in_{args.generations}_generations.npy"
BEST_SCORES_PATH      = f"best_fitness_scores_in_{args.generations}_generations.npy"
ADJ_MAT               = np.array([list(map(lambda x : int(x), list(filter(lambda x: x != '', \
									[x for x in row.replace('\n', '').split(" ")])))) for row in args.adj_mat.readlines()])
COLORS                = np.array(args.colors)





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
	pop = population(COLORS, args.chromosomes, ADJ_MAT)





	# ------------ testing design patterns
	# -----------------------------------------
	print(f"\n\t🧬 third gene of third chromosome with length {len(pop[2])} is ::: {pop[2].gene_objects[2].allele}")







	# ------------ loading best chromosomes and best fitness scores
	# --------------------------------------------------------------------
	if os.path.exists(f"utils/+{BEST_CHROMOSOMES_PATH}") and os.path.exists(f"utils/+{BEST_SCORES_PATH}"): # load saved chromosomes and fitness scores
		print(f"\n‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗loading best chromosomes and fitness scores in each generation‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗\n")
		best_chromosomes = np.load(BEST_CHROMOSOMES_PATH)
		best_fitness_scores = np.load(BEST_SCORES_PATH)
	





	# ------------ running a genetic process to solve our problem
	# --------------------------------------------------------------------
	else: # create a genetic process
		print(f"\n‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗finding minimum valid colors for each node through a genetic process‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗\n")
		app = genetic_process(generations=args.generations, population=pop, parents=args.parents, selection_method=args.selection_method, 
						  	  crossover_method=args.crossover_method, mutation_method=args.mutation_method, replacement_method=args.replacement_method,
						  	  mutation_rate=args.mutation_rate, crossover_rate=args.crossover_rate)
		app.run() # run the process
		app.save() # save best chromosomes and fitness scores
		app.plot() # plot fitness score in each generation after finishing the process 
		app.draw() # draw the colored graph with best chromosome
		best_chromosomes = app.best_chromosomes # all best chromosomes in every generation
		best_fitness_scores = app.best_fitness_scores # all best fitness scores in every generation






	# ------------ logging statistical infos of the genetic process  
	# ------------------------------------------------------------------
	for i in range(len(best_fitness_scores)):
		print(f"\t🔬 Generation {i} Chromosome --- {best_chromosomes[i]}")
	print()
	print('\t▶ Average accepted score                   = ', np.mean(best_fitness_scores))
	print('\t▶ Median score for accepted fitness scores = ', np.median(best_fitness_scores))
