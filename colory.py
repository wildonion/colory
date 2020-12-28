




# coding: utf-8

'''
	Codded By : 
 █     █░ ██▓ ██▓    ▓█████▄  ▒█████   ███▄    █  ██▓ ▒█████   ███▄    █ 
▓█░ █ ░█░▓██▒▓██▒    ▒██▀ ██▌▒██▒  ██▒ ██ ▀█   █ ▓██▒▒██▒  ██▒ ██ ▀█   █ 
▒█░ █ ░█ ▒██▒▒██░    ░██   █▌▒██░  ██▒▓██  ▀█ ██▒▒██▒▒██░  ██▒▓██  ▀█ ██▒
░█░ █ ░█ ░██░▒██░    ░▓█▄   ▌▒██   ██░▓██▒  ▐▌██▒░██░▒██   ██░▓██▒  ▐▌██▒
░░██▒██▓ ░██░░██████▒░▒████▓ ░ ████▓▒░▒██░   ▓██░░██░░ ████▓▒░▒██░   ▓██
 -------------------------------------------------------------------------------------------------
| Feature Selection and Dimensionality Reduction using Genetic Algorithm For Breast Cancer Dataset
|-------------------------------------------------------------------------------------------------
|
| USAGE : 
|			python colory.py --chromosomes 200 --colors 3 \
|				 			 --generation 3 --parents 10 --selection-method roulette_wheel \
|				 			 --crossover-method multi_point --mutation-method flipping --mutation-rate 0.20
|			
|
|
|
|
| AVAILABLE GA OPERATORS METHODS: :
| 									selection methods -> roulette_wheel, tournament, rank
|									crossover methods -> multi_point, single_point, two_point
|									mutation methods  -> filipping, reversing, interchanging
|
|
|
'''


import networkx as nx
import os
import numpy as np
import argparse
from _evolver import population, genetic_process



# ------------ argument options
# ------------------------------
parser = argparse.ArgumentParser(description='Graph Coloring using Genetic Algorithm')
parser.add_argument('--chromosomes', action='store', type=int, help='The number of chromosomes in a population', required=True)
parser.add_argument('--colors', action='store', type=int, help='The number of minimum colors to paint with', required=True)
parser.add_argument('--generation', action='store', type=int, help='The number of generation', required=True)
parser.add_argument('--parents', action='store', type=int, help='The number of parents to mate for offspring', required=True)
parser.add_argument('--selection-method', action='store', type=str, help='Selection method for crossover operation', required=True)
parser.add_argument('--crossover-method', action='store', type=str, help='Crossover method to generate offspring', required=True)
parser.add_argument('--mutation-method', action='store', type=str, help='Mutation method to mutate offspring', required=True)
parser.add_argument('--mutation-rate', action='store', type=float, help='Mutation rate', required=True)
args = parser.parse_args()

BEST_CHROMOSOME_PATH = f"best_chromo_in_{args.generation}_generations.npy"
BEST_SCORE_PATH = f"best_scores_in_{args.generation}_generations.npy"





# ------------ initialize the population using defined arguments
# --------------------------------------------------------------------
pop = population(args.chromosomes, args.colors)
# print(f"second index of third population genes with length {len(pop[3])} is ::: {pop[3].genes_objects[2].allele}") # >>>>>>>>>>>> testing design pattern


if os.path.exists(BEST_CHROMOSOME_PATH) and os.path.exists(BEST_SCORE_PATH): # load saved chromosomes and scores
	best_chromosomes = np.load(BEST_CHROMOSOME_PATH)
	best_scores = np.load(BEST_SCORE_PATH)
else: # create a genetic process
	app = genetic_process(generation=args.generation, population=pop, parents=args.parents, selection_method=args.selection_method, 
					  	  crossover_method=args.crossover_method, mutation_method=args.mutation_method, mutation_rate=args.mutation_rate)
	app.run() # run the process
	app.plot(lib="plotly") # plot the result
	app.save() # save best chromosomes and scores
	best_chromosomes = app.best_chromosomes # all best chromosomes in every generation
	best_scores = app.best_scores # all best scores in every generation




# ------------ logging score infos after finishing the genetic process 
# ----------------------------------------------------------------------------
for i in range(len(best_scores)):
	print(f"🔬 Generation {i} Score --- {best_scores[i]:.0%}")
print()
print('▶ Average accepted score = ', np.mean(best_scores))
print('▶ Median score for accepted scores = ', np.median(best_scores))



# ------------ plotting colored graph using networkx
# --------------------------------------------------------------------