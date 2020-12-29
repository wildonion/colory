




# coding: utf-8

'''
	Codded By : 
 █     █░ ██▓ ██▓    ▓█████▄  ▒█████   ███▄    █  ██▓ ▒█████   ███▄    █ 
▓█░ █ ░█░▓██▒▓██▒    ▒██▀ ██▌▒██▒  ██▒ ██ ▀█   █ ▓██▒▒██▒  ██▒ ██ ▀█   █ 
▒█░ █ ░█ ▒██▒▒██░    ░██   █▌▒██░  ██▒▓██  ▀█ ██▒▒██▒▒██░  ██▒▓██  ▀█ ██▒
░█░ █ ░█ ░██░▒██░    ░▓█▄   ▌▒██   ██░▓██▒  ▐▌██▒░██░▒██   ██░▓██▒  ▐▌██▒
░░██▒██▓ ░██░░██████▒░▒████▓ ░ ████▓▒░▒██░   ▓██░░██░░ ████▓▒░▒██░   ▓██
 -------------------------------------------------------------------------------------------------
| Graph Coloring With Minimum Number of Color Using Genetic Algorithm
|-------------------------------------------------------------------------------------------------
|
| USAGE : 
|			python colory.py --adj-mat adj_mat.txt --colors o r b \
|							 --chromosomes 200 --generations 15 --parents 10 \
|							 --selection-method roulette_wheel --crossover-method multi_point \
|							 --mutation-method flipping --mutation-rate 0.20 --corssover-rate 0.80
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
import os, sys
import numpy as np
import argparse
from _evolver import population, genetic_process





# ------------ argument options
# -------------------------------
parser = argparse.ArgumentParser(description='Graph Coloring using Genetic Algorithm')
parser.add_argument('--adj-mat', help='Path to adjacency matrix', type=argparse.FileType('r', encoding='UTF-8'), required=True)
parser.add_argument('--colors', action='store', nargs="+", help='Name of colors separated by space', required=True)
parser.add_argument('--chromosomes', action='store', type=int, help='The number of total chromosomes in a population', required=True)
parser.add_argument('--generations', action='store', type=int, help='The number of generations', required=True)
parser.add_argument('--parents', action='store', type=int, help='The number of parents to mate for breeding offspring', required=True)
parser.add_argument('--selection-method', action='store', type=str, help='Selection method for crossover operation', required=True)
parser.add_argument('--crossover-method', action='store', type=str, help='Crossover method to generate offspring', required=True)
parser.add_argument('--mutation-method', action='store', type=str, help='Mutation method to mutate offspring', required=True)
parser.add_argument('--mutation-rate', action='store', type=float, help='Mutation rate', required=True)
parser.add_argument('--crossover-rate', action='store', type=float, help='Crossover rate', required=True)
args = parser.parse_args()






# ------------ defining the paths and building adjacency matrix
# -------------------------------------------------------------------------
BEST_CHROMOSOME_PATH = f"best_chromo_in_{args.generation}_generations.npy"
BEST_SCORE_PATH      = f"best_scores_in_{args.generation}_generations.npy"
ADJ_MAT              = np.array([list(map(lambda x : int(x), list(filter(lambda x: x != '', \
									[x for x in row.replace('\n', '').split(" ")])))) for row in args.adj_mat.readlines()])






# ------------ it's not possible to paint the graph using available colors
# --------------------------------------------------------------------------------------
if len(args.colors) == len(ADJ_MAT) or len(args.colors) > len(ADJ_MAT):
	print("\t❌ no need to use GA, you can paint each node with a specific color\n")
	sys.exit(1)





# ------------ we can paint the graph using available colors
# ----------------------------------------------------------------
elif len(args.colors) < len(ADJ_MAT):
	print(f"\n‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗solving the problem by using GA‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗‗\n")





	# ------------ initialize the population using defined arguments
	# --------------------------------------------------------------------
	pop = population(args.chromosomes, args.colors)






	# ------------ testing design patterns
	# -----------------------------------------
	print(f"second index of third population genes with length {len(pop[3])} is ::: {pop[3].genes_objects[2].allele}")







	# ------------ loading best chromosomes and best scores
	# --------------------------------------------------------------------
	if os.path.exists(BEST_CHROMOSOME_PATH) and os.path.exists(BEST_SCORE_PATH): # load saved chromosomes and scores
		best_chromosomes = np.load(BEST_CHROMOSOME_PATH)
		best_scores = np.load(BEST_SCORE_PATH)
	





	# ------------ running a genetic process to solve our problem
	# --------------------------------------------------------------------
	else: # create a genetic process
		app = genetic_process(generations=args.generations, population=pop, parents=args.parents, selection_method=args.selection_method, 
						  	  crossover_method=args.crossover_method, mutation_method=args.mutation_method, 
						  	  mutation_rate=args.mutation_rate, crossover_rate=args.crossover_rate)
		app.run() # run the process
		app.plot(lib="plotly") # plot the result using selected library
		app.save() # save best chromosomes and scores
		best_chromosomes = app.best_chromosomes # all best chromosomes in every generation
		best_scores = app.best_scores # all best scores in every generation






	# ------------ logging statistical infos of the genetic process  
	# ------------------------------------------------------------------
	for i in range(len(best_scores)):
		print(f"🔬 Generation {i} Score --- {best_scores[i]:.0%}")
	print()
	print('▶ Average accepted score = ', np.mean(best_scores))
	print('▶ Median score for accepted scores = ', np.median(best_scores))







# ------------ assign each node a specific color and plot the colored graph using networkx
# -----------------------------------------------------------------------------------------------