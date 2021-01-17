

## 🔧 Setup

```console
pip install -r requirements.txt
```

> ⚠️ A default adjacency matrix is filled out in a text file called `adj_mat.txt` for each structure inside `matrices` folder. You can create another one with the same structure of the default one from your graph and put it in its related directory.

## 💻 Usage Example

```console
python colory.py --adj-mat utils/matrices/tree/adj_mat.txt --colors orange red blue green --chromosomes 50 --generations 20 --parents 30 --selection-method tournament --crossover-method 3-point --mutation-method creep --alpha-rate 0.30 --mutation-rate 0.20 --crossover-rate 0.80

```

> ⚠️ Run `python colory.py --help` for argument details.

## 📋 Procedures

#### 📌 Encoding & Initial Population

> Chromosome Representation of a Colored Graph
<p align="center">
    <img src="https://github.com/wildonion/colory/blob/main/utils/coloring_chromo.png">
</p>

We represent the graph using adjacency matrices. The user gives the colors.

There is no need to use **GA** if the number of colors is more than or equal to the number of nodes. If so, an error will pop up, showing that we do not need to use **GA**.


The encoding that we used is discrete. We have an array filled with the numbers mapped to each color, and the length of the array is equal to the nodes' number.

⚠️ Note that the user have to use the full name of the colors as input. Refer to the setup section.

We use a uniform random function to create the initial population using the given adjacency matrices and colors. The user gives the total number of the population.

#### 📌 Objective Functions

The obtained fitness function consists of two parts: The first part detects the `invalid_genes` and tries to minimize the number of them. In the second part, we have a variable called `m_prime`, which is the minimum number of colors used for coloring the graph. We need to use constant values for both parts of the objective function to show which part is more important. Alpha and beta are the constant variables. Assume that the sum of Alpha and beta is equal to one. The issue is to color the graph using the minimum valid colors. In order to find the proper multi-objective function, both `invalid_genes` and `m_prime` should have the same range. To fix this, we need to normalize the values as shown below:

`invalid_normalized = invalid_genes/ total edges`

`m_prime_normalized = (minimum_number_of_colors - 1) / colors`

`total_fitness = (alpha*invalid_normalized) + (beta*m_prime_normalized)`


[Fitness Function](https://github.com/wildonion/colory/blob/e6e94342b2c72e49019bbc9ff3f7a22580e6eea4/_evolver/_ga.py#L121)

[Objective Function part 1](https://github.com/wildonion/colory/blob/e6e94342b2c72e49019bbc9ff3f7a22580e6eea4/_evolver/_ga.py#L132)

[Objective Function part 2](https://github.com/wildonion/colory/blob/e6e94342b2c72e49019bbc9ff3f7a22580e6eea4/_evolver/_ga.py#L143)

⚠️ Depending on the method that we have used during the process, the fitness function could be equal to the objective function or inverted.

#### 📌 Optimization

We have a two-part goal for this problem. The first and the essential part is to find the valid coloring for the graph in which none of the connected nodes have the same color; the second part is to color the graph using the minimum number of colors. Since **GA** is a stochastic search algorithm based on natural competition principles between individuals, it is possible not to reach the intended minimum colors for each graph depending on the problem inputs. Nonetheless, the achieved result is that the minimum number of colors reached up to the assumed generation with the number of chromosomes and parents using each operation's selected methods. The code is entirely operative.

#### 📌 Genetic Operators

Order of operations are as follows, respectively:

Selection

* As selection methods, we have used the roulette wheel (**FPS**), rank, and tournament.

 * In the roulette wheel (**FPS**), the fitness function is equal to the inversed objective function.
 * In tournament and ranked based selection, the fitness function is equal to the objective function.
 * In the tournament method, we used equal probability distribution for each gene in the selected chromosome.

Crossover

* As crossover methods, we have used t-point and uniform.

Mutation

* As mutation methods, we have use swap, creep, and inversion.

Replacement

* As replacement methods, we have used steady-state, generational gap, and generational method using elitism depending on the proportion of the replaced population defined by the alpha rate, which should be initialized by the user.

⚠️ The user should initialize crossover and mutation rates. As we have used the **GA**, the crossover rate should be much higher than the mutation rate.


#### 📌 Stopping Criterion

Since the total fitness function never equals zero (The minimum number of colors is equal to one.), we considered the stopping criterion to reach the maximum number of generations.

#### 📌 Coloring the Graph

To color the graph, we have defined two **NumPy** arrays and have stored the valid chromosomes and the valid chromosomes with the minimum number of colors in each generation called `total_generations_valid_chromosomes` and `total_generations_minimum_colors_valid_chromosomes`, respectively.  Using these two arrays, we plot the colored graph based on the our coloring problem inputs' size.

[Draw Function](https://github.com/wildonion/colory/blob/e6e94342b2c72e49019bbc9ff3f7a22580e6eea4/_evolver/_ga.py#L774)

## 📊 Results

We ran the **GA** for different test-cases and store their results inside [`utils/results`](https://github.com/wildonion/colory/blob/main/utils/results/) folder. Below are two of the test-cases that our algorithm reached state-of-the-art of solving coloring problem. 🙂

### Test-case-1

[Argument Options](https://github.com/wildonion/colory/blob/main/utils/results/test-case-1/arguments.txt)

[Best Fitness Scores in 60 Generations](https://github.com/wildonion/colory/blob/main/utils/results/test-case-1/best_fitness_scores_in_60_generations.npy)

[Valid Chromosomes with Minimum Colors in 60 Generations](https://github.com/wildonion/colory/blob/main/utils/results/test-case-1/minimum_colors_valid_chromosomes_in_60_generations.npy)

[Valid Chromosomes in 60 Generations](https://github.com/wildonion/colory/blob/main/utils/results/test-case-1/valid_chromosomes_in_60_generations.npy)

> Statistical Logs
<p align="center">
    <img src="https://github.com/wildonion/colory/blob/main/utils/results/test-case-1/stat_log.png">
</p>

> Fitness Generations
<p align="center">
    <img src="https://github.com/wildonion/colory/blob/main/utils/results/test-case-1/fitness_generations.png">
</p>

> Colored Graph using Valid Chromosomes
<p align="center">
    <img src="https://github.com/wildonion/colory/blob/main/utils/results/test-case-1/colored_graph_using_valid_chromosomes_after_60_generations.png">
</p>

### Test-case-2

[Argument Options](https://github.com/wildonion/colory/blob/main/utils/results/test-case-2/arguments.txt)

[Best Fitness Scores in 60 Generations](https://github.com/wildonion/colory/blob/main/utils/results/test-case-2/best_fitness_scores_in_60_generations.npy)

[Valid Chromosomes with Minimum Colors in 60 Generations](https://github.com/wildonion/colory/blob/main/utils/results/test-case-2/minimum_colors_valid_chromosomes_in_60_generations.npy)

[Valid Chromosomes in 60 Generations](https://github.com/wildonion/colory/blob/main/utils/results/test-case-2/valid_chromosomes_in_60_generations.npy)

> Statistical Logs
<p align="center">
    <img src="https://github.com/wildonion/colory/blob/main/utils/results/test-case-2/stat_log.png">
</p>

> Fitness Generations
<p align="center">
    <img src="https://github.com/wildonion/colory/blob/main/utils/results/test-case-2/fitness_generations.png">
</p>

> Colored Graph using Valid Chromosomes with Minimum Colors
<p align="center">
    <img src="https://github.com/wildonion/colory/blob/main/utils/results/test-case-2/colored_graph_using_minimum_colors_after_60_generations.png">
</p>
