

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

#### 📌 Objective Functions

#### 📌 Optimization

#### 📌 Genetic Operators

#### 📌 Stopping Criterion


## 📊 Results

We ran the GA for different test-cases and store their results inside [`utils/results`](https://github.com/wildonion/colory/blob/main/utils/results/) folder. Below are two of the test-cases that our algorithm reached state-of-the-art of solving coloring problem. 🙂

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
