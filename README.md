


#### 【Graph Coloring With Minimum Number of Colors Using Genetic Algorithm】


## 🔧 Setup

```console
pip install requirements.txt
```

## 💻 Usage

```console
python3 colory.py --adj-mat utils/adj_mat.txt --colors o r b --chromosomes 200 --generations 15 \
				  --parents 10 --selection-method roulette_wheel --crossover-method multi_point \ 
				  --mutation-method flipping --mutation-rate 0.20 --crossover-rate 0.80

```

> ⚠️ Run `python3 colory.py --help` for argument details.

> ⚠️ A default adjacency matrix is filled out in a text file called `adj_mat.txt` inside `utils` folder. You can create another one with the same structure of the default one from your graph and put it in there.

## 📋 Procedures