


#### 【Graph Coloring With Minimum Number of Colors Using Genetic Algorithm】

---

## 🔧 Setup

```console
pip install -r requirements.txt
```

> ⚠️ A default adjacency matrix is filled out in a text file called `adj_mat.txt` inside `utils` folder. You can create another one with the same structure of the default one from your graph and put it in there.

## 💻 Usage

```console
python colory.py --adj-mat utils/adj_mat.txt --colors o r b --chromosomes 50 --generations 15 --parents 30 --selection-method tournament --crossover-method 3_point --mutation-method creep --replacement-method generational_elitism --mutation-rate 0.20 --crossover-rate 0.80

```

> ⚠️ Run `python3 colory.py --help` for argument details.

## 📋 Procedures

#### 📌 Encoding

#### 📌 Objective functions