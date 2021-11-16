# Percolation with homeostatic plasticity

Code to reproduce the main findings of the paper *"Percolation in networks equipped with local homeostatic plasticity"*.

The code uses the libraries `networkx` and `seaborn`, which can be installed by running
``` bash
pip install networkx
```
``` bash
pip install seaborn
```
## Example of successive instances

`FIG2.py` runs 3 successive instances of the damage-response process on a given directed network. 

The syntax is:

``` python
python3 FIG2.py [path_edgelist] [weight_distribution] [samples] [outstring] [show_plot]
```

where:

* `path_edgelist`:  `\path\to\edgelist`. 

* `weight_distribution`: assigns a weight to each link extracted from a weight distribution. The options availables are `uniform` and `gauss`.

* `samples`: number of subdivision for the x-axis.

* `outstring`: string for the output file. The output file is always generated and stored locally with the `.npy` format.

* `Plot`: `plot` for visualizing the result of the simulation, `noplot` otherwise. 

To reproduce the results shown in Figure 2 run:

``` python
python3 FIG2.py datasets/RG_N5000_k4.txt gauss 5000 out_fig2 plot
```

Note that `RG_N5000_k4.txt` has the same topology of the network used in Figure 2 in the paper, but the size is 5k nodes instead of 50k.

## Evolution of the Weak Largest Connected Component

`FIG3.py` runs 20 successive instances of both the damage process only and the damage-response process on a given directed network and returns the evolution of the weak largest connected component for the two cases. 

Syntax and input parameters are the same of `FIG2.py`.

To reproduce the results shown in Figure 3b run:

``` python
python3 FIG3.py datasets/mouse_kasthuri_directed.txt uniform 5000 out_fig3 plot
```

## Analytical prediction of average weight, average degree and weight-degree correlation

`FIG4.py` runs 50 single instances of the damage-response process at increasing threshold values and returns the values of average weight, average degree and weight degree covariance. Those quantities can be evluated analitically for a single instance of the process. 

Syntax and input parameters are the same of `FIG2.py`.

To reproduce the results shown in Figure 4b and Figure 6a run:

``` python
python3 FIG4.py datasets/RG_N5000_k4.txt uniform 5000 out_fig4 plot
```
## Evolution of the leading eigenvalue

`FIG5.py` runs 50 single instances and 50 successive instances of the damage-response process and returns the leading eigenvalue of the adjacency matrix. In case of single instances the model provides and analytical approximation.

Syntax and input parameters are the same of `FIG2.py`.

To reproduce the results shown in Figure 5a and Figure 5b run:

``` python
python3 FIG5.py datasets/RG_N400_k10.txt uniform 3000 out_fig5 plot
```

Note that `RG_N5000_k4.txt` has the same topology of the network used in Figure 5 in the paper, but the size is 400 nodes instead of 1000. Furthermore the default number of instances of the Montecarlo simulation are 25 instead of 200.