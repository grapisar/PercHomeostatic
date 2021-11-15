# Percolation with homeostatic plasticity

Code to reproduce the main findings of the paper *"Percolation in networks equipped with local homeostatic plasticity"*.

## Example of successive instances

`FIG2.py` runs three successive instances of the damage-response process on a given directed network. 

The syntax is:

``` python
python3 FIG2.py [path_edgelist] [weight_distribution] [samples] [outstring] [show_plot]
```

where:

* `path_edgelist`:  `\path\to\edgelist`. 

* `weight_distribution`: assigns a weight to each link extracted from a weight distribution. The options availables are `uniform` and `gauss`.

* `samples`: number of subdivision for the x-axis.

* `outstring`: string for the output file.

* `Plot`: `plot` for visualizing the result of the simulation, `noplot` otherwise.

To reproduce the results shown in Figure 2 run:

``` python
python3 FIG2.py datasets/RG_N5000_k4.txt gauss 10000 out_fig2 plot
```

`RG_N5000_k4.txt` has the same topology of the network used in Figure 2 in the paper, but the size is 5k nodes instead of 50k.
