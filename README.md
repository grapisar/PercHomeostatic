# Percolation with homeostatic plasticity

Code to reproduce the findings of the paper *"Percolation in networks equipped with local homeostatic plasticity"*.

## Single Instance

`example_simple.py` runs three successive instances of the damage-response process on a given directed network. The default setting reproduces the findings shown in Fig.2, but on a smaller network, which is a N=5k nodes Random Regular graph with degree k=4.

The syntax is:

``` python
python3 example_single.py [path_edgelist] [weight_distribution] [samples] [outstring] [show_plot]
```

where:

* `path_edgelist`:  `\path\to\edgelist`. 

* `weight_distribution`: assigns a weight to each link extracted from a weight distribution. The options availables are `uniform` and `gauss`.

* `samples`: number of subdivision for the x-axis.

* `outstring`: string for the output file.

* `Plot`: `plot` for visualizing the result of the simulation, `noplot` otherwise.

### Example

``` python
python3 example_single.py RG_N5000_k4.txt gauss 10000 out plot
```
