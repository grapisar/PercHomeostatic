# Percolation with homeostatic plasticity

Code to reproduce the findings on percolation in networks equipped with local homeostatic response mechanism.

## Single Instance

`example_single.py` runs a single instance of the damage-removal process on a given directed network.
The syntax is:

``` python
python3 example_single.py [path_edgelist] [weight_distribution] [samples] [outstring] [show_plot]
```

where:

* `path_edgelist`:  `\path\to\edgelist`.

* `weight_distribution`: assigns a weight to each link extracted from a weight distribution. The options availables are `uniform` and `gauss`.

* `samples`: number of subdivision for the x-axis.

* `outstring`: string for the output file.

* `plot`: write PLOT if you want to visualize the results of the simulation.

### Example

``` python
python3 example_single.py RG_N5000_k4.txt gauss 10000 out PLOT
```