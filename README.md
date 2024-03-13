# MolDynMD-jl
Simple experiments in molecular dynamics written in Julia.

## Dependencies of the Code

Install the following Julia package if not installd already

```
using Pkg
Pkg.add("Plots")
```

The code to run a simple simulation of <sup>1</sup>H<sup>35</sup>Cl is in the `run_hcl.qmd` Quarto Markdown file. To run it, you need to install Quarto and Quarto support for Julia. Find the instructions at:

[Using Julia](https://quarto.org/docs/computations/julia.html)

Once this is done, on macOS, execute Quarto from this repo's folder with

```
quarto preview /ABSOLUTE/PATH/TO/run_hcl.qmd
```

### On Making Quarto Work

Note: This requires that Quarto can find Jupyter and the IJulia Jupyter kernel. To install it, ake the following steps 

Install Jupyter Lab in a conda environment:

```
conda create -n MolDynMD-jl python=3.11
conda activate MolDynMD-jl
conda install jupyterlab
```

Then, with the 
