# The Semi-Lagrangian Multiscale Reconstruction Method (SLMsR)

This is the prototypical Julia implementation of the semil-Lagrangian basis reconstruction method *SLMsR*.

The aim is to demonstrate for the passive tracer transport problem with diffusion that it is
possible to reconstruct subgrid-features that can not be resolved 
in coarse-grid simulations locally into basis functions even when the
problem is dominated by advection. Standard multiscale methods for often fail in
this case.


The SLMsR code requires:

* A Linux distribution (we used Debian)
* *[Julia](www.julialang.com)* v1.3	
* *[Paraview](https://www.paraview.org/)* for the visualization

Get the code, change to the code directory and start the Julia REPL:

```
git clone https://github.com/konsim83/SLMsR-2D.git
cd SLMsR-2D
julia
```

To add all necessary Julia packages just run the <code>startup.jl</code> script 
in the root folder:

```@julia
include("startup.jl")
```



### Running the Examples

You can run (the first of the five) example codes through

```@julia
include("src/test1.jl")
```

This can take some time. Be patient since the code is not very efficient (yet).

The program will compute a low resolution solution, a high resolution reference solution
and a multiscale solution.

After the computations finished there will be a <code>data/</code> folder containing 
vtk files and other error files. You can inspect the files using Paraview.

You can inspect the relative (mean square) error of the solutions 
(or of the derivative) with respect to the reference solution using 

```@julia
Vis.plotErrorL2(errSTD, errMSR, t, "L2")
```

or 

```@julia
Vis.plotErrorL2(errSTD_H1, errMSR_H1, t, "H1")
```


---
**NOTE**

It is also possible to compare the data sets point-wise (for example using a Python
programmable filter) since they have been mapped to a common fine grid.
 
---