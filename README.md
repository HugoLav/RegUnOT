# RegUnOT

Code for the article "Entropy minimization with respect to branching Brownian motion is regularized unbalanced optimal transport" by A. Baradat and H. Lavenant

Article available at [PlaceHolder arxiv]. 

Implementation in Julia of the solver for the Regularized Unbalanced Optimal Transport Problem. It reproduces the figure that can be found in the article. The parameters and the boundary conditions are hard-coded in the file `generate_figures.jl`. Commenting relevant lines from 13 to 31 enables to choose to simulate OT, ROT, UOT or RUOT.

Then running `figure.tex` reproduces Figure 1.1 of the article. 