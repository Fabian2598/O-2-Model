2d and 3d O(2) model and 2d Ising model with the Worm Algorithm.
The programs compute the energy, the magnetic susceptibility and the specific heat. In
the 2d O(2) model the spin stiffness is computed as well.


The .ipynb files are testing versions written in python. They do work, however they are not as
fast as the C++ versions.

In order to run the C++ versions one has to compile 2DXY.cpp or 3DXY.cpp for a specific value of
L. Such a value can be changed in the .cpp files. The programs recieve the following inputs: 

-beta_min: starting value of $\beta$ = $\frac{1}{T}$ to run the simulations.

-beta_max: last value of $\beta$.

-Number of betas: Number of points between beta_min and beta_max to run simulations.

-Thermalization: Number of thermalization sweeps (e. g. 10 000).

-Measurements: Number of measurements.

-Step: Number of sweeps between measurements. 

The program will generate a .dat file with all the relevant information of the simulation for each beta.
In order to write files that contain the energy, specific heat, spin stifness, magnetic susceptibility 
and computing time one has to use the datan.cpp program.

The datan program needs the dimension of the model and a file with a list of names of the .dat files.
To easily create such a list one can write ls -1 *.dat > datfiles.txt in the terminal (linux). 
