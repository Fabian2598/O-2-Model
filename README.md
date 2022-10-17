# 2d and 3d O(2) model and 2d Ising model with the Worm Algorithm.
### The programs compute the energy, the magnetic susceptibility and the specific heat. In the 2d O(2) model the spin stiffness is computed as well. 

We work with square lattices of size $L^2$ for 2d and with cubic lattices of size $L^3$ for 3d. The errors are computed by using the jackknife method.

Some basic explanations of the worm algorithm can be found in the PDF files [WAIsingModel.pdf](WAIsingModel.pdf) and [WAXYModel.pdf](WAXYModel.pdf).
Preliminar results are in the [Results](O(2)Model/Results) folder.

The .ipynb files are testing versions written in python. They do work, however they are not as fast as the C++ versions.

In order to run the C++ versions one has to compile [2DXY.cpp](O(2)Model/2DXY.cpp) or [3DXY.cpp](O(2)Model/3DXY.cpp) for a specific value of $L$.
```console
g++ 2DXY.cpp -o 2DXYL8.x
```
Such a value can be changed at the beginning of the .cpp files. 
```cpp
double pi = 3.14159265359;
double E = 0, En = 0, En2 = 0, chi = 0, Cv = 0, rho = 0;
double dchi, dCv, dE, dE2, drho;
double R; //Acceptance ratio
double Wx = 0, Wy = 0;

int L = 8; //This value must be changed to compute simulations for different lattices. 
```
The programs recieve the following inputs: 

* beta_min: starting value of $\beta$ = $\frac{1}{T}$ to run the simulations.

* beta_max: last value of $\beta$.

* Number of betas: Number of points between beta_min and beta_max to run simulations.

* Thermalization: Number of thermalization sweeps (*e. g.* 10 000).

* Measurements: Number of measurements.

* Step: Number of sweeps between measurements. 

Example:
```console
./2DXYL8.x
beta min: 0.1
beta max: 1.4
Number of betas: 30
Thermalization: 10000
Measurements: 10000
Step (sweeps between measurements): 10
``` 

The program will generate a .dat file with all the relevant information of the simulation for each beta. In order to write files that contain the energy, specific heat, spin stiffness, magnetic susceptibility and computing time one has to use the [datan.cpp](O(2)Model/datan.cpp) program.

The [datan.cpp](O(2)Model/datan.cpp) program needs the dimension of the model and a file with a list of names of the .dat files. To easily create such a list one can write for instance

```console
ls -1 *.dat > datfiles.txt
``` 
An example of how to use datan is the following 
```console
ls -1 *.dat > datfiles.txt
g++ datan.cpp -o datan.x
./datan.x
Dimension of the model (2 or 3): 2
File with list of names: datfiles.txt
``` 
The program generates .txt files with the results written in columns. Their detailed structure is

* 2DXY_Cv_L8_Meas10000.txt: $\beta$, $C_V$, $C_V$ error.
* 2DXY_Chi_L8_Meas10000.txt $\beta$, $\chi$, $\chi$ error.
* 2DXY_E_L8_Meas10000.txt: $\beta$, $E$, $E$ error.
* 2DXY_Time_L8_Meas10000: $\beta$, computing time in seconds.
* 2DXY_Rho_L8_Meas10000: $\beta$, $\rho_s$, $\rho_s$ error. 
