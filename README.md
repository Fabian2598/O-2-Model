# Monte Carlo simulations of the 2d and 3d O(2) model and the 2d Ising model.
## The programs compute the energy E, magnetization $M$, magnetic susceptibility $\chi_M$, specific heat $C_V$ and the vortex density $\rho_V$ (for the O(2) model). 
## The files have C++ implementations of the Worm, cluster, Metropolis and heatbath algorithm. 

* The worm algorithm implementation is based on: N. Prokof'ev and B. Svistunov. Worm Algorithms for Classical Statistical Models. Phys. Rev. Lett., 87:160601, (2001).

* The cluster algorithm implementation is based on: R. H. Swendsen and J.-S. Wang. Nonuniversal critical dynamics in Monte Carlo simulations. Phys. Rev. Lett., 58:86-88, (1987). J. Hoshen and R. Kopelman. Percolation and cluster distribution. I. Cluster multiple labeling technique and critical concentration algorithm. Phys. Rev. B, 14:3438-3445, (1976).

* The heatbath algorithm implementation is based on: T. Hattori and H. Nakajima. Improvement of efficiency in generating random U(1) variables with Boltzmann distribution. Nucl. Phys. B Proc. Suppl., 26:635-637, (1992).

We work on square lattices of size $L^d. The errors are computed by using the jackknife method.

Some basic explanations of the worm algorithm can be found in the PDF files [WAIsingModel.pdf](WAIsingModel.pdf) and [WAXYModel.pdf](WAXYModel.pdf).
Results for simulations in equilibrium for the O(2) model are in the [Results](O(2)Model/Results) folder.

### Instructions to compile and run the codes

The organization of the files could be confusing, but the structure to use them is similar. In the following I detailed how to use the codes. 

Each folder has three types of files: something that begins with 2DXY or 3DXY, statistics.h and datan.cpp. To better explain how to run the programs consider the [O(2)Model/Equilibrium/2Dimensions](O(2)Model/Equilibrium/2Dimensions) folder, which has several files that begin with 2DXY, followed by the initials of the algorithm. In addition, a file named statistics.h and datan.cpp are present as well. In order to compile the 2DXY files, one writes, for instance.
```console
g++ 2DXY.cpp -o 2DXYL8.x
```
All 2DXY or 3DXY files have the following line at some part in the beginning
```cpp
constexpr int L = 8; 
```
This can be changed, before compiling, to modify the volume. The statistics.h file has to be in the same directory too, since it cotains all the functions needed to do the statistical analysis. This file is repeated in many folders for convenience, but it is always the same.

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

The program will generate .dat files with all the relevant information of the simulation ($\beta$, $L$, computing time, $E$, $C_V$, $M$, $\chi_M$), which is printed in the terminal as well. If one wants to better understand the structure of such files, it is recommended to revise the last few lines of the 2DXY or 3DXY codes. However, the datan.cpp program takes that information and puts it into separate files for each observable, which makes it convenient. Its compilation does not need of extra files.

To use the program, [O(2)Model/Equilibrium/2Dimensions/datan.cpp](O(2)Model/Equilibrium/2Dimensions/datan.cpp), do the following;
```console
ls -1 *.dat > datfiles.txt 
``` 
This creates a file with a list containing the name of the *.dat files. Then compile and run the program
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
* 2DXY_M_L8_Meas10000.txt $\beta$, $M$, $M$ error.
* 2DXY_E_L8_Meas10000.txt: $\beta$, $E$, $E$ error.
* 2DXY_Time_L8_Meas10000: $\beta$, computing time in seconds.
* 2DXY_Vort_L8_Meas10000: $\beta$, $\rho_V$, $\rho_V$ error. (vortex density) 
* 2DXY_AVort_L8_Meas10000: $\beta$, $\rho_AV$, $\rho_AV$ error. (anti-vortex density. It is equal to -$\rho_V$) 

For the worm algorithm, one has to compile [O(2)Model/Equilibrium/2Dimensions/datanW.cpp](O(2)Model/Equilibrium/2Dimensions/datanW.cpp) instead, because the observables computed are not the same. I did not find a way to determine $\rho_V$ with that algorithm, but the spin sitfness is simple to measure with that algorithm, so it is considered in the program.




Some results are in the [Results](O(2)Model/Results) folder.
