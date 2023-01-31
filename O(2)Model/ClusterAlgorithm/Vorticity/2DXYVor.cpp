#include <iostream>
#include <vector> 
#include <random> 
#include <time.h> 
#include <fstream>
#include <ctime>
#include <cmath> 
#include <algorithm>
#include "statistics.h"

double pi = 3.14159265359;
double E, Vort, Antvort;
double dE, dVort, dAntvort;
double p;//Acceptance ratio and reflected angle
int label = 0;
double phi; //Spin angle

constexpr int L = 10;
constexpr  int maxsize = L*L;
int CLabels[maxsize]; //Cluster labels.
std::vector<std::vector<double>> LatticeX(L, std::vector<double>(L,0)); 
std::vector<std::vector<double>> LatticeY(L, std::vector<double>(L,0));
std::vector<std::vector<double>> LatticeRX(L, std::vector<double>(L,0));
std::vector<std::vector<double>> LatticeRY(L, std::vector<double>(L,0));

std::vector<std::vector<int>> xBonds(L, std::vector<int>(L,0));
std::vector<std::vector<int>> yBonds(L, std::vector<int>(L,0));
std::vector<std::vector<int>> Labels(L, std::vector<int>(L,0)); //Labels of each site.

inline double correction(double a){    
    if (a>pi){
        return (-2*pi + a);
    }
    return a;
}

inline void angle(double x, double y){
    double norm = sqrt(x*x+ y*y);
    //Quadrants I and II//
    if (asin(y/norm) > 1e-10){
        phi = acos(x/norm);
    }
    //Quadrants III and IV
    else if (asin(y/norm) < -1e-10){
        phi = 2*pi - acos(x/norm);
    }
    else{
        phi = 0;
    }
    //This returns an angle between 0 and 2*pi
}

inline void Vorticity(){
    //Vorticity of a plaquette delta(i,i+x) + delta(i+x,i+x+y) + delta(i+x+y, i+y) + delta(i+y, i)
    //delta(i,j) = theta(j) - theta(i)
    for(int i = 0; i<L; i++){
        for(int j = 0; j<L; j++){
            angle(LatticeX[i][j], LatticeY[i][j]); double theta_i = phi;
            angle(LatticeX[i][modulo(j+1,L)], LatticeY[i][modulo(j+1,L)]); double theta_ix = phi;
            angle(LatticeX[modulo(i-1,L)][modulo(j+1,L)], LatticeY[modulo(i-1,L)][modulo(j+1,L)]); double theta_ixy = phi;
            angle(LatticeX[modulo(i-1,L)][j], LatticeY[modulo(i-1,L)][j]); double theta_iy = phi;
            double q = 1/(2*pi) * ( correction(fmodulo(theta_ix-theta_i,2*pi)) + 
            correction(fmodulo(theta_ixy-theta_ix,2*pi)) + 
            correction(fmodulo(theta_iy-theta_ixy,2*pi)) +
            correction(fmodulo(theta_i-theta_iy,2*pi)) );
            if ( absVal(q-1) <= 1e-6) {Vort += 1;}
            else if ( absVal(q+1) <= 1e-6){Antvort += 1;}
        }
    }  
}
    
//-----Computes the energy of a configuration-----//
inline void initialize_lattice(){
    for(int i = 0; i<L; i++){
        for(int j =0; j<L; j++){
            double r_init = rand_range(0,2*pi);
            LatticeX[i][j] = cos(r_init);
            LatticeY[i][j] = sin(r_init);
        }
    }
}

inline void EnergyCA(){
    for (int i = 0; i<L; i++){
        for (int j = 0; j<L; j++){
            E += -LatticeX[i][j]*LatticeX[i][modulo(j+1,L)] - LatticeY[i][j]*LatticeY[i][modulo(j+1,L)]
                 -LatticeX[i][j]*LatticeX[modulo(i+1,L)][j] - LatticeY[i][j]*LatticeY[modulo(i+1,L)][j];
        }
    }
}

//-----Function that generates bonds between neighbouring sites-----//
inline void Bonds(double beta){
    // modifies xBonds (bonds in the x direction)
    //          yBonds (bonds in the y direction)
    //          SpinLatticeR (Lattice with the angles of the reflected spins)
    double rphi = rand_range(0.0, 2*pi); //Angle of the r vector
    double rx = cos(rphi), ry = sin(rphi); //vector associated to the rphi angle
    for(int i = 0; i<L; i++){
        for(int j = 0; j<L; j++){
            double Sdot = LatticeX[i][j]*rx + LatticeY[i][j]*ry; //Dot product of S with r
            LatticeRX[i][j] = LatticeX[i][j] - 2*rx*Sdot; //Reflected spin.
            LatticeRY[i][j] = LatticeY[i][j] - 2*ry*Sdot;
            
            //---Creating bond with the right neigbour---//
            double SNdot = rx*LatticeX[i][modulo(j+1,L)] + ry*LatticeY[i][modulo(j+1,L)]; //Dot product of the right neighbour with r
            p = 0;
            if (SNdot*Sdot >=0){ p = 1-exp(-2*beta*SNdot*Sdot);}
            double R = ((double) rand() / (RAND_MAX));
            if (R<p){xBonds[i][j] = 1;}

            //---Creating bond with the lower neigbour---//
            SNdot = rx*LatticeX[modulo(i+1,L)][j] + ry*LatticeY[modulo(i+1,L)][j];
            p = 0;
            if (SNdot*Sdot >=0){p = 1-exp(-2*beta*SNdot*Sdot);}
            R = ((double) rand() / (RAND_MAX));
            if (R<p){yBonds[i][j] = 1;}
        }
    }
}
   
inline int find(int x){
    int y = x;
    while (CLabels[y] != y){
        y = CLabels[y];
    }
    while (CLabels[x] != x){
        int z = CLabels[x];
        CLabels[x] = y;
        x = z;
    }
    return y;
}

//-----Hoshen-Kopelman algorithm-----//
inline void HoshenKopelman(){
    for(int i = 0; i<L; i++){
        for(int j = 0; j<L; j++){
            std::vector<int> LBond(4), UBond(4), pLabels;
            int bonds = 0, pLabel;
            if (i>0 && yBonds[i-1][j] == 1){
                UBond[bonds] = i - 1;
                LBond[bonds] = j;
                bonds += 1;
                pLabel = find(Labels[UBond[bonds-1]][LBond[bonds-1]]);
                pLabels.push_back(pLabel);
            }
            if (i == L-1 && yBonds[i][j] == 1){
                UBond[bonds] = 0;
                LBond[bonds] = j;
                bonds += 1;
                pLabel = find(Labels[UBond[bonds-1]][LBond[bonds-1]]);
                pLabels.push_back(pLabel);
            }
            if (j >0 && xBonds[i][j-1] == 1){
                UBond[bonds] = i;
                LBond[bonds] = j-1;
                bonds += 1;
                pLabel = find(Labels[UBond[bonds-1]][LBond[bonds-1]]);
                pLabels.push_back(pLabel);
            }
            if (j == L-1 && xBonds[i][j] == 1){
                UBond[bonds] = i;
                LBond[bonds] == 0;
                bonds += 1;
                pLabel = find(Labels[UBond[bonds-1]][LBond[bonds-1]]);
                pLabels.push_back(pLabel);
            }

            if (bonds == 0){
                Labels[i][j] = label;
                CLabels[label] = label;
                label += 1;
            }
            else{

                int minLabel = *std::min_element(pLabels.begin(), pLabels.end());
                Labels[i][j] = minLabel;
                for (int b = 0; b<bonds; b++){
                    pLabel = pLabels[b];
                    CLabels[pLabel] = minLabel;
                    Labels[UBond[b]][LBond[b]] = minLabel;
                }
            }
        }
    }
    for(int i = 0; i<L; i++){
        for(int j = 0; j<L; j++){
          Labels[i][j] =  find(Labels[i][j]);
        }
    }
            
}

inline void flip(){
    std::vector<double> probs;
    for(int i = 0; i<label; i++){
        double R = ((double) rand() / (RAND_MAX));
        probs.push_back(R);
    }
    for(int i = 0; i<L; i++){
        for(int j = 0; j<L; j++){                  
            if (probs[Labels[i][j]] < 0.5){
                LatticeX[i][j] = LatticeRX[i][j];
                LatticeY[i][j] = LatticeRY[i][j];
            }
        }
    }

}

inline void reset(){
    for(int i = 0; i<L; i++){
        for(int j=0; j<L; j++){
            xBonds[i][j] = 0; yBonds[i][j] = 0; Labels[i][j] = 0; LatticeRX[i][j] = 0; LatticeRY[i][j] = 0;
        }
    }
    for(int i = 0; i<L*L; i++){
        CLabels[i] = 0;
    }
    label = 0;
}    

void CA_XY2d(double beta, int Ntherm, int Nmeas, int Nsteps){
    std::vector<double> Energy(Nmeas), VORT(Nmeas), AVORT(Nmeas);
    //Thermalization//
    initialize_lattice();
    for(int i = 0; i<Ntherm; i++){
        Bonds(beta); //Computes the bonds and the reflected lattice.
        HoshenKopelman(); //Identifies the clusters.
        flip(); //Flips the spins with probability 1/2.
        reset();//Resets bonds and labels.
    }
    for(int i = 0; i<Nmeas; i++){
        Bonds(beta);
        HoshenKopelman(); 
        flip(); 
        reset();
        
        EnergyCA(); 
        Vorticity(); 
        Energy[i] = E; 
        VORT[i] = Vort;
        AVORT[i] = Antvort;
        E = 0; Vort = 0; Antvort = 0;
        for(int j = 0; j<Nsteps; j++){
            Bonds(beta);
            HoshenKopelman(); 
            flip(); 
            reset(); 
        }
    }
    E = mean(Energy); dE = Jackknife_error(Energy, 20);
    Vort = mean(VORT)/(L*L); dVort = Jackknife_error(VORT,20)/(L*L);
    Antvort = mean(AVORT)/(L*L); dAntvort = Jackknife_error(AVORT,20)/(L*L);
}

//-----------------------------------------//
int main(){    
srand(time(0));
    int Ntherm, Nmeas, Nsteps, Nbeta;
    double beta_min, beta_max; 
    //---Input data---//
    std::cout << "beta min: "; 
    std::cin >> beta_min;
    std::cout << "beta max: ";
    std::cin >> beta_max;
    std::cout << "Number of betas: ";
    std::cin >> Nbeta;
    std::cout << "Thermalization: ";
    std::cin >> Ntherm;
    std::cout << "Measurements: ";
    std::cin >> Nmeas;
    std::cout << "Step (sweeps between measurements): ";
    std::cin >> Nsteps;
    std::cout << " " << std::endl;
    std::vector<double> Betas(Nbeta);
    if (Nbeta == 1){ 
        Betas = {beta_min}; 
    }
    else{ 
        Betas = linspace(beta_min,beta_max, Nbeta);
    }
    char NameData[50], Data_str[100];

    for (int i = 0; i < Nbeta; i++) {
        int A = Betas[i];
        int beta= (Betas[i] - A) * 10000;
        if (Betas[i] < 1){sprintf(NameData, "2DXY_L%d_Meas%d_b0%d.dat", L, Nmeas,beta);}
        else{beta = beta + (int) Betas[i] * 10000; sprintf(NameData, "2DXY_L%d_Meas%d_b%d.dat", L, Nmeas,beta);}

        std::ofstream Datfile;
        Datfile.open(NameData);
        
        clock_t begin = clock();
        std::cout << "beta = " << Betas[i] << "  T = " << 1/Betas[i] << std::endl;
        CA_XY2d(Betas[i], Ntherm, Nmeas, Nsteps);
        sprintf(Data_str,"%-30d%-30d%-30d%-30d%-30.17g\n",Ntherm,Nmeas,Nsteps,L,Betas[i]);
        Datfile << Data_str;
        sprintf(Data_str, "%-30.17g%-30.17g\n", E, dE);
        Datfile << Data_str;
        sprintf(Data_str, "%-30.17g%-30.17g\n", Vort, dVort);
        Datfile << Data_str;
        sprintf(Data_str, "%-30.17g%-30.17g\n", Antvort, dAntvort);
        Datfile << Data_str;

        std::cout << "E = " << E << " +- " << dE << std::endl;
        std::cout << "V = " << Vort << " +- " << dVort << std::endl;
        std::cout << "A = " << Antvort << " +- " << dAntvort << std::endl;
        E = 0; dE= 0; Vort = 0; dVort=0; Antvort = 0; dAntvort = 0;
        //----Computing time----//
        clock_t end = clock();
        double elapsed_secs = double(end - begin) / CLOCKS_PER_SEC;
        sprintf(Data_str, "%-30.17g", elapsed_secs);
        Datfile << Data_str; 
        std::cout << "Time = " << elapsed_secs << " s" << std::endl;
        std::cout << "------------------------------" << std::endl;
    }
    return 0;       
}
