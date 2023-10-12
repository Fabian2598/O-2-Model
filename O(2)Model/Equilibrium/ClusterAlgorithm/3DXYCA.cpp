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
double E = 0, En2 = 0, M = 0, M2=0, chi = 0, Cv = 0;
double dchi, dM, dM2, dCv, dE, dE2;
double p; //Acceptance ratio.
int label = 0;

constexpr int L = 12;
constexpr  int maxsize = L*L*L;
int CLabels[maxsize]; //Cluster labels.

//We store the spin components.
double LatticeX[L][L][L]; 
double LatticeY[L][L][L]; 
double LatticeRX[L][L][L]; 
double LatticeRY[L][L][L]; 

double xBonds[L][L][L]; 
double yBonds[L][L][L]; 
double zBonds[L][L][L]; 
double Labels[L][L][L]; 


inline void initialize_lattice(){
    for(int i = 0; i<L; i++){
        for(int j =0; j<L; j++){
            for(int k = 0; k<L; k++){
                double r_init = rand_range(0,2*pi);
                LatticeX[i][j][k] = cos(r_init);
                LatticeY[i][j][k] = sin(r_init);
            }
        }
    }
}

inline void EnergyCA(){
    for (int i = 0; i<L; i++){
        for (int j = 0; j<L; j++){
            for(int k = 0; k<L; k++ ){
                E += -LatticeX[i][j][k]*LatticeX[i][modulo(j+1,L)][k] - LatticeY[i][j][k]*LatticeY[i][modulo(j+1,L)][k]
                 -LatticeX[i][j][k]*LatticeX[modulo(i+1,L)][j][k] - LatticeY[i][j][k]*LatticeY[modulo(i+1,L)][j][k]
                 -LatticeX[i][j][k]*LatticeX[i][j][modulo(k+1,L)] - LatticeY[i][j][k]*LatticeY[i][j][modulo(k+1,L)];
            }        
        }
    }
}

inline void Magnetization(){
    double sx=0, sy=0;
    for(int i = 0; i<L; i++){
        for(int j = 0; j<L; j++){
            for(int k =0; k<L; k++){
                sx += LatticeX[i][j][k];
                sy += LatticeY[i][j][k];
            }
        }
    }
    M = sqrt(sx*sx + sy*sy);
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
            for(int k = 0; k<L; k++){
                double Sdot = LatticeX[i][j][k]*rx + LatticeY[i][j][k]*ry; //Dot product of S with r
                LatticeRX[i][j][k] = LatticeX[i][j][k] - 2*rx*Sdot; //Reflected spin.
                LatticeRY[i][j][k] = LatticeY[i][j][k] - 2*ry*Sdot;
            
                //---Creating bond with the right neigbour---//
                double SNdot = rx*LatticeX[i][modulo(j+1,L)][k] + ry*LatticeY[i][modulo(j+1,L)][k]; //Dot product of the right neighbour with r
                p = 0;
                if (SNdot*Sdot >=0){ p = 1-exp(-2*beta*SNdot*Sdot);}
                double R = ((double) rand() / (RAND_MAX));
                if (R<p){xBonds[i][j][k] = 1;}

                //---Creating bond with the lower neigbour---//
                SNdot = rx*LatticeX[modulo(i+1,L)][j][k] + ry*LatticeY[modulo(i+1,L)][j][k];
                p = 0;
                if (SNdot*Sdot >=0){p = 1-exp(-2*beta*SNdot*Sdot);}
                R = ((double) rand() / (RAND_MAX));
                if (R<p){yBonds[i][j][k] = 1;}
                //---Creating bond with the front neigbour---//
                SNdot = rx*LatticeX[i][j][modulo(k+1,L)] + ry*LatticeY[i][j][modulo(k+1,L)];
                p = 0;
                if (SNdot*Sdot >=0){p = 1-exp(-2*beta*SNdot*Sdot);}
                R = ((double) rand() / (RAND_MAX));
                if (R<p){zBonds[i][j][k] = 1;}
            }
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
            for(int k = 0; k<L; k++){
                std::vector<int> LBond(6), UBond(6), TBond(6), pLabels;
                int bonds = 0, pLabel;
                if (i>0 && yBonds[i-1][j][k] == 1){
                    UBond[bonds] = i - 1;
                    LBond[bonds] = j;
                    TBond[bonds] = k;
                    bonds += 1;
                    pLabel = find(Labels[UBond[bonds-1]][LBond[bonds-1]][TBond[bonds-1]]);
                    pLabels.push_back(pLabel);
                }
                if (i == L-1 && yBonds[i][j][k] == 1){
                    UBond[bonds] = 0;
                    LBond[bonds] = j;
                    TBond[bonds] = k;
                    bonds += 1;
                    pLabel = find(Labels[UBond[bonds-1]][LBond[bonds-1]][TBond[bonds-1]]);
                    pLabels.push_back(pLabel);
                }
                if (j >0 && xBonds[i][j-1][k] == 1){
                    UBond[bonds] = i;
                    LBond[bonds] = j-1;
                    TBond[bonds] = k;
                    bonds += 1;
                    pLabel = find(Labels[UBond[bonds-1]][LBond[bonds-1]][TBond[bonds-1]]);
                    pLabels.push_back(pLabel);
                }
                if (j == L-1 && xBonds[i][j][k] == 1){
                    UBond[bonds] = i;
                    LBond[bonds] == 0;
                    TBond[bonds] = k;
                    bonds += 1;
                    pLabel = find(Labels[UBond[bonds-1]][LBond[bonds-1]][TBond[bonds-1]]);
                    pLabels.push_back(pLabel);
                }
                if (k>0 && zBonds[i][j][k-1] == 1){
                    UBond[bonds] = i;
                    LBond[bonds] = j;
                    TBond[bonds] = k-1;
                    bonds += 1;
                    pLabel = find(Labels[UBond[bonds-1]][LBond[bonds-1]][TBond[bonds-1]]);
                    pLabels.push_back(pLabel);
                }
                if (k == L-1 && zBonds[i][j][k] == 1){
                    UBond[bonds] = i;
                    LBond[bonds] = j;
                    TBond[bonds] = 0;
                    bonds += 1;
                    pLabel = find(Labels[UBond[bonds-1]][LBond[bonds-1]][TBond[bonds-1]]);
                    pLabels.push_back(pLabel);
                }

                if (bonds == 0){
                    Labels[i][j][k] = label;
                    CLabels[label] = label;
                    label += 1;
                }
                else{
                    int minLabel = *std::min_element(pLabels.begin(), pLabels.end());
                    Labels[i][j][k] = minLabel;
                    for (int b = 0; b<bonds; b++){
                        pLabel = pLabels[b];
                        CLabels[pLabel] = minLabel;
                        Labels[UBond[b]][LBond[b]][TBond[b]] = minLabel;
                    }
                }
            }
        }
    }
    for(int i = 0; i<L; i++){
        for(int j = 0; j<L; j++){
            for(int k = 0; k<L; k++){
                Labels[i][j][k] = find(Labels[i][j][k]);
            }
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
            for(int k = 0; k<L; k++){                   
                if (probs[Labels[i][j][k]] < 0.5){
                    LatticeX[i][j][k] = LatticeRX[i][j][k];
                    LatticeY[i][j][k] = LatticeRY[i][j][k];
                }
            }
        }
    }

}

inline void reset(){
    for(int i = 0; i<L; i++){
        for(int j=0; j<L; j++){
            for(int k = 0; k<L; k++){
                xBonds[i][j][k] = 0; yBonds[i][j][k] = 0; zBonds[i][j][k] = 0; Labels[i][j][k] = 0; LatticeRX[i][j][k] = 0; LatticeRY[i][j][k] = 0;
            } 
        }
    }
    for(int i = 0; i<label; i++){
        CLabels[i] = 0;
    }
    label = 0;
}    

void CA_XY3d(double beta, int Ntherm, int Nmeas, int Nsteps){
    std::vector<double> Energy(Nmeas), Energy2(Nmeas), Magn(Nmeas), Magn2(Nmeas);
    //Thermalization//
    initialize_lattice();
    for(int i = 0; i<Ntherm; i++){
        Bonds(beta); //Computes the bonds and the reflected lattice.
        HoshenKopelman(); //Identifies the clusters.
        flip(); //Flips the spins with probability 1/2.
        reset();//Resets bonds and labels.
    }
    for(int i = 0; i<Nmeas; i++){
        Bonds(beta); //Computes the bonds and the reflected lattice.
        HoshenKopelman(); //Identifies the clusters.
        flip(); //Flips the spins with probability 1/2.
        reset();//Resets bonds and labels.
        
        EnergyCA(); //We compute the energy
        Magnetization(); //Magnetization
        Energy[i] = E; 
        Energy2[i] = E*E;
        Magn[i] = M;
        Magn2[i] = M*M;
        E = 0; M = 0;
        for(int j = 0; j<Nsteps; j++){
            Bonds(beta); //Computes the bonds and the reflected lattice.
            HoshenKopelman(); //Identifies the clusters.
            flip(); //Flips the spins with probability 1/2.
            reset(); //Resets bonds and labels.
        }
    }
    E = mean(Energy); dE = Jackknife_error(Energy, 20);
    En2 = mean(Energy2); dE2 = Jackknife_error(Energy2, 20);
    Cv = beta * beta * (En2-E*E) /(L*L*L); dCv = absVal(beta * beta * (dE2 + 2*E*dE)/(L*L*L));
    M = mean(Magn); dM = Jackknife_error(Magn, 20);
    M2 = mean(Magn2); dM2 = Jackknife_error(Magn2, 20);
    chi = (M2-M*M) /(L*L*L); dchi = absVal((dM2 + 2*M*dM)/(L*L*L));
}

//-----------------------------------------//
int main(){    
srand(time(0));
    int Ntherm, Nmeas, Nsteps, Nbeta;
    double beta_min, beta_max; 
    //---Input data---//
    std::cout << "L " << L << std::endl;
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
        if (Betas[i] < 1){sprintf(NameData, "3DXY_L%d_Meas%d_b0%d.dat", L, Nmeas,beta);}
        else{beta = beta + (int) Betas[i] * 10000; sprintf(NameData, "3DXY_L%d_Meas%d_b%d.dat", L, Nmeas,beta);}
        
        std::ofstream Datfile;
        Datfile.open(NameData);
        clock_t begin = clock();
        std::cout << "beta = " << Betas[i] << "  T = " << 1/Betas[i] << std::endl;
        CA_XY3d(Betas[i], Ntherm, Nmeas, Nsteps);
        sprintf(Data_str,"%-30d%-30d%-30d%-30d%-30.17g\n",Ntherm,Nmeas,Nsteps,L,Betas[i]);
        Datfile << Data_str;
        sprintf(Data_str, "%-30.17g%-30.17g\n", E, dE);
        Datfile << Data_str;
        sprintf(Data_str, "%-30.17g%-30.17g\n", Cv, dCv);
        Datfile << Data_str;
        sprintf(Data_str, "%-30.17g%-30.17g\n", M, dM);
        Datfile << Data_str;
        sprintf(Data_str, "%-30.17g%-30.17g\n", chi, dchi);
        Datfile << Data_str;
        
        std::cout << "E = " << E << " +- " << dE << std::endl;
        std::cout << "Cv = " << Cv << " +- " << dCv << std::endl;
        std::cout << "M = " << M << " +- " << dM << std::endl;
        std::cout << "Chi = " << chi << " +- " << dchi << std::endl; 
        E = 0; dE= 0; M = 0; dM = 0;
        //----Computing time----//
        clock_t end = clock();
        double elapsed_secs = double(end - begin) / CLOCKS_PER_SEC;
        sprintf(Data_str, "%-30.17g", elapsed_secs);
        Datfile << Data_str; 
        std::cout << "Time = " << elapsed_secs << " s" << std::endl;
        std::cout << "------------------------------" << std::endl;
        Datfile.close();
    }
    
    return 0;       
}
