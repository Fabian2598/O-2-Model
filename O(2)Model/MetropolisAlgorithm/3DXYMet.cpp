#include <iostream>
#include <vector> 
#include <random> 
#include <time.h> 
#include <fstream>
#include <cmath> 
#include "statistics.h"


double pi = 3.14159265359;
double E = 0, Vort = 0, Antvort = 0;
double dE = 0, dVort = 0, dAntvort = 0;
double E_av=0, V_av=0;

const int L = 16;
long double SpinLattice[L][L][L];

//-----Computes the energy of a configuration-----//
void random_init(){
    for(int i = 0; i<L; i++){
        for(int j = 0; j<L; j++){
            for(int k = 0; k<L; k++){
                    SpinLattice[i][j][k] = rand_range(0.0, 2*pi);
            }   
        }
    }
}

inline void EnergyMeas(){
    for (int i = 0; i<L; i++){
        for (int j = 0; j<L; j++){
            for(int k = 0; k<L; k++){   
                E += -cos(SpinLattice[i][j][k] - SpinLattice[i][modulo(j+1,L)][k]) 
                     -cos(SpinLattice[i][j][k] - SpinLattice[modulo(i+1,L)][j][k])
                     -cos(SpinLattice[i][j][k] - SpinLattice[i][j][modulo(k+1,L)]);
            }
        }
    }
}

inline double correction(double a){    
    if (a>pi){
        return (-2*pi + a);
    }
    return a;
}

inline void Vorticity(){
    //Vorticity of a plaquette delta(i,i+x) + delta(i+x,i+x+y) + delta(i+x+y, i+y) + delta(i+y, i)
    //delta(i,j) = theta(j) - theta(i)
    for(int i = 0; i<L; i++){
        for(int j = 0; j<L; j++){
            for(int k = 0; k<L; k++){
                double theta_i = SpinLattice[i][j][k];
                double theta_ix = SpinLattice[i][modulo(j+1,L)][k];
                double theta_ixy = SpinLattice[modulo(i-1,L)][modulo(j+1,L)][k];
                double theta_iy = SpinLattice[modulo(i-1,L)][j][k];
                double theta_iz = SpinLattice[i][j][modulo(k-1,L)];
                double theta_izy = SpinLattice[modulo(i-1,L)][j][modulo(k-1,L)];
                double theta_ixz = SpinLattice[i][modulo(j+1,L)][modulo(k-1,L)];

                //WE GO AROUND THE PLAQUETTES IN CLOCKWISE DIRECTION//
                //Plaquette 1 
                double q = 1/(2*pi) * ( correction(fmodulo(-theta_ix+theta_i,2*pi)) + 
                correction(fmodulo(-theta_ixy+theta_ix,2*pi)) + 
                correction(fmodulo(-theta_iy+theta_ixy,2*pi)) +
                correction(fmodulo(-theta_i+theta_iy,2*pi)) );
                if ( absVal(q-1) <= 1e-6) {Vort += 1;}
                else if ( absVal(q+1) <= 1e-6){Antvort += 1;}

                //Plaquette 2
                q = 1/(2*pi) * ( correction(fmodulo(theta_iz-theta_i,2*pi)) + 
                correction(fmodulo(theta_izy-theta_iz,2*pi)) + 
                correction(fmodulo(theta_iy-theta_izy,2*pi)) +
                correction(fmodulo(theta_i-theta_iy,2*pi)) );
                if ( absVal(q-1) <= 1e-6) {Vort += 1;}
                else if ( absVal(q+1) <= 1e-6){Antvort += 1;}

                //Plaquette 3
                q = 1/(2*pi) * ( correction(fmodulo(theta_ix-theta_i,2*pi)) + 
                correction(fmodulo(theta_ixz-theta_ix,2*pi)) + 
                correction(fmodulo(theta_iz-theta_ixz,2*pi)) +
                correction(fmodulo(theta_i-theta_iz,2*pi)) );
                if ( absVal(q-1) <= 1e-6) {Vort += 1;}
                else if ( absVal(q+1) <= 1e-6){Antvort += 1;}
            }
        }
    }  
}
    
//-----Choose a new angle-----//
inline void Reflection(int i, int j, int k,double beta){
    double phi = SpinLattice[i][j][k];
    //std::cout << "i " << i << " j " << j << " k " << k << " phi " << phi << std::endl; 
    double R = rand_range(0.0, 1.0);
    double theta = rand_range(0.0, 2*pi);//fmodulo(rand_range(0,2*pi)+phi,2*pi);
    double deltaE = cos(phi-SpinLattice[modulo(i-1,L)][j][k]) + cos(phi-SpinLattice[modulo(i+1,L)][j][k])
            +cos(phi-SpinLattice[i][modulo(j-1,L)][k]) + cos(phi-SpinLattice[i][modulo(j+1,L)][k])
            +cos(phi-SpinLattice[i][j][modulo(k-1,L)]) + cos(phi-SpinLattice[i][j][modulo(k+1,L)])
            -cos(theta-SpinLattice[modulo(i-1,L)][j][k]) - cos(theta-SpinLattice[modulo(i+1,L)][j][k])
            -cos(theta-SpinLattice[i][modulo(j-1,L)][k]) - cos(theta-SpinLattice[i][modulo(j+1,L)][k])
            -cos(theta-SpinLattice[i][j][modulo(k-1,L)]) - cos(theta-SpinLattice[i][j][modulo(k+1,L)]);
    //std::cout << " delta E " << deltaE << std::endl;
    double p = exp(-beta*deltaE);
    if (R<p){
        SpinLattice[i][j][k] = theta;
        //std::cout << "accepted" << std::endl;
    }     
}


inline void Met_XY3d(double beta,int Ntherm, int Nmeas, int Nstep){
    std::vector<double> Energy(Nmeas), VORT(Nmeas), AVORT(Nmeas);
    //Random initial condition
    random_init();
    //Thermalization
    for(int i = 0; i<Ntherm*L*L*L; i++){
            Reflection(rand() % L,rand() % L,rand() % L ,beta);  
    }
    //Measurements//
    for(int i = 0; i<Nmeas; i++){
        for(int j = 0; j<L*L*L; j++){
            Reflection(rand() % L,rand() % L,rand() % L ,beta);  
        }
        EnergyMeas(); 
        Vorticity(); 
        Energy[i] = E; 
        VORT[i] = Vort;
        AVORT[i] = Antvort;
        E_av += E;
        V_av += Vort;
        E=0; Vort=0; Antvort = 0;
        //Steps between measurements//
        for(int j = 0; j<Nstep*L*L*L; j++){
            Reflection(rand() % L,rand() % L,rand() % L ,beta);  
        }    
    }
    std::cout << " L " << L << std::endl;
    E_av /= Nmeas;
    V_av /= (Nmeas*L*L*L);
    E = mean(Energy); dE = Jackknife(Energy, {5,10,20,50,100});
    Vort = mean(VORT)/(L*L*L); dVort = Jackknife(VORT,{5,10,20,50,100})/(L*L*L);
    Antvort = mean(AVORT)/(L*L*L); dAntvort = Jackknife(AVORT,{5,10,20,50,100})/(L*L*L);
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
        Met_XY3d(Betas[i], Ntherm, Nmeas, Nsteps);
        sprintf(Data_str,"%-30d%-30d%-30d%-30d%-30.17g\n",Ntherm,Nmeas,Nsteps,L,Betas[i]);
        Datfile << Data_str;
        sprintf(Data_str, "%-30.17g%-30.17g\n", E, dE);
        Datfile << Data_str;
        sprintf(Data_str, "%-30.17g%-30.17g\n", Vort, dVort);
        Datfile << Data_str;
        sprintf(Data_str, "%-30.17g%-30.17g\n", Antvort, dAntvort);
        Datfile << Data_str;

        std::cout << "E = " << E << " +- " << dE << std::endl;
        std::cout << "Eav = " << E_av << std::endl;
        std::cout << "V = " << Vort << " +- " << dVort << std::endl;
        std::cout << "Vav = " << V_av << std::endl;
        std::cout << "A = " << Antvort << " +- " << dAntvort << std::endl;
        E = 0; dE= 0; Vort = 0; dVort=0; Antvort = 0; dAntvort = 0;
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
