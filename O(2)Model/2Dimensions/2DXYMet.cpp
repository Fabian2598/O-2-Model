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
double Vort = 0, Antvort = 0;
double dVort, dAntvort;
double E = 0, En2 = 0, M = 0, M2=0, chi = 0, Cv = 0;
double dchi, dM, dM2, dCv, dE, dE2;
double p;

constexpr int L = 8;
constexpr  int maxsize = L*L;
std::vector<std::vector<double>> SpinLattice(L, std::vector<double>(L,0));
std::vector<std::vector<double>> SpinLatticeR(L, std::vector<double>(L,0));

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
            double q = 1/(2*pi) * ( 
            correction(fmodulo(SpinLattice[i][modulo(j+1,L)]- SpinLattice[i][j],2*pi))   +
            correction(fmodulo(SpinLattice[modulo(i-1,L)][modulo(j+1,L)]- SpinLattice[i][modulo(j+1,L)],2*pi)) +
            correction(fmodulo(SpinLattice[modulo(i-1,L)][j]- SpinLattice[modulo(i-1,L)][modulo(j+1,L)],2*pi))  +
            correction(fmodulo(SpinLattice[i][j]- SpinLattice[modulo(i-1,L)][j],2*pi))
            );
            if ( absVal(q-1) <= 1e-6) {Vort += 1;}
            else if ( absVal(q+1) <= 1e-6){Antvort += 1;}
        }
    }  
}
    
inline void EnergyM(){
    for (int i = 0; i<L; i++){
        for (int j = 0; j<L; j++){
            E += -cos(SpinLattice[i][j]-SpinLattice[i][modulo(j+1,L)]) 
                 - cos(SpinLattice[i][j]-SpinLattice[modulo(i+1,L)][j]);
        }
    }
}

inline void Magnetization(){
    double sx=0, sy=0;
    for(int i = 0; i<L; i++){
        for(int j = 0; j<L; j++){
            sx += cos(SpinLattice[i][j]);
            sy += sin(SpinLattice[i][j]);
        }
    }
    M = sqrt(sx*sx + sy*sy);
}


//-----Chooses a new angle-----//
inline void Reflection(double beta, std::vector<int> site){
    double phi = SpinLattice[site[0]][site[1]];
    double R = ((double) rand() / (RAND_MAX));
    double theta = rand_range(-pi, pi);
    double deltaE = cos(phi-SpinLattice[modulo(site[0]-1,L)][site[1]]) + cos(phi-SpinLattice[modulo(site[0]+1,L)][site[1]])
            +cos(phi-SpinLattice[site[0]][modulo(site[1]-1,L)]) + cos(phi-SpinLattice[site[0]][modulo(site[1]+1,L)])
            -cos(theta-SpinLattice[modulo(site[0]-1,L)][site[1]]) - cos(theta-SpinLattice[modulo(site[0]+1,L)][site[1]])
            -cos(theta-SpinLattice[site[0]][modulo(site[1]-1,L)]) - cos(theta-SpinLattice[site[0]][modulo(site[1]+1,L)]);
    p = 1;
    if (deltaE*beta>0){p = exp(-beta*deltaE);}
    if (R<p){SpinLattice[site[0]][site[1]] = theta;}         
}

inline void Met_XY2d(double beta, int Ntherm,int Nmeas,int Nsteps){
    std::vector<double> VORT(Nmeas), AVORT(Nmeas);
    std::vector<double> Energy(Nmeas), Energy2(Nmeas), Magn(Nmeas), Magn2(Nmeas);
    //Random initial configuration//
    for(int l = 0; l<L; l++){
        for(int m = 0; m<L; m++){
            SpinLattice[l][m] = rand_range(-pi,pi);
        }
    }
    //Thermalization//
    for(int i = 0; i<Ntherm*L*L; i++){
        Reflection(beta, {rand() % L,rand() % L });  
    }
    for(int i = 0; i<Nmeas; i++){
        //One sweep//
        for(int k = 0; k<L*L; k++){
            Reflection(beta,{rand() % L,rand() % L }); 
        } 
        EnergyM();
        Magnetization();
        Vorticity(); 
        Energy[i] = E; 
        Energy2[i] = E*E;
        Magn[i] = M;
        Magn2[i] = M*M;
        VORT[i] = Vort;
        AVORT[i] = Antvort;
        E = 0; M = 0; Vort = 0; Antvort = 0;
        //Decorrelation sweeps//
        for(int k = 0; k<Nsteps*L*L; k++){
            Reflection(beta,{rand() % L,rand() % L }); 
        }
    }
    E = mean(Energy); dE = Jackknife_error(Energy, 20);
    En2 = mean(Energy2); dE2 = Jackknife_error(Energy2, 20);
    Cv = beta * beta * (En2-E*E) /(L*L); dCv = absVal(beta * beta * (dE2 + 2*E*dE)/(L*L));
    M = mean(Magn); dM = Jackknife_error(Magn, 20);
    M2 = mean(Magn2); dM2 = Jackknife_error(Magn2, 20);
    chi = (M2-M*M) /(L*L); dchi = absVal((dM2 - 2*M*dM)/(L*L));
    Vort = mean(VORT)/(L*L); dVort = Jackknife_error(VORT,20)/(L*L);
    Antvort = mean(AVORT)/(L*L); dAntvort = Jackknife_error(AVORT,20)/(L*L);
}



//-----------------------------------------//
int main(){    
srand(time(0));
    int Ntherm, Nmeas, Nsteps, Nbeta;
    double beta_min, beta_max; 
    //---Input data---//
    std::cout << "Metropolis" << std::endl;
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
    char NameData[1000], Data_str[1000];
    for (int i = 0; i < Nbeta; i++) {
        int A = Betas[i];
        int beta= (Betas[i] - A) * 10000;
        if (Betas[i] < 1){sprintf(NameData, "2DXY_L%d_Meas%d_b0%d.dat", L, Nmeas,beta);}
        else{beta = beta + (int) Betas[i] * 10000; sprintf(NameData, "2DXY_L%d_Meas%d_b%d.dat", L, Nmeas,beta);}
        
        std::ofstream Datfile;
        Datfile.open(NameData);
        clock_t begin = clock();
        std::cout << "beta = " << Betas[i] << "  T = " << 1/Betas[i] << std::endl;
        Met_XY2d(Betas[i], Ntherm, Nmeas, Nsteps);
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
        sprintf(Data_str, "%-30.17g%-30.17g\n", Vort, dVort);
        Datfile << Data_str;
        sprintf(Data_str, "%-30.17g%-30.17g\n", Antvort, dAntvort);
        Datfile << Data_str;
        
        std::cout << "E = " << E << " +- " << dE << std::endl;
        std::cout << "Cv = " << Cv << " +- " << dCv << std::endl;
        std::cout << "M = " << M << " +- " << dM << std::endl;
        std::cout << "Chi = " << chi << " +- " << dchi << std::endl; 
        std::cout << "V = " << Vort << " +- " << dVort << std::endl;
        std::cout << "A = " << Antvort << " +- " << dAntvort << std::endl;
        E = 0; dE= 0; M = 0; dM = 0; Vort = 0; dVort=0; Antvort = 0; dAntvort = 0;
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
