#include <iostream>
#include <vector> 
#include <random> 
#include <time.h> 
#include <fstream>
#include <cmath> 
#include "statistics.h"


double pi = 3.14159265359;
double E = 0, Vort = 0, Antvort = 0;
double dE, dVort, dAntvort;

constexpr int L = 8;
constexpr  int maxsize = L*L*L;
double SpinLattice[L][L][L];

inline void random_init(){
    for(int i = 0; i<L; i++){
        for(int j = 0; j<L; j++){
            for(int k = 0; k<L; k++){
                    SpinLattice[i][j][k] = rand_range(0.0, 2*pi);
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
    double R = rand_range(0.0, 1.0);
    double theta = rand_range(0.0, 2*pi);
    double deltaE = cos(phi-SpinLattice[modulo(i-1,L)][j][k]) + cos(phi-SpinLattice[modulo(i+1,L)][j][k])
            +cos(phi-SpinLattice[i][modulo(j-1,L)][k]) + cos(phi-SpinLattice[i][modulo(j+1,L)][k])
            +cos(phi-SpinLattice[i][j][modulo(k-1,L)]) + cos(phi-SpinLattice[i][j][modulo(k+1,L)])
            -cos(theta-SpinLattice[modulo(i-1,L)][j][k]) - cos(theta-SpinLattice[modulo(i+1,L)][j][k])
            -cos(theta-SpinLattice[i][modulo(j-1,L)][k]) - cos(theta-SpinLattice[i][modulo(j+1,L)][k])
            -cos(theta-SpinLattice[i][j][modulo(k-1,L)]) - cos(theta-SpinLattice[i][j][modulo(k+1,L)]);
    double p = exp(-beta*deltaE);
    if (R<p){
        SpinLattice[i][j][k] = theta;
    }     
}

int T_to_String(double Tf){
    int TF = (Tf - (int) Tf)*1000;
    if (Tf>=1){TF = TF + (int) Tf*1000;}
    return TF;
}

inline void Met_XY3dV2(int tauQ, int Ntherm, int Nmeas, double Ti,double Tf){
    // TauQ --> inverse cooling rating
    // Ntherm --> Thermalization steps
    std::vector<double> VORT, AVORT;
    //double Tc = 0.89; 
    int tmax, MedSize;
    if (Tf == 0){tmax = 2*tauQ - 1; MedSize = 2*tauQ;}
    else {tmax = 2*tauQ; MedSize = 2*tauQ + 1;}
    std::vector<std::vector<double>> Mediciones(MedSize, std::vector<double>(Nmeas,0));
    for(int k = 0; k<Nmeas; k++){
        double beta = 1/(Ti);
        //Random initial condition//
        random_init();
        //Thermalization//
        for(int i = 0; i<Ntherm*L*L*L; i++){
            Reflection(rand() % L,rand() % L,rand() % L ,beta);  
        }
        //Cooling// 
        for(int i = 0; i<=tmax; i++){
            beta = 1.0/ ( (Tf-Ti)/(2.0*tauQ) * i + Ti  );
            for(int i = 0; i<L*L*L; i++){
                Reflection(rand() % L,rand() % L,rand() % L ,beta);  
            }
                Vorticity(); 
                Mediciones[i][k] = Vort;
                Vort=0; Antvort = 0;
        }     
    }
    char NameData[50], Data_str[100];
    sprintf(NameData, "3DXY_L%d_tdep_tauQ%d_Ti%d_Tf%d.txt", L, tauQ, T_to_String(Ti), T_to_String(Tf));
    std::ofstream Datfile;
    Datfile.open(NameData);
    for(int i=0;i<=tmax;i++ ){
        Vort = mean(Mediciones[i])/(L*L*L); dVort = Jackknife_error(Mediciones[i],20)/(L*L*L);
        sprintf(Data_str, "%-30d%-30.17g%-30.17g\n",i, Vort, dVort);
        Datfile << Data_str;
    }
    Datfile.close();
}

//-----------------------------------------//
int main(){    
    srand(time(0));
    int Ntherm, Nmeas, tauQ_min;
    double Ti, Tf; 
    //---Input data---//
    std::cout << "-----Metropolis 3DXY-----" << std::endl;
    std::cout << "L = " << L << std::endl;
    std::cout << "tauQ: "; 
    std::cin >> tauQ_min;
    std::cout << "Ti: ";
    std::cin >> Ti;
    std::cout << "Tf: ";
    std::cin >> Tf;
    std::cout << "Thermalization: ";
    std::cin >> Ntherm;
    std::cout << "Measurements: ";
    std::cin >> Nmeas;
    std::cout << " " << std::endl;

    clock_t begin = clock();
    Met_XY3dV2(tauQ_min, Ntherm, Nmeas, Ti,Tf);
    clock_t end = clock();
    double elapsed_secs = double(end - begin) / CLOCKS_PER_SEC; 
    std::cout << "Time = " << elapsed_secs << " s" << std::endl;
    
 
    return 0;
}
