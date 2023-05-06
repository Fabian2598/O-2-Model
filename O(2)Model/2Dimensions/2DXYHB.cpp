#include <iostream>
#include <vector> 
#include <random> 
#include <time.h> 
#include <fstream>
#include <cmath> 
#include "statistics.h"


double pi = 3.14159265359;
double Vort = 0, Antvort = 0;
double dVort, dAntvort;
double E = 0, En2 = 0, M = 0, M2=0, chi = 0, Cv = 0;
double dchi, dM, dM2, dCv, dE, dE2;
double delta = 0, alpha = 0, beta_a = 0, h_a = 0, G_a = 0, g_a = 0;
double a_const = 0.798953686083986, eps = 0.001, aVal = 0;
double tx;

const int L = 16;
double SpinLattice[L][L];

//-----Computes the energy of a configuration-----//
void random_init(){
    for(int i = 0; i<L; i++){
        for(int j = 0; j<L; j++){
                SpinLattice[i][j] = rand_range(-pi, pi);
              
        }
    }
}

inline void EnergyMeas(){
    for (int i = 0; i<L; i++){
        for (int j = 0; j<L; j++){  
                E += -cos(SpinLattice[i][j] - SpinLattice[i][modulo(j+1,L)]) 
                     -cos(SpinLattice[i][j] - SpinLattice[modulo(i+1,L)][j]);
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
    
//delta function//    
inline void delta_function(){
    std::vector<double> v1{0,aVal-a_const};
    double maxval = *std::max_element(v1.begin(), v1.end());
    delta = 0.35*maxval + 1.03*sqrt(maxval);
}

//alpha function//
inline void alpha_function(){
    delta_function(); //We compute delta(a).
    std::vector<double> v1{sqrt(eps*aVal),delta};
    double maxval = *std::max_element(v1.begin(), v1.end());
    std::vector<double> v2{sqrt(aVal*(2-eps)),maxval};
    alpha = *std::min_element(v2.begin(), v2.end());
}

//beta function//
inline void beta_function(){
    std::vector<double> v1{alpha*alpha/aVal, ( cosh(pi*alpha)-1 )/(exp(2*aVal)-1) };
    beta_a = *std::max_element(v1.begin(), v1.end())-1.0;
}

//h function//
inline void h_function(double x){
    alpha_function(); beta_function();
    h_a = 2/alpha * atanh( sqrt((1+beta_a)/(1-beta_a)) * tan(    (2*x-1) * atan( tanh( pi*alpha/2 )  * sqrt((1-beta_a)/(1+beta_a)) ) ) );
}

//G function//
inline void G_function(double x){
    G_a = 1.0 - cos(x) - 1.0/aVal * log(1.0 + (cosh(alpha*x)-1.0)/(1.0+beta_a) );
}

//g function//
inline void g_function(double x){
    h_function(x);
    G_function(h_a);
    g_a = exp(-aVal*G_a);
}

inline void angle(double x, double y){
    if (x>1e-6){ tx = atan(y/x); }
    else if(x<-1e-6 && y>1e-6){ tx = pi + atan(y/x); }
    else if(x<-1e-6 && y<-1e-6){ tx = -pi + atan(y/x);}

    else if (absVal(x)<1e-6){ tx = pi/2 * y/absVal(y); }
    else if (absVal(y)<1e-6){ tx = pi * (1-x/absVal(x))/2; } 
}

inline void HB_sweep(double beta){
    double Sx, Sy;
    double w1, w2;
    for(int i = 0; i<L; i++){
        for(int j = 0; j<L; j++){
            for(int k = 0; k<L; k++){
                Sx = cos(SpinLattice[modulo(i+1,L)][j]) + cos(SpinLattice[modulo(i-1,L)][j]) 
                + cos(SpinLattice[i][modulo(j+1,L)]) + cos(SpinLattice[i][modulo(j-1,L)]); 
                Sy = sin(SpinLattice[modulo(i+1,L)][j]) + sin(SpinLattice[modulo(i-1,L)][j]) 
                + sin(SpinLattice[i][modulo(j+1,L)]) + sin(SpinLattice[i][modulo(j-1,L)]);
                aVal = beta*sqrt(Sx*Sx + Sy*Sy);
                angle(Sx,Sy); //Now we have tx from -pi to pi 
                bool boolean = true;
                while (boolean == true){
                    w1 = rand_range(0.0,1.0); w2 = rand_range(0.0,1.0);
                    g_function(w1);
                    if (w2 <= g_a){
                        angle( cos(h_a+tx),sin(h_a+tx)  );
                        SpinLattice[i][j] = tx;
                        boolean = false;
                    }
                }
                Sx = 0.0; Sy = 0.0;      
            }
        }
    }
}

inline void HB_XY2d(double beta,int Ntherm, int Nmeas, int Nstep){
    std::vector<double> VORT(Nmeas), AVORT(Nmeas);
    std::vector<double> Energy(Nmeas), Energy2(Nmeas), Magn(Nmeas), Magn2(Nmeas);
    //Random initial condition
    random_init();
    //Thermalization
    for(int i = 0; i<Ntherm; i++){
            HB_sweep(beta);  
    }
    //Measurements//
    for(int i = 0; i<Nmeas; i++){
        HB_sweep(beta);
        EnergyMeas(); 
        Magnetization();
        Vorticity(); 
        Energy[i] = E; 
        Energy2[i] = E*E;
        Magn[i] = M;
        Magn2[i] = M*M;
        VORT[i] = Vort;
        AVORT[i] = Antvort;
        E = 0; M = 0; Vort = 0; Antvort = 0;
        //Steps between measurements//
        for(int j = 0; j<Nstep; j++){
            HB_sweep(beta);
        }    
    }
    
    E = mean(Energy); dE = Jackknife(Energy, {5,10,20,50});
    En2 = mean(Energy2); dE2 = Jackknife(Energy2, {5,10,20,50});
    Cv = beta * beta * (En2-E*E) /(L*L); dCv = absVal(beta * beta * (dE2 + 2*E*dE)/(L*L));

    M = mean(Magn); dM = Jackknife(Magn, {5,10,20,50});
    M2 = mean(Magn2); dM2 = Jackknife(Magn2, {5,10,20,50});
    chi = (M2-M*M) /(L*L); dchi = absVal((dM2 - 2*M*dM)/(L*L));

    Vort = mean(VORT)/(L*L); dVort = Jackknife(VORT,{5,10,20,50})/(L*L);
    Antvort = mean(AVORT)/(L*L); dAntvort = Jackknife(AVORT,{5,10,20,50})/(L*L);
}

//-----------------------------------------//
int main(){   
    srand(time(0));
    int Ntherm, Nmeas, Nsteps, Nbeta;
    double beta_min, beta_max; 
    //---Input data---//
    std::cout << "Heat Bath. L=" << L << std::endl;
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
        HB_XY2d(Betas[i], Ntherm, Nmeas, Nsteps);
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
