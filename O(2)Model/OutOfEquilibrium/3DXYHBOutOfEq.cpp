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
double delta = 0, alpha = 0, beta_a = 0, h_a = 0, G_a = 0, g_a = 0;
double a_const = 0.798953686083986, eps = 0.001, aVal = 0;
double tx;

constexpr int L = 8;
constexpr  int maxsize = L*L*L;
double SpinLattice[L][L][L];

//Arrays needed for cluster algorithm//
int CLabels[maxsize]; 
double phi = 0;
double p;//Acceptance ratio
int label = 0;

double LatticeX[L][L][L]; 
double LatticeY[L][L][L]; 
double LatticeRX[L][L][L]; 
double LatticeRY[L][L][L]; 

double xBonds[L][L][L]; 
double yBonds[L][L][L]; 
double zBonds[L][L][L]; 
double Labels[L][L][L]; 
//----------------------------------//

inline void random_init(){
    for(int i = 0; i<L; i++){
        for(int j = 0; j<L; j++){
            for(int k = 0; k<L; k++){
                    //SpinLattice[i][j][k] = rand_range(0.0, 2*pi);
                    double r = rand_range(0.0, 2*pi);
                    LatticeX[i][j][k] = cos(r);
                    LatticeY[i][j][k] = sin(r);
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

//--------------Cluster algorithm--------------//
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


//-----Function that generates bonds between neighbouring sites-----//
inline void Bonds(double beta){
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

//---------------------------------------------//

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

inline void angleHB(double x, double y){
    if (x>1e-6){ tx = atan(y/x); }
    else if(x<-1e-6 && y>1e-6){ tx = pi + atan(y/x); }
    else if(x<-1e-6 && y<-1e-6){ tx = -pi + atan(y/x);}

    else if (absVal(x)<1e-6){ tx = pi/2 * y/absVal(y); }
    else if (absVal(y)<1e-6){ tx = pi * (1-x/absVal(x))/2; } 
}

//Heatbath sweep
inline void HB_sweep(double beta){
    double Sx, Sy;
    double w1, w2;
    for(int i = 0; i<L; i++){
        for(int j = 0; j<L; j++){
            for(int k = 0; k<L; k++){
                Sx = cos(SpinLattice[modulo(i+1,L)][j][k]) + cos(SpinLattice[modulo(i-1,L)][j][k]) 
                + cos(SpinLattice[i][modulo(j+1,L)][k]) + cos(SpinLattice[i][modulo(j-1,L)][k]) 
                + cos(SpinLattice[i][j][modulo(k+1,L)]) + cos(SpinLattice[i][j][modulo(k-1,L)]);
                Sy = sin(SpinLattice[modulo(i+1,L)][j][k]) + sin(SpinLattice[modulo(i-1,L)][j][k]) 
                + sin(SpinLattice[i][modulo(j+1,L)][k]) + sin(SpinLattice[i][modulo(j-1,L)][k]) 
                + sin(SpinLattice[i][j][modulo(k+1,L)]) + sin(SpinLattice[i][j][modulo(k-1,L)]);
                aVal = beta*sqrt(Sx*Sx + Sy*Sy);
                angleHB(Sx,Sy); //Now we have tx from -pi to pi 
                bool boolean = true;
                while (boolean == true){
                    w1 = rand_range(0.0,1.0); w2 = rand_range(0.0,1.0);
                    g_function(w1);
                    if (w2 <= g_a){
                        angleHB( cos(h_a+tx),sin(h_a+tx)  );
                        SpinLattice[i][j][k] = tx;
                        boolean = false;
                    }
                }
                Sx = 0.0; Sy = 0.0;      
            }
        }
    }
}


int T_to_String(double Tf){
    int TF = (Tf - (int) Tf)*1000;
    if (Tf>=1){TF = TF + (int) Tf*1000;}
    return TF;
}

inline void HB_XY3dV2(int tauQ, int Ntherm, int Nmeas, double Ti,double Tf){
    // TauQ --> inverse cooling rating
    // Ntherm --> Thermalization steps
    std::vector<double> VORT, AVORT;
    //double Tc = 4.40334654...; 
    int tmax, MedSize;
    if (Tf == 0){tmax = 2*tauQ - 1; MedSize = 2*tauQ;}
    else {tmax = 2*tauQ; MedSize = 2*tauQ + 1;}
    std::vector<std::vector<double>> Mediciones(MedSize, std::vector<double>(Nmeas,0));
    for(int k = 0; k<Nmeas; k++){
        double beta = 1/(Ti);
        //Random initial condition//
        random_init();
        //Thermalization//
        for(int i = 0; i<Ntherm; i++){
            Bonds(beta); 
            HoshenKopelman(); 
            flip(); 
            reset();
        }

        //We transform the spin components of the cluster algorithm to angles for heatbath//
        for(int i = 0; i<L; i++){
            for(int j = 0; j<L; j++){
                for (int k = 0; k<L; k++){
                    angle(LatticeX[i][j][k], LatticeY[i][j][k]); SpinLattice[i][j][k] = phi;
                }
            }
        }
        //--------------------------------------------------------------------------------//

        //Cooling with Heatbath// 
        for(int i = 0; i<=tmax; i++){
            beta = 1.0/ ( (Tf-Ti)/(2.0*tauQ) * i + Ti  );
            HB_sweep(beta);
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
    std::cout << "-----Heatbath 3DXY-----" << std::endl;
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
    HB_XY3dV2(tauQ_min, Ntherm, Nmeas, Ti,Tf);
    clock_t end = clock();
    double elapsed_secs = double(end - begin) / CLOCKS_PER_SEC; 
    std::cout << "Time = " << elapsed_secs << " s" << std::endl;
    
 
    return 0;
}
