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
double Vort=0, Antvort=0;
double dVort, dAntvort;
double p;//Acceptance ratio
int label = 0;
double phi, beta; //Spin angle

constexpr int L = 8;
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

inline void random_init(){
    for(int i = 0; i<L; i++){
        for(int j = 0; j<L; j++){
            for(int k = 0; k<L; k++){
                    double r = rand_range(0, 2*pi);
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
    //Vorticity of a plaquette:1/(2pi) * ( delta(i,i+x) + delta(i+x,i+x+y) + delta(i+x+y, i+y) + delta(i+y, i)
    //delta(i,j) = theta(j) - theta(i) )
    for(int i = 0; i<L; i++){
        for(int j = 0; j<L; j++){
            for (int k = 0; k<L; k++){
                angle(LatticeX[i][j][k], LatticeY[i][j][k]); double theta_i = phi;
                angle(LatticeX[i][modulo(j+1,L)][k], LatticeY[i][modulo(j+1,L)][k]); double theta_ix = phi;
                angle(LatticeX[modulo(i-1,L)][modulo(j+1,L)][k], LatticeY[modulo(i-1,L)][modulo(j+1,L)][k]); double theta_ixy = phi;
                angle(LatticeX[modulo(i-1,L)][j][k], LatticeY[modulo(i-1,L)][j][k]); double theta_iy = phi;
                angle(LatticeX[i][j][modulo(k-1,L)], LatticeY[i][j][modulo(k-1,L)]); double theta_iz = phi;
                angle(LatticeX[modulo(i-1,L)][j][modulo(k-1,L)], LatticeY[modulo(i-1,L)][j][modulo(k-1,L)]); double theta_izy = phi;
                angle(LatticeX[i][modulo(j+1,L)][modulo(k-1,L)], LatticeY[i][modulo(j+1,L)][modulo(k-1,L)]); double theta_ixz = phi;

                //WE GO AROUND THE PLAQUETTES IN CLOCWKISE DIRECTION//
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


//-----Function that generates bonds between neighbouring sites-----//
inline void Bonds(){
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

int T_to_String(double Tf){
    int TF = (Tf - (int) Tf)*1000;
    if (Tf>=1){TF = TF + (int) Tf*1000;}
    return TF;
}

inline void Cluster_XY3d(int tauQ, int Ntherm, int Nmeas, double Ti,double Tf){
    // TauQ --> inverse cooling rating
    // Ntherm --> Thermalization steps

    //double 2*Tc = 4.40334654...; 
    int tmax, MedSize;
    if (Tf == 0){tmax = 2*tauQ - 1; MedSize = 2*tauQ;}
    else {tmax = 2*tauQ; MedSize = 2*tauQ + 1;}
    std::vector<std::vector<double>> Mediciones(MedSize, std::vector<double>(Nmeas,0));
    for(int k = 0; k<Nmeas; k++){
        beta = 1/Ti;
        //Random initial condition//
        random_init();
        //Thermalization//
        for(int i = 0; i<Ntherm; i++){
            Bonds(); 
            HoshenKopelman(); 
            flip(); 
            reset();
        }
        //Cooling// 
        for(int i = 0; i<=tmax; i++){
            beta = 1.0/ ( (Tf-Ti)/(2.0*tauQ) * i + Ti  );
            Bonds(); 
            HoshenKopelman(); 
            flip(); 
            reset();

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
    std::cout << "-----Cluster algorithm 3DXY-----" << std::endl;
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
    Cluster_XY3d(tauQ_min, Ntherm, Nmeas, Ti,Tf); 
    clock_t end = clock();
    double elapsed_secs = double(end - begin) / CLOCKS_PER_SEC; 
    std::cout << "Time = " << elapsed_secs << " s" << std::endl;
    
    return 0;
}