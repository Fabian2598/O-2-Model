#include <iostream>
#include <vector> 
#include <random> 
#include <time.h> 
#include <fstream>
#include <ctime>
#include <cmath> 
#include <algorithm>
#include <unordered_set>
#include "statistics.h"


double pi = 3.14159265359;
double Vort, Antvort, strdensity;
double dVort, dAntvort, dstrdensity;
double p;//Acceptance ratio
int label = 0, PlaqLabel = 0;
double phi; //Spin angle

constexpr int L = 8;
constexpr  int maxsize = L*L*L;
int CLabels[maxsize]; //Cluster labels.
int PlaqLabels[maxsize*3]; //Plaquette labels

//We store the spin components.
std::vector<std::vector<std::vector<double>>> LatticeX(L, std::vector<std::vector<double>> (L, std::vector<double>(L,0))); 
std::vector<std::vector<std::vector<double>>> LatticeY(L, std::vector<std::vector<double>> (L, std::vector<double>(L,0))); 
std::vector<std::vector<std::vector<double>>> LatticeRX(L, std::vector<std::vector<double>> (L, std::vector<double>(L,0))); 
std::vector<std::vector<std::vector<double>>> LatticeRY(L, std::vector<std::vector<double>> (L, std::vector<double>(L,0))); 


std::vector<std::vector<std::vector<int>>> xBonds(L, std::vector<std::vector<int> > (L, std::vector<int>(L,0)));
std::vector<std::vector<std::vector<int>>> yBonds(L, std::vector<std::vector<int> > (L, std::vector<int>(L,0)));
std::vector<std::vector<std::vector<int>>> zBonds(L, std::vector<std::vector<int> > (L, std::vector<int>(L,0)));
std::vector<std::vector<std::vector<int>>> Labels(L, std::vector<std::vector<int> > (L, std::vector<int>(L,0)));

//Array of dimension 3 * L*L*L where we store the vorticity of each plaquette
//VortexP[0][i][j][k] --> Plaquette 1 at (i,j,k), VortexP[1][i][j][k] --> Plaquette 2 at (i,j,k), VortexP[2][i][j][k] --> Plaquette 3 at (i,j,k)
//See figures in my notes to understand this.

//MarkerPIn MarkerPOut --> We mark those plaquettes with a string. One end of the string enters the plaquette and the other exits.
//We must mark both sides of the plaquettes, otherwise closed string wouldn't form.
//LabelP --> Label of each plaquette. It tells us to which string belongs each plaquette.
std::vector<std::vector<std::vector<std::vector<int>>>> VortexP(3, std::vector<std::vector<std::vector<int>>>(L, std::vector<std::vector<int>>(L, std::vector<int>(L,0))));
std::vector<std::vector<std::vector<std::vector<int>>>> MarkerPIn(3, std::vector<std::vector<std::vector<int>>>(L, std::vector<std::vector<int>>(L, std::vector<int>(L,0))));
std::vector<std::vector<std::vector<std::vector<int>>>> MarkerPOut(3, std::vector<std::vector<std::vector<int>>>(L, std::vector<std::vector<int>>(L, std::vector<int>(L,0))));
std::vector<std::vector<std::vector<std::vector<int>>>> LabelP(3, std::vector<std::vector<std::vector<int>>>(L, std::vector<std::vector<int>>(L, std::vector<int>(L,0))));



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
                //We associate 3 plaquettes with each site.
                angle(LatticeX[i][j][k], LatticeY[i][j][k]); double theta_i = phi;
                angle(LatticeX[i][modulo(j+1,L)][k], LatticeY[i][modulo(j+1,L)][k]); double theta_ix = phi;
                angle(LatticeX[modulo(i-1,L)][modulo(j+1,L)][k], LatticeY[modulo(i-1,L)][modulo(j+1,L)][k]); double theta_ixy = phi;
                angle(LatticeX[modulo(i-1,L)][j][k], LatticeY[modulo(i-1,L)][j][k]); double theta_iy = phi;
                angle(LatticeX[i][j][modulo(k-1,L)], LatticeY[i][j][modulo(k-1,L)]); double theta_iz = phi;
                angle(LatticeX[modulo(i-1,L)][j][modulo(k-1,L)], LatticeY[modulo(i-1,L)][j][modulo(k-1,L)]); double theta_izy = phi;
                angle(LatticeX[i][modulo(j+1,L)][modulo(k-1,L)], LatticeY[i][modulo(j+1,L)][modulo(k-1,L)]); double theta_ixz = phi;

                //Plaquette 1
                double q = 1/(2*pi) * ( correction(fmodulo(theta_ix-theta_i,2*pi)) + 
                correction(fmodulo(theta_ixy-theta_ix,2*pi)) + 
                correction(fmodulo(theta_iy-theta_ixy,2*pi)) +
                correction(fmodulo(theta_i-theta_iy,2*pi)) );
                if ( absVal(q-1) <= 1e-6) {Vort += 1; VortexP[0][i][j][k] = 1;}
                else if ( absVal(q+1) <= 1e-6){Antvort += 1; VortexP[0][i][j][k] = 1;}
                

                //Plaquette 2
                q = 1/(2*pi) * ( correction(fmodulo(theta_iz-theta_i,2*pi)) + 
                correction(fmodulo(theta_izy-theta_iz,2*pi)) + 
                correction(fmodulo(theta_iy-theta_izy,2*pi)) +
                correction(fmodulo(theta_i-theta_iy,2*pi)) );
                if ( absVal(q-1) <= 1e-6) {Vort += 1; VortexP[1][i][j][k] = 1;}
                else if ( absVal(q+1) <= 1e-6){Antvort += 1; VortexP[1][i][j][k] = 1;}

                //Plaquette 3
                q = 1/(2*pi) * ( correction(fmodulo(theta_ix-theta_i,2*pi)) + 
                correction(fmodulo(theta_ixz-theta_ix,2*pi)) + 
                correction(fmodulo(theta_iz-theta_ixz,2*pi)) +
                correction(fmodulo(theta_i-theta_iz,2*pi)) );
                if ( absVal(q-1) <= 1e-6) {Vort += 1; VortexP[2][i][j][k] = 1;}
                else if ( absVal(q+1) <= 1e-6){Antvort += 1; VortexP[2][i][j][k] = 1;}
            }
        }
    }  
}

inline int findPlaq(int x){
    while (PlaqLabels[x] != x){
        x = PlaqLabels[x];
    }
    return x;
}

inline void make_strings(){
    int a, b;
    for(int i = 0; i<L; i++){
        for(int j =0; j<L; j++){
            for(int k = 0; k<L; k++){
                //Plaquette 1
                for (int plaquette = 0; plaquette<3; plaquette++){
                    if (plaquette == 0){a=1; b=2;}
                    else if (plaquette == 1){a=0; b=2;}
                    else if(plaquette == 2){a=0; b=1;}
                    //If the plaquette has vorticity and it's not marked.
                    if (VortexP[plaquette][i][j][k] != 0 && MarkerPIn[plaquette][i][j][k] == 0){
                        std::vector<int> indices, plaqueta, xcoord, ycoord, zcoord;
                        //Plaquette identifiers. The plaquettes that belong to the site (i,j,k) are identified with a and b.
                        //0 --> VortexP[a][i][j][k], 1 --> VortexP[b][i][j][k], 2 --> VortexP[0][i][j][k-1],
                        //3 --> VortexP[1][modulo(j+1,L)][k], 4 --> VortexP[2][modulo(i-1,L)][j][k].
                        int plaquettes_avail = 0; //We revise whether any of the other five plaquettes of 
                        // the cube at (i,j,k) have vorticity different from zero.
                        if(VortexP[a][i][j][k] !=0){ indices.push_back(0); plaqueta.push_back(a); xcoord.push_back(i); ycoord.push_back(j);
                        zcoord.push_back(k); plaquettes_avail+=1;}
                        if(VortexP[b][i][j][k] !=0){ indices.push_back(1); plaqueta.push_back(b); xcoord.push_back(i); ycoord.push_back(j);
                        zcoord.push_back(k); plaquettes_avail+=1;}
                        if(VortexP[0][i][j][modulo(k-1,L)] !=0){ indices.push_back(2); plaqueta.push_back(0); xcoord.push_back(i); ycoord.push_back(j);
                        zcoord.push_back(modulo(k-1,L)); plaquettes_avail+=1;}
                        if(VortexP[1][i][modulo(j+1,L)][k] !=0){ indices.push_back(3); plaqueta.push_back(1); xcoord.push_back(i); ycoord.push_back(modulo(j+1,L));
                        zcoord.push_back(k); plaquettes_avail+=1;}
                        if(VortexP[2][modulo(i-1,L)][j][k] !=0){ indices.push_back(4); plaqueta.push_back(2); xcoord.push_back(modulo(i-1,L)); ycoord.push_back(j);
                        zcoord.push_back(k); plaquettes_avail+=1;}
                    
                        MarkerPIn[plaquette][i][j][k] = 1; //We mark the current plaquette as visited by a string
                        //if Plaquette at site i, j, k can't be connected with another plaquette we make a new label.
                        //std::cout <<"plaquetas disponibles " << plaquettes_avail << std::endl;
                        if (plaquettes_avail == 0){
                        LabelP[plaquette][i][j][k] = PlaqLabel; 
                        PlaqLabels[PlaqLabel] = PlaqLabel; 
                        PlaqLabel += 1;
                        }
                        //Plaquette 1 can be connected with other plaquettes.
                        else{
                            //We randomly choose one of the available plaquettes.
                            int index = rand() % plaquettes_avail; //indices[index] gives us the plaquette identifier
                            //std::cout <<"identificador " << indices[index] << std::endl;
                            //We only proceed if the plaquette is not marked.
                            bool proceed = false;
                            //In plaquettes, i.e. the ones that belong to the site (i,j,k).
                            if (indices[index] == 0 || indices[index] == 1){
                                //if plaqueta[index] is free to be connected we mark it and we proceed,
                                //otherwise we continue with the next plaquette.
                                if (MarkerPIn[plaqueta[index]][xcoord[index]][ycoord[index]][zcoord[index]]  == 0 ){
                                    MarkerPIn[plaqueta[index]][xcoord[index]][ycoord[index]][zcoord[index]]  = 1; //We mark one face of the plaquette.
                                    proceed = true;
                                }
                            }
                            //Out plaquettes, i.e. the ones that don't belong to the site (i,j,k) but to its neighbours.
                            else{
                                //Plaquette is free to be connected
                                if (MarkerPOut[plaqueta[index]][xcoord[index]][ycoord[index]][zcoord[index]] == 0){
                                    MarkerPOut[plaqueta[index]][xcoord[index]][ycoord[index]][zcoord[index]] = 1; //We mark one face of the plaquette.
                                    proceed = true;
                                }
                            
                            }
                            //std::cout << "proceed " << proceed << std::endl;
                            if (proceed == true){
                                std::vector<int> pLabels(2); 
                                pLabels[0] = findPlaq(LabelP[plaquette][i][j][k]);
                                pLabels[1] = findPlaq(LabelP[plaqueta[index]][xcoord[index]][ycoord[index]][zcoord[index]]);
                                //std::cout << "current site label " << pLabels[0] << std::endl;
                                //std::cout << "connect plaquette label " << pLabels[1] << std::endl;
                                //If both plaquettes don't belong to any string, then we generate a new label.
                                if (pLabels[0] == 0 && pLabels[1] == 0){
                                    LabelP[plaquette][i][j][k] = PlaqLabel; //Current plaquette
                                    LabelP[plaqueta[index]][xcoord[index]][ycoord[index]][zcoord[index]] = PlaqLabel; //Plaquette to be connected.
                                    PlaqLabels[PlaqLabel] = PlaqLabel;
                                    PlaqLabel += 1; 
                                }
                                //If one of the plaquettes belongs to a string, we relabel both plaquettes with the minimum proper label
                                //(this is part of the HK algorithm)
                                else if ((pLabels[0] == 0 && pLabels[1] != 0) || (pLabels[0] != 0 && pLabels[1] == 0)){
                                    //std::cout << "aqui " << std::endl;
                                    int maxLabel = *std::max_element(pLabels.begin(), pLabels.end()); //This has to be different from zero.
                                    LabelP[plaquette][i][j][k] = maxLabel;
                                    LabelP[plaqueta[index]][xcoord[index]][ycoord[index]][zcoord[index]] = maxLabel;
                                    if (pLabels[0] != 0){PlaqLabels[pLabels[0]] = maxLabel;} //std::cout << "plabel " << pLabels[0] << std::endl;}
                                    else if(pLabels[1] != 0){PlaqLabels[pLabels[1]] = maxLabel;} //std::cout << "plabel " << pLabels[1] << std::endl;}
                                }
                                else{
                                    //std::cout << "aqui 2" << std::endl;
                                    //both plaquettes belong to a string, then we must join the string.
                                    int minLabel = *std::min_element(pLabels.begin(), pLabels.end());
                                    LabelP[plaquette][i][j][k] = minLabel;
                                    LabelP[plaqueta[index]][xcoord[index]][ycoord[index]][zcoord[index]] = minLabel;
                                    PlaqLabels[pLabels[1]] = minLabel; 
                                    
                                }
                            }
                            
                        }
                        //std::cout << "" << std::endl;
                    }
                    //Plaquettes for loop
                    //std::cout << LabelP[plaquette][i][j][k] << std::endl;
                }

            //Volume for loop
            }
        }
    }
}


inline void string_density(){
    std::unordered_set<int> PLabelsSet(PlaqLabels,PlaqLabels + PlaqLabel);
    strdensity = 1.0*PLabelsSet.size()/(L*L*L);
    //std::cout << PLabelsSet.size() << std::endl;
    //std::cout << strdensity << std::endl;
    /*for (int i = 0; i<3*L*L*L; i++){
        if (i >= PlaqLabel){
            std::cout << "Debe ser cero " << PlaqLabels[i] << std::endl;
        }
        else{
            std::cout <<  PlaqLabels[i] << std::endl;
        }
    }
    */   
}


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
                VortexP[0][i][j][k] = 0; VortexP[1][i][j][k] = 0; VortexP[2][i][j][k] = 0;
                MarkerPIn[0][i][j][k] = 0; MarkerPIn[1][i][j][k] = 0; MarkerPIn[2][i][j][k] = 0;
                MarkerPOut[0][i][j][k] = 0; MarkerPOut[1][i][j][k] = 0; MarkerPOut[2][i][j][k] = 0;
                LabelP[0][i][j][k] = 0; LabelP[1][i][j][k] = 0; LabelP[2][i][j][k] = 0;
            } 
        }
    }
    for(int i = 0; i<label; i++){
        CLabels[i] = 0;
    }
    label = 0;
    for(int i = 0; i<3*L*L*L; i++){
        PlaqLabels[i] = 0;
    }
    PlaqLabel = 0;
}    

void CA_XY3d(double beta, int Ntherm, int Nmeas, int Nsteps){
    std::vector<double> SDensity(Nmeas), VORT(Nmeas), AVORT(Nmeas);
    //Thermalization//
    initialize_lattice();
    for(int i = 0; i<Ntherm; i++){
        Bonds(beta); 
        HoshenKopelman(); 
        flip(); 
        reset();
    }
    for(int i = 0; i<Nmeas; i++){
        Bonds(beta);
        HoshenKopelman(); 
        flip(); 

        Vorticity(); 
        make_strings();
        string_density();
        SDensity[i] = strdensity;
        VORT[i] = Vort;
        AVORT[i] = Antvort;
        Vort = 0; Antvort = 0; strdensity = 0;
        reset();
        for(int j = 0; j<Nsteps; j++){
            Bonds(beta);
            HoshenKopelman(); 
            flip(); 
            reset(); 
        }
    }
    strdensity = mean(SDensity); dstrdensity = Jackknife_error(SDensity,20);
    Vort = mean(VORT)/(L*L*L); dVort = Jackknife_error(VORT,20)/(L*L*L);
    Antvort = mean(AVORT)/(L*L*L); dAntvort = Jackknife_error(AVORT,20)/(L*L*L);
}

//-----------------------------------------//
int main(){    
srand(time(0));
    int Ntherm, Nmeas, Nsteps, Nbeta;
    double beta_min, beta_max; 
    //---Input data---//
    std::cout << "L= " << L << std::endl;  
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
        sprintf(Data_str, "%-30.17g%-30.17g\n", strdensity, dstrdensity);
        Datfile << Data_str;
        sprintf(Data_str, "%-30.17g%-30.17g\n", Vort, dVort);
        Datfile << Data_str;
        sprintf(Data_str, "%-30.17g%-30.17g\n", Antvort, dAntvort);
        Datfile << Data_str;
        std::cout << "String density = " << strdensity << " +- " << dstrdensity << std::endl;
        std::cout << "V = " << Vort << " +- " << dVort << std::endl;
        std::cout << "A = " << Antvort << " +- " << dAntvort << std::endl;
        Vort = 0; dVort=0; Antvort = 0; dAntvort = 0; strdensity = 0; dstrdensity=0;
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