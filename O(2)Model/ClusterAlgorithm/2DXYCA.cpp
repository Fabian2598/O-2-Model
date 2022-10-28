#include <iostream>
#include <vector> 
#include <random> 
#include <time.h> 
#include <fstream>
#include <ctime>
#include <cmath> 
#include <algorithm>


double pi = 3.14159265359;
double E = 0, En2 = 0, chi = 0, Cv = 0, rho = 0;
double dchi, dCv, dE, dE2, drho=0;
double p, reflected_phi;//Acceptance ratio and reflected angle
int label = 0;

constexpr int L = 8;
constexpr  int maxsize = L*L;
int CLabels[maxsize]; //Cluster labels.
std::vector<std::vector<double>> SpinLattice(L, std::vector<double>(L,0));
std::vector<std::vector<double>> SpinLatticeR(L, std::vector<double>(L,0));

std::vector<std::vector<int>> xBonds(L, std::vector<int>(L,0));
std::vector<std::vector<int>> yBonds(L, std::vector<int>(L,0));
std::vector<std::vector<int>> Labels(L, std::vector<int>(L,0)); //Labels of each site.



void printmatrix(std::vector<std::vector<int>> A){
    for (int i = 0; i<L ; i++){
        for(int j = 0; j<L; j++){
            std::cout << A[i][j] << " ";
        }   
        std::cout << "" << std::endl;
    }
}

void printlabels(int max){
    for(int i = 0; i<max ; i++){
        std::cout << CLabels[i] << " ";
    }
    std::cout << "" << std::endl;
}

double mean(std::vector<double> x) {
    double prom = 0;
    for (double i : x) {
        prom += i;
    }
    prom = prom / x.size();
    return prom;
}

double absVal(double z){
    if (z < 0){
        return -z;
    }
    else{
        return z;
    }
}

inline double rand_range(double a, double b){
    //generates a random double number in the inteval [a,b] a = min, b = max
    double cociente = ((double) rand() / (RAND_MAX));
    double x = (b-a) * cociente + a;
    return x;
}

//----------Jackknife---------//
std::vector<double> samples_mean(std::vector<double> dat, int bin) {
    std::vector<double> samples_mean(bin);
    int dat_bin = dat.size() / bin;
    double prom = 0;
    for (int i = 0; i < bin; i++) {
        for (int k = 0; k < bin; k++) {
            for (int j = k * dat_bin; j < k * dat_bin + dat_bin; j++) {
                if (k != i) {
                    prom += dat[j];
                }
            }
        }
        prom = prom / (dat.size() - dat_bin);
        samples_mean[i] = prom;
        prom = 0;
    }
    return samples_mean;
}

double Jackknife_error(std::vector<double> dat, int bin) {
    double error = 0;
    std::vector<double> sm = samples_mean(dat, bin);
    double normal_mean = mean(dat);
    for (int m = 0; m < bin; m++) {
        error += pow((sm[m] - normal_mean), 2);
    }
    error = sqrt(error * (bin - 1) / bin);
    return error;
}


double Jackknife(std::vector<double> dat, std::vector<int> bins) {
    std::vector<double> errores(bins.size());
    double error;
    for (int i = 0; i < bins.size(); i++) {
        errores[i] = Jackknife_error(dat, bins[i]);
    }
    error = *std::max_element(errores.begin(), errores.end());
    return error;
}
//-------------End of Jackknife--------------//

//---------------Linspace (similar to python)----------------------//
std::vector<double> linspace(double min, double max, int n) {
    std::vector<double> linspace;
    double h = (max - min) / (n - 1);
    for (int i = 0; i < n; ++i) {
        linspace.insert(linspace.begin() + i, min + i * h); 
    }
    return linspace;
}

//n modulus m 
int modulo(int n, int m) {
    int res;
    if (n < 0) {
        res = (n + m) % m;
    }
    else {
        res = n % m;
    }
    return res;
}

//-----Computes the energy of a configuration-----//
inline void EnergyCA(){
    for (int i = 0; i<L; i++){
        for (int j = 0; j<L; j++){
            E += -cos(SpinLattice[i][j] - SpinLattice[i][modulo(j+1,L)]) - cos(SpinLattice[i][j] - SpinLattice[modulo(i+1,L)][j]);
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
    chi = sqrt(sx*sx + sy*sy)/(L*L);
}

//-----Reflects a vector with respect to the Wolff line-----//
inline void Reflection(std::vector<int> site, double rphi){
    //site --> site of the SpinLattice where to be reflected
    //r --> angle of a 2D random vector 
    // returns the reflected angle with respect to the Wolff line
    double phi = SpinLattice[site[0]][site[1]];
    double sx = cos(phi), sy = sin(phi), rx = cos(rphi), ry = sin(rphi); //Vector components
    double sr_x = sx-2*( sx*rx + sy*ry )*rx, sr_y = sy-2*( sx*rx + sy*ry )*ry; //Reflected components
    //We compute the angle depending on the quadrant where the reflected spin is directed
    if (sr_x> 0 && sr_y>0){ reflected_phi = atan(sr_y/sr_x); } //First quadrant
    else if (sr_x<0 && sr_y>0){reflected_phi = pi - atan(-sr_y/sr_x);} //Second quadrant
    else if (sr_x<0 && sr_y<0){reflected_phi = pi + atan(sr_y/sr_x);} //Third quadrant
    else if (sr_x>0 && sr_y<0){reflected_phi =  2*pi - atan(-sr_y/sr_x); }//Fourth quadrant
    else if (absVal(sr_x) <= 1e-10){
        if (sr_y>0){reflected_phi = pi/2;}
        else{reflected_phi = 3*pi/2;}
    }
    else if (absVal(sr_y) <= 1e-10){
        if (sr_x>0){reflected_phi = 0;}
        else{reflected_phi = pi;}
    }
}

//-----Function that generates bonds between neighbouring sites-----//
inline void Bonds(double beta){
    // modifies xBonds (bonds in the x direction)
    //          yBonds (bonds in the y direction)
    //          SpinLatticeR (Lattice with the angles of the reflected spins)
    double rphi = rand_range(0.0, 2*pi);
    double rv_x = cos(rphi), rv_y = sin(rphi); //vector associated to the rphi angle
    for(int i = 0; i<L; i++){
        for(int j = 0; j<L; j++){
            Reflection({i,j}, rphi);//We reflect the spin
            SpinLatticeR[i][j] = reflected_phi; 
            double sx = cos(SpinLattice[i][j]), sy = sin(SpinLattice[i][j]);
            double s_dot_r = rv_x*sx + rv_y*sy; //dot product between r and s.

            //---Creating bond with the right neigbour---//
            double sv_x = cos(SpinLattice[i][modulo(j+1,L)]), sv_y = sin(SpinLattice[i][modulo(j+1,L)]); //Neighbour spin
            
            double sv_dot_r = rv_x*sv_x + rv_y*sv_y;
            if (s_dot_r*sv_dot_r >=0){
                p = 1-exp(-2*beta*s_dot_r*sv_dot_r);
            }
            else{
                p = 0;
            }
            double R = ((double) rand() / (RAND_MAX));
            if (R<p){xBonds[i][j] = 1;}
            //---Creating bond with the lower neigbour---//
            sv_x = cos(SpinLattice[modulo(i+1,L)][j]); sv_y = sin(SpinLattice[modulo(i+1,L)][j]);
            sv_dot_r = rv_x*sv_x + rv_y*sv_y;
            if (s_dot_r*sv_dot_r >=0){
                p = 1-exp(-2*beta*s_dot_r*sv_dot_r);
            }
            else{
                p = 0;
            }
            R = ((double) rand() / (RAND_MAX));
            if (R<p){yBonds[i][j] = 1;}
        }
    }

}
   
//-----Hoshen-Kopelman algorithm-----//
inline void HoshenKopelman(){
    //xBonds --> If xBonds[i][j] = 1 then there is a bond between the sites
    //[i][j] and [i][(j+1)%L] (positive x direction). If xBonds[i][j] = 0 there 
    //isn't any bond.
    //yBonds --> If yBonds[i,j] = 1 then there is a bond between the sites 
    //[i][j] and [(i+1)%L][j] (negative y direction). If yBonds[i][j] = 0 there 
    //isn't any bond.
    //i --> row, j --> column.
    //modifies: CLabels --> Different labels of the clusters
    //          Labels --> L x L matrix with the labels at each site
    for(int i = 0; i<L; i++){
        for(int j = 0; j<L; j++){
            std::vector<int> row, column;
            int bonds = 0;
            if (xBonds[i][j] == 1){
                bonds += 1; row.push_back(i); column.push_back(modulo(j+1,L));
            }
            if (xBonds[i][modulo(j-1,L)] == 1){
                bonds += 1; row.push_back(i); column.push_back(modulo(j-1,L));
            }
            if (yBonds[i][j] == 1){
                bonds += 1; row.push_back(modulo(i+1,L)); column.push_back(j);
            }
            if (yBonds[modulo(i-1,L)][j] == 1){
                bonds += 1; row.push_back(modulo(i-1,L)); column.push_back(j);
            }

            std::vector<int> NLabels;
            //Neighbors that are in a cluster already.
            if (bonds != 0){
                for(int k=0; k<bonds; k++){
                    if (Labels[row[k]][column[k]] != 0){NLabels.push_back(CLabels[Labels[row[k]][column[k]]-1]);}
                }
            }
            

            if (Labels[i][j] == 0){
                if (NLabels.size() == 0){
                    label += 1;
                    CLabels[label-1] = label;
                    Labels[i][j] = label;
                    if (bonds != 0){
                        for(int k=0; k<bonds; k++){Labels[row[k]][column[k]] = label;}
                    }
                    
                }
                else{
                    //Finding the minimum label
                    int min_label = NLabels[0];
                    for(int b = 1; b<NLabels.size(); b++ ){
                        if (min_label > NLabels[b]){ min_label = NLabels[b]; }
                    }
                    //-------------------------//
                    Labels[i][j] = min_label;
                    CLabels[Labels[i][j]-1] = min_label;
                    for (int k = 0; k<bonds; k++){
                        if (Labels[row[k]][column[k]] == 0) {Labels[row[k]][column[k]] = Labels[i][j];}
                    }
                    for (int k = 0; k<label; k++){
                        if (std::find(NLabels.begin(), NLabels.end(), CLabels[k]) != NLabels.end()){
                            CLabels[k] = min_label;
                        }
                    }
                }
            }

            else{
                if (NLabels.size() == 0){
                    if (bonds != 0){
                        for(int k=0; k<bonds; k++){Labels[row[k]][column[k]] = Labels[i][j];}
                    }
                }
                else{
                    //Finding the minimum label
                    int min_label = NLabels[0];
                    for(int b = 1; b<NLabels.size(); b++ ){
                        if (min_label > NLabels[b]){ min_label = NLabels[b]; }
                    }
                    //-------------------------//
                    CLabels[Labels[i][j]-1] = min_label;
                    for (int k = 0; k<bonds; k++){
                        if (Labels[row[k]][column[k]] == 0) {Labels[row[k]][column[k]] = Labels[i][j];}
                    }
                    for (int k = 0; k<label; k++){
                        if (std::find(NLabels.begin(), NLabels.end(), CLabels[k]) != NLabels.end()){
                            CLabels[k] = min_label;
                        }
                    }

                }

            }
            
        }
    }
}

inline int indexfinder(int element){
    int i = 0;
    while (i < L*L){
        if (CLabels[i] == element){
            break;
        }
        i++;
    }
    return i;
}

inline void flip(){
    std::vector<double> probs;
    int index;
    for(int i = 0; i<label; i++){
        double R = ((double) rand() / (RAND_MAX));
        probs.push_back(R);
    }
    for(int i = 0; i<L; i++){
        for(int j = 0; j<L; j++){
            index = indexfinder(CLabels[Labels[i][j]-1]);                     
            if (probs[index] < 0.5){
                SpinLattice[i][j] = SpinLatticeR[i][j];
            }
        }
    }

}

inline void reset(){
    for(int i = 0; i<L; i++){
        for(int j=0; j<L; j++){
            xBonds[i][j] = 0; yBonds[i][j] = 0; Labels[i][j] = 0; SpinLatticeR[i][j] = 0;
        }
    }
    for(int i = 0; i<L*L; i++){
        CLabels[i] = 0;
    }
    label = 0;
}    

void CA_XY2d(double beta, int Ntherm, int Nmeas, int Nsteps){
    std::vector<double> Energy(Nmeas), Energy2(Nmeas), Magn(Nmeas);
    //Thermalization//
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
        Magn[i] = chi;
        E = 0; chi = 0;
        for(int j = 0; j<Nsteps; j++){
            Bonds(beta); //Computes the bonds and the reflected lattice.
            HoshenKopelman(); //Identifies the clusters.
            flip(); //Flips the spins with probability 1/2.
            reset(); //Resets bonds and labels.
        }
    }
    E = mean(Energy); dE = Jackknife_error(Energy, 20);
    En2 = mean(Energy2); dE2 = Jackknife_error(Energy2, 20);
    Cv = beta * beta * (En2-E*E) /(L*L); dCv = absVal(beta * beta * (dE2 + 2*E*dE)/(L*L));
    chi = mean(Magn); dchi = Jackknife_error(Magn, 20);
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
        else{sprintf(NameData, "2DXY_L%d_Meas%d_b%d.dat", L, Nmeas,beta);}
        
        std::ofstream Datfile;
        Datfile.open(NameData);
        clock_t begin = clock();
        std::cout << "beta = " << Betas[i] << "  T = " << 1/Betas[i] << std::endl;
        CA_XY2d(Betas[i], Ntherm, Nmeas, Nsteps);
        sprintf(Data_str,"%-30d%-30d%-30d%-30d%-30.17g\n",Ntherm,Nmeas,Nsteps,L,Betas[i]);
        Datfile << Data_str;
        sprintf(Data_str, "%-30.17g%-30.17g\n", E, dE);
        Datfile << Data_str;
        sprintf(Data_str, "%-30.17g%-30.17g\n", Cv, dCv);
        Datfile << Data_str;
        sprintf(Data_str, "%-30.17g%-30.17g\n", chi, dchi);
        Datfile << Data_str;
        sprintf(Data_str, "%-30.17g%-30.17g\n", rho, drho);
        Datfile << Data_str;
        
        std::cout << "E = " << E << " +- " << dE << std::endl;
        std::cout << "Cv = " << Cv << " +- " << dCv << std::endl;
        std::cout << "Chi = " << chi << " +- " << dchi << std::endl; 
        E = 0; dE= 0; Cv = 0; dCv = 0; chi = 0; dchi = 0; rho = 0; drho = 0;
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