#include <iostream>
#include <vector> 
#include <random> 
#include <time.h> 
#include <fstream>
#include <ctime>
#include <cmath> 
#include <algorithm>


double pi = 3.14159265359;
double E = 0, En2 = 0, M = 0, M2=0, chi = 0, Cv = 0, Vort = 0, Antvort = 0, beta;
double dchi, dM, dM2, dCv, dE, dE2, dVort, dAntvort;
double p, reflected_phi;//Acceptance ratio and reflected angle
int label = 0;

constexpr int L = 24;
constexpr  int maxsize = L*L;
int CLabels[maxsize]; //Cluster labels.
std::vector<std::vector<double>> SpinLattice(L, std::vector<double>(L,0));
std::vector<std::vector<double>> SpinLatticeR(L, std::vector<double>(L,0));

std::vector<std::vector<int>> xBonds(L, std::vector<int>(L,0));
std::vector<std::vector<int>> yBonds(L, std::vector<int>(L,0));
std::vector<std::vector<int>> Labels(L, std::vector<int>(L,0)); //Labels of each site.


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

//---------------Logspace (similar to python)----------------------//
std::vector<int> logspace(double min, double max, int n) {
    std::vector<int> logspace(n);
    double h = (max*1.0 - min*1.0) / (n - 1);
    for (int i = 0; i < n; ++i) {
        logspace[i] = (int) pow(10.0, min + i * h); 
    }
    return logspace;
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

//n modulus m 
double fmodulo(double n, double m) {
    double res;
    if (n < 0) {
        res = fmod(n + m, m);
    }
    else {
        res = fmod(n,m);
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
inline void Bonds(){
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
   
inline int find(int x){
    while (CLabels[x] != x){
        x = CLabels[x];
    }
    return x;
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
            index = indexfinder(CLabels[find(Labels[i][j])]);    //antes era find(Labels[i][j]) -1                  
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

//V1 --> Measures the vortex density for a specific time tauQ + t. There isn't a fast quench mode this.
//V2 --> Measures the vortex density for every time in the range (0,2t auQ). mode = 1 slow quench, = 2 fast quench.
//------Version 1 --------//
inline void CA_XY2dV1(int tauQ, int t,int Ntherm,int Nmeas){
    // TauQ --> inverse cooling rating
    // Ntherm --> Thermalization steps
    // t --> we measure the density of defects at t+tauQ-1 t<tauQ
    std::vector<double>  VORT, AVORT;
    double Tc = 0.89;
    for(int k=0; k<Nmeas; k++){
        beta = 1/(2*Tc);
        //Random initial configuration
        for(int l = 0; l<L; l++){
            for(int m = 0; m<L; m++){
                SpinLattice[l][m] = rand_range(-pi,pi);
            }
        }
        //Thermalization//
        for(int i = 0; i<Ntherm; i++){
            Bonds(); 
            HoshenKopelman(); 
            flip(); 
            reset();
        }
        for(int i = 0; i<tauQ + t ; i++){
            beta = 1.0/( Tc*(1.0-( i*1.0-tauQ*1.0  )/(tauQ*1.0)) );
            Bonds(); 
            HoshenKopelman(); 
            flip(); 
            reset(); 
            if (i==tauQ + t - 1){
                Vorticity(); 
                VORT.push_back(Vort);
                AVORT.push_back(Antvort);
                Vort = 0; Antvort = 0;
            }
        }
    }
    Vort = mean(VORT)/(L*L); dVort = Jackknife_error(VORT,20)/(L*L);
    Antvort = mean(AVORT)/(L*L); dAntvort = Jackknife_error(AVORT,20)/(L*L);
}

//------Version 2 --------//
inline void CA_XY2dV2(int tauQ, int Ntherm, int Nmeas,int mode, double T){
    // TauQ --> inverse cooling rating
    // Ntherm --> Thermalization steps
    std::vector<double> VORT, AVORT;
    double Tc = 0.89; 
    std::vector<std::vector<double>> Mediciones(2*tauQ, std::vector<double>(Nmeas,0));
    for(int k = 0; k<Nmeas; k++){
        beta = 1/(2*Tc);
        //Random initial condition//
        for(int l = 0; l<L; l++){
            for(int m = 0; m<L; m++){
                SpinLattice[l][m] = rand_range(-pi,pi);
            }
        }
        //Slow quench//
        if (mode == 1){
            //Thermalization//
            for(int i = 0; i<Ntherm; i++){
                Bonds(); 
                HoshenKopelman(); 
                flip(); 
                reset(); 
            }
            for(int i = 0; i<2*tauQ; i++){
                beta = 1.0/( Tc*(1.0-( i*1.0-tauQ*1.0  )/(tauQ*1.0)) );
                Bonds(); 
                HoshenKopelman(); 
                flip(); 
                reset();
                Vorticity(); 
                Mediciones[i][k] = Vort;
                Vort=0, Antvort = 0;
            }   
        }
        //Fast quench//
        else if (mode == 2){
            beta = 1/T;
            for(int i = 0; i<2*tauQ; i++){
                Bonds(); 
                HoshenKopelman(); 
                flip(); 
                reset(); 
                Vorticity(); 
                Mediciones[i][k] = Vort;
                Vort=0, Antvort = 0;
            }  
        }
       else{std::cout << "choose a valid mode (1 or 2)" << std::endl;} 
    }
    char NameData[50], Data_str[100];
    if (mode == 1){sprintf(NameData, "2DXY_L%d_tdep_tauQ%d.txt", L, tauQ);}
    else if (mode == 2){
        float z =  T*pow(10,floor(log10(T)));
        int y = (int) z;
        sprintf(NameData, "2DXY_L%d_tdep_T%0*d.txt", L, (int) abs(floor(log10(T)))+1,y);
        }
    std::ofstream Datfile;
    Datfile.open(NameData);
    for(int i=0;i<2*tauQ;i++ ){
        Vort = mean(Mediciones[i])/(L*L); dVort = Jackknife_error(Mediciones[i],20)/(L*L);
        sprintf(Data_str, "%-30d%-30.17g%-30.17g\n",i, Vort, dVort);
        Datfile << Data_str;
    }
    Datfile.close();
}
//-----------------------------------------//

int main(){    
    srand(time(0));
    int Ntherm=1, Nmeas=1, Ntau=1, version, t=1, tauQ=1;
    double tauQ_min=1, tauQ_max=1, T=1; 
    //---Input data---//
    std::cout << "-----Cluster algorithm-----" << std::endl;
    std::cout << "Mode: " <<  std::endl;
    std::cout << "V1 --> Slow quench, measures the vortex density for a specific time tauQ + t"<<  std::endl;
    std::cout << "V2 --> Slow quench, measures the vortex density for every time in the range (0,2 tauQ)"<<  std::endl;
    std::cout << "V3 --> Fast quench, measures the vortex density for every time in the range (0,2 tauQ) at temperature T" <<  std::endl;
    std::cout << "V4 --> Slow quench, measures the vortex density at t = 2 tauQ" <<  std::endl;
    std::cin >> version;
    if (version == 1 || version == 4){
        std::cout << "tauQ min (order of magnitude): "; 
        std::cin >> tauQ_min;
        std::cout << "tauQ max (order of magnitude): ";
        std::cin >> tauQ_max;
        std::cout << "Number of taus: ";
        std::cin >> Ntau;
        std::cout << "Instant of time to measure (t<tauQ min): ";
        std::cin >> t;
        std::cout << "Thermalization: ";
        std::cin >> Ntherm;
    }
    if (version == 2){
        std::cout << "Thermalization: ";
        std::cin >> Ntherm;
    }
    if (version == 2 || version == 3){
        std::cout << "tauQ: ";
        std::cin >> tauQ;
    }
    std::cout << "Measurements: ";
    std::cin >> Nmeas;
    if (version == 3){
        std::cout << "Temperature: ";
        std::cin >> T;
    }
    std::cout << " " << std::endl;
    //-------TauQ-------//
    std::vector<int> TAUS(Ntau);
    if (Ntau == 1){ 
        TAUS = {(int) tauQ_min}; 
    }
    else{ 
        TAUS = logspace(tauQ_min,tauQ_max, Ntau); //Logarithmic step.
    }
    //--------------------//
    char NameData[50], Data_str[100];
    if (version == 1 || version == 4){
        for (int i = 0; i < Ntau; i++) {
            if (version == 4){t = TAUS[i];}
            sprintf(NameData, "2DXY_L%d_tauQ%d_t%d.dat", L, TAUS[i],t);
            std::ofstream Datfile;
            Datfile.open(NameData); 

            clock_t begin = clock();
            std::cout << "tauQ = " << TAUS[i] << " t = " << t <<  std::endl;
            CA_XY2dV1(TAUS[i], t,Ntherm, Nmeas);   
            clock_t end = clock(); 

            sprintf(Data_str,"%-30d%-30d%-30d%-30d%-30d\n",Ntherm,TAUS[i],Nmeas,L,t);
            Datfile << Data_str; 
            sprintf(Data_str, "%-30.17g%-30.17g\n", Vort, dVort);
            Datfile << Data_str;
            std::cout << "V = " << Vort << " +- " << dVort << std::endl;
            std::cout << "A = " << Antvort << " +- " << dAntvort << std::endl;
            Vort = 0; dVort=0; Antvort = 0; dAntvort = 0;
            //----Computing time----// 
            double elapsed_secs = double(end - begin) / CLOCKS_PER_SEC; 
            sprintf(Data_str, "%-30.17g", elapsed_secs);
            Datfile << Data_str; 
            std::cout << "Time = " << elapsed_secs << " s" << std::endl;
            std::cout << "------------------------------" << std::endl;
            Datfile.close();
        }
    }

    else if (version == 2){
        clock_t begin = clock();
        CA_XY2dV2(tauQ, Ntherm, Nmeas, 1, T);//slow quench
        clock_t end = clock();
        double elapsed_secs = double(end - begin) / CLOCKS_PER_SEC; 
        std::cout << "Time = " << elapsed_secs << " s" << std::endl;
    }
    else if (version == 3){
        clock_t begin = clock();
        CA_XY2dV2(tauQ, Ntherm, Nmeas, 2, T);//Fast quench
        clock_t end = clock();
        double elapsed_secs = double(end - begin) / CLOCKS_PER_SEC; 
        std::cout << "Time = " << elapsed_secs << " s" << std::endl;
    }
 
    return 0;
}