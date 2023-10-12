#include <iostream>
#include <vector> 
#include <random> 
#include <time.h> 
#include <fstream>
#include <ctime>
#include <cmath> 
#include <algorithm>


double pi = 3.14159265359;
double Vort = 0, Antvort = 0;
double dVort, dAntvort;
double p, beta;

constexpr int L = 100;
constexpr  int maxsize = L*L;
std::vector<std::vector<double>> SpinLattice(L, std::vector<double>(L,0));
std::vector<std::vector<double>> SpinLatticeR(L, std::vector<double>(L,0));



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
    


//-----Chooses a new angle-----//
inline void Reflection(std::vector<int> site){
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

//V1 --> Measures the vortex density for a specific time tauQ + t. There isn't a fast quench mode this.
//V2 --> Measures the vortex density for every time in the range (0,2t auQ). mode = 1 slow quench, = 2 fast quench.
//Version 1//
inline void Met_XY2dV1(int tauQ, int t,int Ntherm,int Nmeas){
    // TauQ --> inverse cooling rating
    // Ntherm --> Thermalization steps
    // t --> we measure the density of defects at t+tauQ-1 t<tauQ
    std::vector<double> VORT, AVORT;
    double Tc = 0.89;
    if (t>tauQ){std::cout << "t must be smaller than tauQ" << std::endl;}
    for(int k = 0; k<Nmeas; k++){
        beta = 1/(2*Tc);
        //Random initial configuration//
        for(int l = 0; l<L; l++){
            for(int m = 0; m<L; m++){
                SpinLattice[l][m] = rand_range(-pi,pi);
            }
        }
        //Thermalization//
        for(int i = 0; i<Ntherm*L*L; i++){
            Reflection({rand() % L,rand() % L });  
        }
        for(int i = 0; i<tauQ+t; i++){
            beta = 1.0/( Tc*(1.0-( i*1.0-tauQ*1.0  )/(tauQ*1.0)) );
            for(int l = 0; l<L; l++){
                for(int m = 0; m<L; m++){
                    Reflection({rand() % L,rand() % L }); 
                }
            }
            if (i == tauQ+t-1 ){
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

//Version 2//
inline void Met_XY2dV2(int tauQ, int Ntherm, int Nmeas,int mode, double T){
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
            for(int i = 0; i<Ntherm*L*L; i++){
                Reflection({rand() % L,rand() % L });  
            }
            for(int i = 0; i<2*tauQ; i++){
                beta = 1.0/( Tc*(1.0-( i*1.0-tauQ*1.0  )/(tauQ*1.0)) );
                for(int l = 0; l<L; l++){
                    for(int m = 0; m<L; m++){
                        Reflection({rand() % L,rand() % L }); 
                    }
                }
                Vorticity(); 
                Mediciones[i][k] = Vort;
                Vort=0; Antvort = 0;
            }   
        }
        //Fast quench//
        else if (mode == 2){
            beta = 1/T;
            for(int i = 0; i<2*tauQ; i++){
                for(int l = 0; l<L; l++){
                    for(int m = 0; m<L; m++){
                        Reflection({rand() % L,rand() % L }); 
                    }
                }
                Vorticity(); 
                Mediciones[i][k] = Vort;
                Vort=0; Antvort = 0;
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
    std::cout << "-----Metropolis-----" << std::endl;
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
            Met_XY2dV1(TAUS[i], t,Ntherm, Nmeas);   
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
        Met_XY2dV2(tauQ, Ntherm, Nmeas, 1, T);//slow quench
        clock_t end = clock();
        double elapsed_secs = double(end - begin) / CLOCKS_PER_SEC; 
        std::cout << "Time = " << elapsed_secs << " s" << std::endl;
    }
    else if (version == 3){
        clock_t begin = clock();
        Met_XY2dV2(tauQ, Ntherm, Nmeas, 2, T);//Fast quench
        clock_t end = clock();
        double elapsed_secs = double(end - begin) / CLOCKS_PER_SEC; 
        std::cout << "Time = " << elapsed_secs << " s" << std::endl;
    }
 
    return 0;
}
