#include <iostream>
#include <vector> 
#include <random> 
#include <time.h> 
#include <fstream>
#include <ctime>
#include <cmath> 
#include <algorithm>


double pi = 3.14159265359;
double E = 0, En = 0, En2 = 0, chi = 0, Cv = 0;
double dchi, dCv, dE, dE2;
double R; //Acceptance ratio

int L = 32;
std::vector<int> Ira{ rand() % L, rand() % L};
std::vector<int> Masha = {Ira[0], Ira[1]};

std::vector<std::vector<int>> FluxX(L, std::vector<int>(L,0));
std::vector<std::vector<int>> FluxY(L, std::vector<int>(L,0));


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

//-----------------Bessel functions stuff-----------------//

//Fourth order numerical integrator//
//b>a
double Inth4(std::vector<double> f ,std::vector<double> x){
    //For an odd number of points has slightly better precision.
    double f_int = 0, h = x[1] - x[0];
    int N = f.size();
    if (N % 2 == 0){    
        for(int i=2; i < N-1; i=i+2){
            f_int += (f[i+1] + 4*f[i] + f[i-1]);       
        }     
        f_int = f_int*h/3;
        f_int += (h/48)*(23*f[0] + 26*f[1] - 4*f[2] + 6*f[3] - 3*f[4]);
    }
    else{
        for (int i=1; i < N-1 ; i=i+2){
            f_int += (f[i+1] + 4*f[i] + f[i-1]);
        }
        f_int = f_int*h/3;
    }
    return f_int;
}

//Modified Bessel function of the first kind//
// z--> argument, n--> order, steps--> number of
//elements on the lattice when we integrate.
double BesselInu(double z, int n, int steps){
    std::vector<double> x = linspace(0,pi,steps);
    std::vector<double> f(steps);
    for (int i = 0; i<steps; i++){
        f[i] = exp(z*cos(x[i]))*cos(n*x[i])/pi;
    }
    double Inu = Inth4(f, x);
    return Inu;
}

//Derivatives of Bessel functions
//order --> number of derivative (1 or 2)
double BesselInup(double z, int n, int steps, int order){
    double Inup;
    if (order == 1){
        Inup = n * BesselInu(z, n, steps) / z + BesselInu(z,n+1,steps);
    }
    else if(order == 2){
        Inup = (BesselInu(z, n-2, steps) + 2*BesselInu(z, n, steps) + BesselInu(z, n+2, steps)) / 4.0;
    }
    return Inup;
}

struct BESSELTAB {
    std::vector<double> Inu;
    std::vector<double> Inu1;
    std::vector<double> Inu2;
};

struct BESSEL {
    double Inu;
    double Inu1;
    double Inu2;
};

BESSELTAB BesselITable(double beta,int max_nu){
    //Generates a table (list) with the values of the modified Bessel function//
    //of the first kind.//
    std::vector<double> Inu(max_nu+1), Inu1(max_nu+1), Inu2(max_nu+1);
    BESSELTAB Table;
    for (int i = 0; i < max_nu+1; i++){
        Inu[i] = BesselInu(beta, i, 1001);
        Inu1[i] = BesselInup(beta, i, 1001, 1);
        Inu2[i] = BesselInup(beta, i, 1001, 2);
    }
    Table.Inu = Inu;
    Table.Inu1 = Inu1;
    Table.Inu2 = Inu2;
    return Table;
}

BESSEL BesselTableCaller(double beta, int Jij, std::vector<double> Inu, std::vector<double> Inu1, std::vector<double> Inu2){
    BESSEL functions;
    Jij = abs(Jij); 
    int max_nu = Inu.size()-1;
    double Iv, Iv1, Iv2;
    if (Jij <= max_nu){
        Iv = Inu[Jij];
        Iv1 = Inu1[Jij];
        Iv2 = Inu2[Jij];
    }
    else{
        Iv = BesselInu(beta, Jij, 1001);
        Iv1 = BesselInup(beta, Jij, 1001, 1);
        Iv2 = BesselInup(beta, Jij, 1001, 2); 
    }
    functions.Inu = Iv;
    functions.Inu1 = Iv1;
    functions.Inu2 = Iv2;
    return functions;
}


//-----Computes the energy-----//

inline void EnergyWA(double beta, std::vector<double> Inu, std::vector<double> Inu1, std::vector<double> Inu2){
    struct BESSEL Tab;
    double Iv, Iv1, Iv2;
    for(int i = 0; i < L; i++){
        for(int j = 0; j < L; j++){
            Tab = BesselTableCaller(beta, FluxX[i][j], Inu, Inu1, Inu2);
            Iv = Tab.Inu; Iv1 = Tab.Inu1; Iv2 = Tab.Inu2;
            En -= Iv1/Iv; En2 += Iv2/Iv - (Iv1/Iv)*(Iv1/Iv);
            Tab = BesselTableCaller(beta, FluxY[i][j], Inu, Inu1, Inu2);
            Iv = Tab.Inu; Iv1 = Tab.Inu1; Iv2 = Tab.Inu2;
            En -= Iv1/Iv; En2 += Iv2/Iv - (Iv1/Iv)*(Iv1/Iv);
        }
    }
    En2 += En*En;
}

//Acceptance Ratio//
inline void AccRatio(int mu, double beta, std::vector<double> Inu, std::vector<double> Inu1, std::vector<double> Inu2){
    //This function computes the acceptance ratio based on the flux between
    //two neighbouring sites.

    struct BESSEL Tab;
    double IvFin, IvIni;
    int FluxIni, FluxFin;
    if (mu == 0){FluxIni = FluxX[Masha[0]][Masha[1]];}
    else if (mu == 1){FluxIni = -FluxX[Masha[0]][modulo(Masha[1]-1,L)];}
    else if (mu == 2){FluxIni = -FluxY[modulo(Masha[0]-1,L)][Masha[1]];}
    else if (mu == 3){FluxIni = FluxY[Masha[0]][Masha[1]];}
    FluxFin = FluxIni + 1;
    Tab = BesselTableCaller(beta, FluxIni, Inu, Inu1, Inu2);
    IvIni = Tab.Inu;
    Tab = BesselTableCaller(beta, FluxFin, Inu, Inu1, Inu2);
    IvFin = Tab.Inu;
    R = IvFin/IvIni;
}

inline void WA_Sweep(double beta, std::vector<double> Inu, std::vector<double> Inu1, std::vector<double> Inu2){
// mu = 0 --> x direction, mu = 1 --> -x direction //
// mu = 2 --> y direction, mu = 3 --> -y direction // 
    int mu = (rand() % 4);  
    double r = ((double) rand() / (RAND_MAX));
    AccRatio(mu, beta, Inu, Inu1, Inu2); //Computes R
    if (r < R){
        if (mu == 0){
            FluxX[Masha[0]][Masha[1]] += 1;
            Masha[1] = modulo(Masha[1]+1,L);
        }
        else if (mu == 1){
            FluxX[Masha[0]][modulo(Masha[1]-1,L)] -= 1;
            Masha[1] = modulo(Masha[1]-1,L);
        }
        else if (mu == 2){
            FluxY[modulo(Masha[0]-1,L)][Masha[1]] -= 1;
            Masha[0] = modulo(Masha[0]-1,L);
        }
        else if (mu == 3){
            FluxY[Masha[0]][Masha[1]] += 1;
            Masha[0] = modulo(Masha[0]+1,L);
        }
    }
}

void WA_XY2d(double beta, int Ntherm, int Nmeas, int Nsteps){
    //Worm algorithm for one value of beta//
    struct BESSELTAB Table = BesselITable(beta,10); //This computes the modified Bessel functions and its derivatives up to nu = 10.
    std::vector<double> Inu = Table.Inu, Inu1 = Table.Inu1, Inu2 = Table.Inu2;

    std::vector<double> Energy(Nmeas), Energy2(Nmeas), Chi(Nmeas);
    double tau = 0;
    //-----Thermalization------//
    for(int i = 0; i < Ntherm; i++){
        if (Masha == Ira){
            Ira = { rand() % L, rand() % L };
            Masha = {Ira[0],Ira[1]};
            tau = 0;
        }
        else{
            tau += 1;
        }
        WA_Sweep(beta, Inu, Inu1, Inu2); 
    }
    //------Measurements------//
    int count = 0;
    while (count < Nmeas){
        if (Masha == Ira){
            EnergyWA(beta, Inu, Inu1, Inu2);
            Energy[count] = En; Energy2[count] = En2;
            En = 0, En2 = 0;
            Chi[count] = tau;       
            Ira = { rand() % L, rand() % L };
            Masha = {Ira[0],Ira[1]};
            tau = 0;
            count += 1;
        }
        else{
            tau += 1;
        }     
        WA_Sweep(beta, Inu, Inu1, Inu2);
    //------Decorrelation steps------//
        for (int j = 0; j<Nsteps; j++){
            if (Masha == Ira){
                Ira = { rand() % L, rand() % L};
                Masha = {Ira[0],Ira[1]};
            }
        WA_Sweep(beta, Inu, Inu1, Inu2);
        }
    }
        E = mean(Energy);
        dE = Jackknife_error(Energy, 20);

        En2 = mean(Energy2);
        dE2 = Jackknife_error(Energy2, 20);

        Cv = beta * beta * (En2-E*E) /(L*L);
        dCv = absVal(beta * beta * (dE2 + 2*E*dE)/(L*L));

        chi = mean(Chi);
        dchi = Jackknife_error(Chi, 20);
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
    std::cout << "Step (sweeps between measurements: ";
    std::cin >> Nsteps;
    std::cout << " " << std::endl;

    std::vector<double> Betas = linspace(beta_min,beta_max, Nbeta);
    char NameCv[50], Cv_str[50], NameChi[50], Chi_str[50], NameE[50], E_str[50], NameTime[50], Time_str[50];
    sprintf(NameE, "2DXY_E_L%d_Meas%d.txt", L, Nmeas);
    sprintf(NameCv, "2DXY_Cv_L%d_Meas%d.txt", L, Nmeas);
    sprintf(NameChi, "2DXY_Chi_L%d_Meas%d.txt", L, Nmeas);
    sprintf(NameTime, "2DXY_Time_L%d_Meas%d.txt", L, Nmeas);
    std::ofstream Cvfile; std::ofstream Chifile; std::ofstream Efile; std::ofstream Timefile;
    Cvfile.open(NameCv); Chifile.open(NameChi); Efile.open(NameE); Timefile.open(NameTime);

    for (int i = 0; i < Nbeta; i++) {
        clock_t begin = clock();
        std::cout << "beta = " << Betas[i] << "  T = " << 1/Betas[i] << std::endl;
        WA_XY2d(Betas[i], Ntherm, Nmeas, Nsteps);
        sprintf(E_str, "%-30.17g%-30.17g%-30.17g\n", Betas[i], E, dE);
        sprintf(Cv_str, "%-30.17g%-30.17g%-30.17g\n", Betas[i], Cv, dCv);
        sprintf(Chi_str, "%-30.17g%-30.17g%-30.17g\n", Betas[i], chi, dchi);
        Efile << E_str;
        Cvfile << Cv_str;
        Chifile << Chi_str;
        std::cout << "E = " << E << " +- " << dE << std::endl;
        std::cout << "Cv = " << Cv << " +- " << dCv << std::endl;
        std::cout << "Chi = " << chi << " +- " << dchi << std::endl; 
        E = 0; dE= 0; Cv = 0; dCv = 0; chi = 0; dchi = 0;
        //----Computing time----//
        clock_t end = clock();
        double elapsed_secs = double(end - begin) / CLOCKS_PER_SEC;
        sprintf(Time_str, "%-30.17g%-30.17g\n", Betas[i],elapsed_secs);
        Timefile << Time_str;
        std::cout << "Time = " << elapsed_secs << " s" << std::endl;
        std::cout << "------------------------------" << std::endl;
    }
    Efile.close();
    Cvfile.close();
    Chifile.close();
    Timefile.close();
    return 0;      
}