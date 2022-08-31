#include <iostream>
#include <vector> 
#include <random> 
#include <time.h> 
#include <fstream>
#include <ctime>
#include <cmath> 
#include <algorithm>


double mean(std::vector<double> x) {
    double prom = 0;
    for (double i : x) {
        prom += i;
    }
    prom = prom / x.size();
    return prom;
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

std::vector<int> factores(int N) {
    std::vector<int> fact;
    for (int i = 2; i < N + 1; i++) {
        if (N % i == 0) {
            fact.push_back(i);
        }
    }
    return fact;
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

//Generates a random numbre between x_min and x_max
double x_rand(double x_min, double x_max) {
    double x;
    double cociente;
    int i = 0;
    while (i != 1) {
        cociente = (((double)rand()) / (double)rand());
        if (cociente <= 1) {
            x = (x_max - x_min) * cociente + x_min;
            return x;
        }
    }
    return 0;
}

//---------------Linspace (similar to python)----------------------//
std::vector<double> linspace(double min, double max, int n) {
    std::vector<double> linspace;
    double h = (max - min) / (n - 1);
    for (int i = 0; i < n; ++i) {
        linspace.insert(linspace.begin() + i, min + i * h);
    }
    return linspace;
}

//n módulo m 
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

//Matrix of zeros initializer//
std::vector<std::vector<int>> ZeroMat(int L) {
    std::vector<std::vector<int>> Matrix;
    Matrix.resize(L, std::vector<int>(L, 0));
    for (int i = 0; i < L; i++) {
        for (int j = 0; j < L; j++) {
            Matrix[i][j] = 0;
        }
    }
    return Matrix;
}

//Acceptance ratio//
double AccRatio(bool CreateBond, int mu, int L, double beta, std::vector<int> Masha, std::vector<std::vector<int>> BondsX, std::vector<std::vector<int>> BondsY){
    // mu = 0 --> x direction, mu = 1 --> -x direction //
    // mu = 2 --> y direction, mu = 3 --> -y direction //
    // If we create a new bond then R = beta / (nb+1)  //
    // If we delete a bond then R = nb / beta          //
    double R, nb;
    if (mu == 0){nb = BondsX[modulo(Masha[0],L)][modulo(Masha[1],L)];} //Number of bonds at Masha in the x direction

    else if (mu == 1){ nb = BondsX[modulo(Masha[0],L)][modulo(Masha[1]-1,L)];} //Number of bonds at Masha in the -x direction

    else if (mu == 2){ nb = BondsY[modulo(Masha[0]-1,L)][modulo(Masha[1],L)];} //Number of bonds at Masha in the y direction

    else if (mu == 3){ nb = BondsY[modulo(Masha[0],L)][modulo(Masha[1],L)];} //Number of bonds at Masha in the -y direction

    if (CreateBond == true){R = beta / (nb + 1);}
    else{R = nb / beta;}
    
    return R;
}

//-----Data structure for a Worm algorithm sweep-----//
struct SWEEP {
    std::vector<std::vector<int>> BondsX;
    std::vector<std::vector<int>> BondsY;
    std::vector<int> Masha;
    double Nb;
};

//----Sweep of the Worm Algorithm----//
SWEEP WA_Sweep(int L, double beta, double Nb, std::vector<int> Masha, std::vector<std::vector<int>> BondsX, std::vector<std::vector<int>> BondsY){
    SWEEP sweep;
    //------------Bond creation------------//
    bool CreateBond = false;
    if ((rand() % 2) == 0){CreateBond = true;}
    // mu = 0 --> x direction, mu = 1 --> -x direction //
    // mu = 2 --> y direction, mu = 3 --> -y direction //
    int mu = (rand() % 4);
    double R = AccRatio(CreateBond, mu, L, beta, Masha, BondsX, BondsY), r = x_rand(0, 1), nb;
    if(CreateBond == true){nb = 1.0;}
    else{nb = -1.0;}
    if (r < R){
        if (mu == 0){
            BondsX[modulo(Masha[0],L)][modulo(Masha[1],L)] += nb;
            Masha[1] = modulo((Masha[1]+1),L); //Move Masha one site to the right
        }
        else if (mu == 1){
            BondsX[modulo(Masha[0],L)][modulo(Masha[1]-1,L)] += nb;
            Masha[1] = modulo((Masha[1]-1),L); //Move Masha one site to the left
        }
        else if (mu == 2){
            BondsY[modulo(Masha[0]-1,L)][modulo(Masha[1],L)] += nb; 
            Masha[0] = modulo(Masha[0]-1,L); //Move Masha one site up
        }
        else if (mu == 3){
            BondsY[modulo(Masha[0],L)][modulo(Masha[1],L)] += nb; 
            Masha[0] = modulo(Masha[0]+1,L); //Move Masha one site down
        }
        Nb += nb; //Modify the total number of bonds.
    } 
    sweep.Masha = Masha; sweep.Nb = Nb; sweep.BondsX = BondsX; sweep.BondsY = BondsY;
    return sweep;
}

//-----Data structure for a Worm algorithm sweep-----//
struct RESULTS {
    double Cv;
    double dCv;
    double Chi;
    double dChi;
};

RESULTS WA_IsingModel(double beta, int L, int Ntherm, int Nmeas, int Nsteps){
    RESULTS results;
    SWEEP sweep;
    std::vector<int> Ira{ rand() % L, rand() % L };
    std::vector<int> Masha = {Ira[0], Ira[1]};
    std::vector<std::vector<int>> BondsX = ZeroMat(L), BondsY = ZeroMat(L);
    double Z = 0.0, Nb = 0.0;
    std::vector<double> G(L*L, 0.0), CHI(Nmeas), NB(Nmeas), NB2_NB(Nmeas);
    double Chi, Cv, dChi, dCv;
    //----Thermalization sweeps----//
    for (int i = 0; i<Ntherm; i++){
        if (Masha == Ira){
            Ira = { rand() % L, rand() % L };
            Masha = {Ira[0],Ira[1]};
        }
        sweep = WA_Sweep(L,beta,Nb,Masha,BondsX,BondsY);
        Masha = sweep.Masha; Nb = sweep.Nb; BondsX = sweep.BondsX; BondsY = sweep.BondsY;
    }

    int Count = 0;
    //----Measurements----//
    while (Count <Nmeas){
        int SiteIra = Ira[1] + Ira[0]*L;
        int SiteMasha = Masha[1] + Masha[0]*L;
        G[abs(SiteIra-SiteMasha)] += 1;
        if (Masha == Ira){
            Z += 1; 
            Ira = { rand() % L, rand() % L };
            Masha = {Ira[0],Ira[1]};
            NB[Count] = Nb * 1.0;
            NB2_NB[Count] = (Nb * Nb * 1.0 - Nb * 1.0);
            CHI[Count]= mean(G)/Z * L * L;
            Count += 1;
        }
        sweep = WA_Sweep(L,beta,Nb,Masha,BondsX,BondsY);
        Masha = sweep.Masha; Nb = sweep.Nb; BondsX = sweep.BondsX; BondsY = sweep.BondsY;
        //----Decorrelation steps-----// ¿?
        for (int j=0; j<Nsteps; j++){
            int SiteIra = Ira[1] + Ira[0]*L;
            int SiteMasha = Masha[1] + Masha[0]*L;
            G[abs(SiteIra-SiteMasha)] += 1;
            if (Masha == Ira){
                Z += 1; 
                Ira = { rand() % L, rand() % L };
                Masha = {Ira[0],Ira[1]};
            }
            sweep = WA_Sweep(L,beta,Nb,Masha,BondsX,BondsY);
            Masha = sweep.Masha; Nb = sweep.Nb; BondsX = sweep.BondsX; BondsY = sweep.BondsY;
        }
    }
    //----Computing averages and errors----/
    Chi = mean(CHI)*beta;
    dChi = Jackknife(CHI,factores(Nmeas))*beta;
    Cv = (mean(NB2_NB) - mean(NB)*mean(NB))/(L*L);
    dCv =  Jackknife(NB2_NB,factores(Nmeas)) + 2*mean(NB)*Jackknife(NB,factores(Nmeas)) / (L*L);
    results.Chi = Chi; results.dChi = dChi; results.Cv = Cv; results.dCv = dCv;
    return results;
}

int main() {
    srand(time(0));
    int L, Ntherm, Nmeas, Nsteps, Nbeta;
    double beta_min, beta_max; 
    //***Input data***//
    std::cout << "L: ";
    std::cin >> L;
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

    RESULTS results;
    double Chi, Cv, dChi, dCv;
    std::vector<double> Betas = linspace(beta_min,beta_max, Nbeta);
    

    char NameCv[50], Cv_str[50], NameChi[50], Chi_str[50];
    sprintf(NameCv, "CvL%d_Meas%d.txt", L, Nmeas);
    sprintf(NameChi, "ChiL%d_Meas%d.txt", L, Nmeas);
    std::ofstream Cvfile; std::ofstream Chifile;
    Cvfile.open(NameCv); Chifile.open(NameChi);

    for (int i = 0; i < Nbeta; i++) {
        std::cout << "beta = " << Betas[i] << "  T = " << 1/Betas[i] << std::endl;
        results = WA_IsingModel(Betas[i], L, Ntherm, Nmeas, Nsteps);
        Chi = results.Chi; dChi = results.dChi; Cv = results.Cv; dCv = results.dCv;
        sprintf(Cv_str, "%-30.17g%-30.17g%-30.17g\n", Betas[i], Cv, dCv);
        sprintf(Chi_str, "%-30.17g%-30.17g%-30.17g\n", Betas[i], Chi, dChi);
        Cvfile << Cv_str;
        Chifile << Chi_str;
        
        std::cout << "Cv = " << Cv << " +- " << dCv << std::endl;
        std::cout << "Chi = " << Chi << " +- " << dChi << std::endl;
        std::cout << "------------------------------" << std::endl;
    }
    Cvfile.close();
    Chifile.close();

    return 0;      
}
   

