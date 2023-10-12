#ifndef STATISTICS_H_INCLUDED
#define STATISTICS_H_INCLUDED
#include <vector> 
#include <cmath>
#include <algorithm>

//mean of a vector
template <typename T>
double mean(std::vector<T> x){ 
    double prom = 0;
    for (T i : x) {
        prom += i*1.0;
    }   
    prom = prom / x.size();
    return prom;
}
//absolute value
template <typename T>
inline T absVal(T z){
    if (z < 0){
        return -z;
    }
    else{
        return z;
    }
}
//random double number in the inteval [a,b] a = min, b = max
inline double rand_range(double a, double b){
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
template <typename T>
std::vector<double> linspace(T min, T max, int n) {
    std::vector<double> linspace;
    double h = (1.0*max - 1.0*min) / (n - 1);
    for (int i = 0; i < n; ++i) {
        linspace.insert(linspace.begin() + i, min*1.0 + i * h); 
    }
    return linspace;
}
//---------------Logspace (similar to python)----------------------//
template <typename T>
std::vector<int> logspace(T min, T max, int n) {
    std::vector<int> logspace(n);
    double h = (max*1.0 - min*1.0) / (n - 1);
    for (int i = 0; i < n; ++i) {
        logspace[i] = (int) pow(10.0, min*1.0 + i * h); 
    }
    return logspace;
}

//n modulus m 
inline int modulo(int n, int m) {
    if (n < 0) {
        return (n + m) % m;
    }
    else {
        return n % m;
    }
}
//n modulus m 
inline double fmodulo(double n, double m) {
    if (n < 0) {
        return fmod(n + m, m);
    }
    else {
        return fmod(n,m);
    }
}

#endif