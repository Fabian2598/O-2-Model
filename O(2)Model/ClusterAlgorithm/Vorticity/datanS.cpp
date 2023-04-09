#include <iostream>
#include <vector> 
#include <fstream>
#include <cmath> 
#include <algorithm>

//inFile.is_open() is True if inFile is open. However, when we extract the data of inFile, 
//if there is a white space or a line break, the method is_open() closes the file. To prevent
//this we must use a while loop with the inFile.good() method

//----------.dat files structure for the cluster algorithm-------------//
//For the vorticity the structure is
// Ntherm, Nmeas, Nstep, L, Beta
// StringDensity, dStringDensity
//StringLength, dStringLength
// Vortex, dVortex
// Antivortex, dAntivortex
// computing time (in seconds)
//------------------------------------------------------------------//

int L, Nmeas, Ntherm, Nstep;
double beta, rhoS, rhoLenght, t, vort, avort; 
double drhoS, drhoLength, dvort, davort;

int main(){
    char Files[100];
    char name[100];
    std::cout << "File with list of names: ";
    std::cin >> Files;
    std::cout << " " << std::endl;

    double values;
    std::vector<std::string> Names;
    std::ifstream inFile; 
    inFile.open(Files);
    int cont = 0; 
    if (inFile.is_open() ) {  
        while ( inFile.good() ) {
        inFile >> name; 
        Names.push_back(name);    
        cont += 1;
        } 
    } 
    inFile.close();

    std::ofstream RhoSfile; std::ofstream RhoLenfile;  std::ofstream Timefile; std::ofstream Vortfile; std::ofstream Avortfile;

    for(int i = 0; i<cont-1; i++){
        std::vector<double> DATA; 
        //We extract the data from the .dat files//
        std::ifstream Data;
        Data.open(Names[i]);
        if (Data.is_open() ) {  
            while ( Data.good() ) {
            Data >> values; 
            DATA.push_back(values);    
            } 
        } 
        Data.close(); 

        Ntherm = (int) DATA[0]; Nmeas = (int) DATA[1]; Nstep = (int) DATA[2]; L = (int) DATA[3];
        beta = DATA[4]; rhoS = DATA[5]; drhoS = DATA[6]; rhoLenght = DATA[7]; drhoLength = DATA[8];
        vort = DATA[9]; dvort = DATA[10]; 
        avort = DATA[11]; davort = DATA[12]; t = DATA[13]; 
        char NameVort[100], Vort_str[100], NameRhoS[100], RhoS_str[100], NameTime[100], Time_str[100];
        char NameAvort[100], Avort_str[100], NameRhoSLen[100], RhoSLen_str[100];

        if (i == 0){
            sprintf(NameRhoS, "%dDXY_Sdensity_L%d_Meas%d.txt", 3,L, Nmeas);
            sprintf(NameRhoSLen, "%dDXY_SLen_L%d_Meas%d.txt", 3,L, Nmeas);
            sprintf(NameVort, "%dDXY_Vort_L%d_Meas%d.txt", 3, L, Nmeas);
            sprintf(NameAvort, "%dDXY_Avort_L%d_Meas%d.txt", 3, L, Nmeas);
            sprintf(NameTime, "%dDXY_Time_L%d_Meas%d.txt", 3, L, Nmeas);
            Vortfile.open(NameVort); Avortfile.open(NameAvort); RhoSfile.open(NameRhoS); RhoLenfile.open(NameRhoSLen); Timefile.open(NameTime);
        }
        sprintf(RhoS_str, "%-30.17g%-30.17g%-30.17g\n", beta, rhoS, drhoS);
        sprintf(RhoSLen_str, "%-30.17g%-30.17g%-30.17g\n", beta, rhoLenght, drhoLength);
        sprintf(Vort_str, "%-30.17g%-30.17g%-30.17g\n", beta, vort, dvort);
        sprintf(Avort_str, "%-30.17g%-30.17g%-30.17g\n", beta, avort, davort);
        sprintf(Time_str, "%-30.17g%-30.17g\n", beta, t);
        RhoSfile << RhoS_str;
        RhoLenfile << RhoSLen_str;
        Vortfile << Vort_str;
        Avortfile << Avort_str;
        Timefile << Time_str;
        std::cout << "------------------------------------------------------------------------" << std::endl;
        std::cout << "| beta = " << beta << "  T = " << 1/beta << "  L = " << L << std::endl;
        std::cout << "| Ntherm = " << Ntherm << "  Nmeas = " << Nmeas << "  Nstep = " << Nstep << std::endl;
        std::cout << "| rhoS = " << rhoS << " +- " << drhoS << std::endl;
        std::cout << "| String length = " << rhoLenght << " +- " << drhoLength << std::endl;
        std::cout << "| Vorticity = " << vort << " +- " << dvort << std::endl;
        std::cout << "| Antivorticity = " << avort << " +- " << davort << std::endl;
        std::cout << "| Time = " << t << " seconds" << std::endl; 
        std::cout << "------------------------------------------------------------------------" << std::endl; 
    }
    RhoSfile.close(); RhoLenfile.close(); Timefile.close(); Vortfile.close(); Avortfile.close();
    return 0;
}