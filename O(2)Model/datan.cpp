#include <iostream>
#include <vector> 
#include <fstream>
#include <cmath> 
#include <algorithm>

//inFile.is_open() is True if inFile is open. However, when we extract the data of inFile, 
//if there is a white space or a line break, the method is_open() closes the file. To prevent
//this we must use a while loop with the inFile.good() method

//----------.dat files structure-------------//
//For 2D XY model:   
// Ntherm, Nmeas, Nstep, L, Beta
// E, dE
// Cv, dCv
// chi, dchi
// rho, drho
// computing time (in seconds)
//-------------------------------------------//
//For 3D XY model:   
// Ntherm, Nmeas, Nstep, L, Beta
// E, dE
// Cv, dCv
// chi, dchi
// computing time (in seconds)
//-------------------------------------------//
int L, Nmeas, Ntherm, Nstep;
double beta, E, Cv, chi, rho, t; 
double dE, dCv, dchi, drho;

int main(){
    char Files[50];
    char name[50];
    int d;
    std::cout << "Dimension of the model (2 or 3): ";
    std::cin >> d;
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

    std::ofstream Cvfile; std::ofstream Chifile; std::ofstream Efile; std::ofstream Timefile;  std::ofstream Rhofile;

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
        beta = DATA[4]; E = DATA[5]; dE = DATA[6]; Cv = DATA[7]; dCv = DATA[8]; 
        chi = DATA[9]; dchi = DATA[10]; 
        if (d == 2){
            rho = DATA[11]; drho = DATA[12]; t = DATA[13];
        }
        else{
            t = DATA[11];
        }
        //------------------------------------//

        //Now we generate the files with all the information//
        char NameCv[50], Cv_str[50], NameChi[50], Chi_str[50], NameE[50], E_str[50], NameTime[50], Time_str[50];
        char NameRho[50], Rho_str[50];

        if (i == 0){
            sprintf(NameE, "%dDXY_E_L%d_Meas%d.txt", d,L, Nmeas);
            sprintf(NameCv, "%dDXY_Cv_L%d_Meas%d.txt", d, L, Nmeas);
            sprintf(NameChi, "%dDXY_Chi_L%d_Meas%d.txt", d, L, Nmeas);
            sprintf(NameTime, "%dDXY_Time_L%d_Meas%d.txt", d, L, Nmeas);
            Cvfile.open(NameCv); Chifile.open(NameChi); Efile.open(NameE); Timefile.open(NameTime); 
            if (d == 2){sprintf(NameRho, "2DXY_Rho_L%d_Meas%d.txt", L, Nmeas); Rhofile.open(NameRho);}
        }
        sprintf(E_str, "%-30.17g%-30.17g%-30.17g\n", beta, E, dE);
        sprintf(Cv_str, "%-30.17g%-30.17g%-30.17g\n", beta, Cv, dCv);
        sprintf(Chi_str, "%-30.17g%-30.17g%-30.17g\n", beta, chi, dchi);
        sprintf(Time_str, "%-30.17g%-30.17g\n", beta, t);
        if (d==2){sprintf(Rho_str, "%-30.17g%-30.17g%-30.17g\n", beta, rho, drho);Rhofile << Rho_str;}
        Efile << E_str;
        Cvfile << Cv_str;
        Chifile << Chi_str;
        Timefile << Time_str;
        std::cout << "------------------------------------------------------------------------" << std::endl;
        std::cout << "| beta = " << beta << "  T = " << 1/beta << "  L = " << L << std::endl;
        std::cout << "| Ntherm = " << Ntherm << "  Nmeas = " << Nmeas << "  Nstep = " << Nstep << std::endl;
        std::cout << "| E = " << E << " +- " << dE  << std::endl;
        std::cout << "| Cv = " << Cv << " +- " << dCv << std::endl;
        std::cout << "| Chi = " << chi << " +- " << dchi << std::endl; 
        if ( d==2 ){std::cout << "| Rho = " << rho << " +- " << drho << std::endl;}
        std::cout << "| Time = " << t << " seconds" << std::endl; 
        std::cout << "------------------------------------------------------------------------" << std::endl; 
    }

    Efile.close();
    Cvfile.close();
    Chifile.close();
    Timefile.close();
    if (d==2){Rhofile.close();}

    return 0;
}