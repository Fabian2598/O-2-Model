#include <iostream>
#include <vector> 
#include <fstream>
#include <cmath> 
#include <algorithm>

//inFile.is_open() is True if inFile is open. However, when we extract the data of inFile, 
//if there is a white space or a line break, the method is_open() closes the file. To prevent
//this we must use a while loop with the inFile.good() method

//-----------------------------.dat files----------------------------//
//For the vorticity the structure is
// Ntherm, tauQ, Nmeas, L, t  
// Vortex, dVortex
// computing time (in seconds)
//------------------------------------------------------------------//

int L, Nmeas, Ntherm, tauQ, t;
double vort, dvort, comp_time; 

int main(){
    char Files[100];
    char name[100];
    int mode, d;
    std::cout << "Dimension of the model (2 or 3): " << std::endl;
    std::cin >> d;
    std::cout << "Measuring time mode (1 --> t<2tauQ, 2 --> t=2tauQ) " << std::endl;
    std::cin >> mode;
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

    std::ofstream Timefile, Vortfile; 

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
        Ntherm = (int) DATA[0]; tauQ = (int) DATA[1]; Nmeas = (int) DATA[2]; L = (int) DATA[3];
        t = DATA[4]; vort = DATA[5]; dvort = DATA[6];
        comp_time = DATA[7]; //Computing time
        //------------------------------------//

        //Now we generate the files with all the information//
        char NameVort[100], Vort_str[100], NameTime[100], Time_str[100];
        if (i == 0){
            if (mode == 1){
                sprintf(NameVort, "%dDXY_Vort_L%d_Meas%d_t%d.txt", d, L, Nmeas,t);
                sprintf(NameTime, "%dDXY_Time_L%d_Meas%d_t%d.txt", d, L, Nmeas,t);
            }
            else if (mode == 2){
                sprintf(NameVort, "%dDXY_Vort_L%d_Meas%d_remainder.txt", d, L,Nmeas);
                sprintf(NameTime, "%dDXY_Time_L%d_Meas%d_remainder.txt", d, L,Nmeas);
            }
            Vortfile.open(NameVort); Timefile.open(NameTime); 
        }
        sprintf(Vort_str, "%-30d%-30.17g%-30.17g\n", tauQ, vort, dvort);
        sprintf(Time_str, "%-30d%-30.17g\n", tauQ, comp_time);
        Vortfile << Vort_str;
        Timefile << Time_str;
        std::cout << "------------------------------------------------------------------------" << std::endl;
        std::cout << "| tauQ = " << tauQ << "  t = " << t << "  L = " << L << std::endl;
        std::cout << "| Ntherm = " << Ntherm << "  Nmeas = " << Nmeas << std::endl;
        std::cout << "| Vort/V = " << vort << " +- " << dvort  << std::endl;
        std::cout << "| Time = " << comp_time << " seconds" << std::endl; 
        std::cout << "------------------------------------------------------------------------" << std::endl;   
        }     
    Vortfile.close();
    Timefile.close();
    return 0;
}