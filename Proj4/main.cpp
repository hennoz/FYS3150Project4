#include <iostream>
#include <random>
#include <fstream>
#include <iomanip>
#include "analytical.h"
#include "metropolissampling.h"
#include <mpi.h>

using namespace std;
ofstream outfile;

void output( int dim, int MCcycles, double T, double *ExpectVal ) {
    double norm = 1/(( double ) MCcycles );
    double meanE  = ExpectVal[0]*norm;
    double meanE2 = ExpectVal[1]*norm;
    double meanM  = ExpectVal[2]*norm;
    double meanM2 = ExpectVal[3]*norm;
    double absM   = ExpectVal[4]*norm;
    double varE   = ( meanE2 - meanE*meanE )/( dim*dim );
    double varM   = ( meanM2 - meanM*meanM )/( dim*dim );

    outfile.open("info.txt");
    outfile << setw(15) << setprecision(8) << T;
    outfile << setw(15) << setprecision(8) << meanE;
    outfile << setw(15) << setprecision(8) << meanE2;
    outfile << setw(15) << setprecision(8) << meanM;
    outfile << setw(15) << setprecision(8) << meanM2;
    outfile << setw(15) << setprecision(8) << absM;
    outfile << setw(15) << setprecision(8) << varM/T;
    outfile << setw(15) << setprecision(8) << varE/T/T;
} // end output function


int main()
{
//     PROJECT 4b)
    int dim = 2;
    //  To get money results, 1e7 MCcycles is needed
    int MCcycles = 1e6;
    double T = 1.0;
    MetropolisSampling ( dim, MCcycles, T );

    cout << endl;
    cout << "<E> (analytic)             = " << E_() << endl;
    cout << "|M| (analytic)             = " << absM() << endl;
    cout << "Heat capacity (analytic)   = " << CV()/( dim*dim ) << endl;
    cout << "Susceptibility (analytic)  = " << xhi()/( dim*dim ) << endl;

//    //  Project 4c)
//    int dim = 20;
//    int MCcycles = 1e5;
//    double T = 2.4;
//    MetropolisSampling( dim, MCcycles, T );

    output( dim, MCcycles, T, ExpectVal);
    return 0;
}
