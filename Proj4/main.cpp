#include <iostream>
#include <random>
#include <fstream>
#include <iomanip>
#include "analytical.h"
#include "metropolissampling.h"
#include <mpi.h>
#include <chrono>

using namespace std;
ofstream outfile;

void output( int dim, double T, double *ExpectVal, int MCcycles, double timing );

int main( int argc, char *argv[] )
{
    string filename;
    int dim = 2;
    int MCcycles = 1e5;
    double InitialTemp = 2.0;
    double FinalTemp = 3.0;
    double TimeStep = 0.1;
    double *ExpectVal = new double[5];
    double *TotalExpectVal = new double[5];
    for( int i = 0; i < 5; i++ ) TotalExpectVal[i] = 0;
    double timing;
    chrono::high_resolution_clock::time_point t1;
    chrono::high_resolution_clock::time_point t2;

    int nProcs, my_rank;

    MPI_Init (&argc, &argv);
    MPI_Comm_size (MPI_COMM_WORLD, &nProcs);
    MPI_Comm_rank (MPI_COMM_WORLD, &my_rank);

    MPI_Bcast (&dim, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast (&InitialTemp, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast (&FinalTemp, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast (&TimeStep, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);

    int cycleInterval = MCcycles/nProcs;
    int loopStart = my_rank*cycleInterval;
    int loopStop = (my_rank+1)*cycleInterval;

    if (my_rank == 0){
        outfile.open("Lattice" + to_string(dim) + "Cycles" + to_string(MCcycles) + ".txt");
        outfile << setw(15) << setprecision(8) << "T";
        outfile << setw(15) << setprecision(8) << "E";
        outfile << setw(15) << setprecision(8) << "E2";
        outfile << setw(15) << setprecision(8) << "M";
        outfile << setw(15) << setprecision(8) << "M2";
        outfile << setw(15) << setprecision(8) << "M abs";
        outfile << setw(15) << setprecision(8) << "Susceptibility";
        outfile << setw(15) << setprecision(8) << "Heat capacity";
        outfile << setw(15) << setprecision(8) << "Run time";
        outfile << setw(15) << setprecision(8) << "No. of cycles" << endl;
    }
    for ( double T = InitialTemp; T <= FinalTemp; T+=TimeStep) {
        if (my_rank==0) t1 = chrono::high_resolution_clock::now();
        MetropolisSampling( dim, MCcycles, T, ExpectVal);
        for ( int i = 0; i < 5; i++ ) {
            MPI_Reduce(&ExpectVal[i], &TotalExpectVal[i], 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
        }
            if (my_rank==0) {
            t2 = chrono::high_resolution_clock::now();
            chrono::duration<double> time_span = std::chrono::duration_cast<chrono::duration<double>>(t2 - t1);
            timing = time_span.count();
            output(dim, T, TotalExpectVal, MCcycles, timing);
            cout << "T = " << T << " done...\n";
        }
    }
    outfile.close();

//------------------------------------------------------------------------
////     PROJECT 4b)
//    int dim = 2;
//    //  To get money results, 1e7 MCcycles is needed
//    int MCcycles = 1e6;
//    double T = 1.0;
//    MetropolisSampling ( dim, MCcycles, T, ExpectVal );
//    cout << endl;
//    cout << "<E> (analytic)             = " << E_() << endl;
//    cout << "|M| (analytic)             = " << absM() << endl;
//    cout << "Heat capacity (analytic)   = " << CV()/( dim*dim ) << endl;
//    cout << "Susceptibility (analytic)  = " << xhi()/( dim*dim ) << endl;
//--------------------------------------------------------------------------

    delete [] ExpectVal;
    delete [] TotalExpectVal;
    return 0;
}
void output( int dim, double T, double *ExpectVal, int MCcycles, double timing ) {
  for( int i = 0; i < 5; i++ ) ExpectVal[i] /= MCcycles;
  double E_variance = (ExpectVal[1] - ExpectVal[0]*ExpectVal[0])/dim/dim;
  double M_variance = (ExpectVal[3] - ExpectVal[2]*ExpectVal[2])/dim/dim;
  outfile << setw(15) << setprecision(8) << T;
  outfile << setw(15) << setprecision(8) << ExpectVal[0]/dim/dim;
  outfile << setw(15) << setprecision(8) << ExpectVal[1]/dim/dim;
  outfile << setw(15) << setprecision(8) << ExpectVal[2]/dim/dim;
  outfile << setw(15) << setprecision(8) << ExpectVal[3]/dim/dim;
  outfile << setw(15) << setprecision(8) << ExpectVal[4]/dim/dim;
  outfile << setw(15) << setprecision(8) << M_variance/T;
  outfile << setw(15) << setprecision(8) << E_variance/(T*T);
  outfile << setw(15) << setprecision(8) << timing;
  outfile << setw(15) << setprecision(8) << MCcycles << endl;
}
