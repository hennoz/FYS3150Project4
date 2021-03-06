#include <iostream>
#include <random>
#include <fstream>
#include <iomanip>
#include "initializelattice.h"
using namespace std;

void MetropolisSampling ( int dim, int MCcycles, int loopStart, int loopStop, double T, double *ExpectVal,
                          int ordered, double &Accept ) {
    //  Random number
    random_device rd;
    mt19937_64 gen( rd() );
    //  Uniform "REAL" distribution for x \in [0,1]
    uniform_real_distribution<double>RandomNumberGenerator( 0, 1 );
    //  Initialize lattice spin values
    uniform_int_distribution<int>SpinDistribution( 0, dim-1 );

    double **SpinMatrix = new double*[dim];
    for ( int i = 0; i < dim; i++ ) SpinMatrix[i] = new double[dim];
    double E = 0;
    double M = 0;
    for ( int i = 0; i < 5; i++ ) ExpectVal[i] = 0;

    InitializeLattice ( dim, SpinMatrix, E, M, ordered );
    double *EnergyDifference = new double[17];

      //  4c -- 4e
    //  Setup energy and magnetization vector
    //  for calculations of acceptratio
//    vector<vector<int>> EMV;      //  Initiate Energy,Magnetization vector
//    EMV.resize(2);
//    EMV[0].push_back( (int) E );  //  Append first values
//    EMV[1].push_back( (int) M );
    //  --
    for ( int i = 0; i < 17; i += 4 ) EnergyDifference[i] = exp( -( i-8 )/T );

    //  Start Monte Carlo cycle
//    int numOfAccepts = 0;                           //  Initiate counter
    for ( int cycle = loopStart; cycle <= loopStop; cycle++ ) {
        for ( int x = 0; x < dim; x++ ) {           //  Sweep lattice in x-direction
            for ( int y = 0; y < dim; y++ ) {       //  Sweep lattice in y-direction
                int ix = SpinDistribution( gen );   //  Random int \in[0,1]
                int iy = SpinDistribution( gen );   //  Random int \in[0,1]
                /*Edge elements are 'ignored' using the next couple of lines,
                  this is quite demanding for large lattice.
                  An alternative way would've been 'mazing, cuz time is money,
                  sadly this was the best I could do*/
                int DeltaE = 2*SpinMatrix[ix][iy]*
                        ( SpinMatrix[ix][PeriodicBoundary( iy, dim, -1 )] +
                        SpinMatrix[PeriodicBoundary( ix, dim, -1 )][iy] +
                        SpinMatrix[ix][PeriodicBoundary( iy, dim, 1 )] +
                        SpinMatrix[PeriodicBoundary( ix, dim, 1 )][iy] );
                //  If r \leq w, accept new configuration, else keep old
                if ( RandomNumberGenerator( gen ) <= EnergyDifference[DeltaE + 8] ) {
                    SpinMatrix[ix][iy] *= -1.0;     // Flip one spin => new configuration
                    M += 2*SpinMatrix[ix][iy];      //  Update magnetization
                    E += DeltaE;                    //  Update energy
//                    numOfAccepts++;                 //  Count number of accepted configurations
                }
            }
        }
//        EMV[0].push_back( E );            //  "push_back" in c++ is like "append" in python
//        EMV[1].push_back( fabs( M ) );

        ExpectVal[0] += E;
        ExpectVal[1] += E*E;
        ExpectVal[2] += M;
        ExpectVal[3] += M*M;
        ExpectVal[4] += fabs( M );

    }
//    //  4d) write to binary file the energies to compute histogram for T=1.0 and T=2.4
    ofstream ofile;
//    ofile.open("/Users/hennoz/FYS3150Project4/energiesT24.bin", ofstream::binary);
//    ofile.write(reinterpret_cast<const char*> (EMV[0].data()), EMV[0].size()*sizeof(int));
//    ofile.close();

//    //  4c) write to binary file
//    for(int i = 0; i < 2; i++){
//        ofile.open("/Users/hennoz/FYS3150Project4/" + to_string(i) + "calibrate" +
//                   to_string(dim) + "Cycles" + to_string(MCcycles) + "Ordered" + to_string(ordered)
//                   + ".bin", ofstream::binary);
//        ofile.write(reinterpret_cast<const char*> (EMV[i].data()), EMV[i].size()*sizeof(int));
//        ofile.close();
//    }

//    Accept = numOfAccepts/double(loopStop);

//    double norm = 1/(( double ) MCcycles );
//    double meanE  = ExpectVal[0]*norm;
//    double meanE2 = ExpectVal[1]*norm;
//    double meanM  = ExpectVal[2]*norm;
//    double meanM2 = ExpectVal[3]*norm;
//    double absM   = ExpectVal[4]*norm;
//    double varE   = ( meanE2 - meanE*meanE )/( dim*dim );
//    double varM   = ( meanM2 - meanM*meanM )/( dim*dim );

//    cout << endl;
//    cout << "L                = " << dim << endl;
//    cout << "T                = " << T << endl;
//    cout << "Number of cycles = " << MCcycles << endl << endl;

//    cout << "Metropolis gives " << endl;
//    cout << "<E>              = " << meanE << endl;
//    cout << "|M|              = " << absM << endl;
//    cout << "Heat capacity    = " << varE << endl;
//    cout << "Susceptibility   = " << varM << endl;

//    ofstream outfile;
//    outfile.open("info.txt");
//    outfile << setw(15) << setprecision(8) << T;
//    outfile << setw(15) << setprecision(8) << meanE;
//    outfile << setw(15) << setprecision(8) << meanE2;
//    outfile << setw(15) << setprecision(8) << meanM;
//    outfile << setw(15) << setprecision(8) << meanM2;
//    outfile << setw(15) << setprecision(8) << absM;
//    outfile << setw(15) << setprecision(8) << varM/T;
//    outfile << setw(15) << setprecision(8) << varE/T/T;
//    outfile.close();

    //  Memory deallocation
    for ( int i = 0; i < dim; i++ ) {
        delete [] SpinMatrix[i];
    }
    delete [] SpinMatrix;
    delete [] EnergyDifference;
}
