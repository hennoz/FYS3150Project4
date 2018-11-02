#include <iostream>
#include <random>
#include <fstream>
#include <iomanip>

using namespace std;

inline int PeriodicBoundary( int i, int limit, int add ) {
    return ( i + limit + add ) % ( limit );
}
//  Initial E, M, and SpinMatrix
//  && => R-value reference
void InitializeLattice ( int dim, double **SpinMatrix, double &E, double &M ) {
    for ( int x = 0; x < dim; x++ ) {
        for ( int y = 0; y < dim; y++ ) {
            SpinMatrix[x][y] = 1.0;
            M += ( double ) SpinMatrix[x][y];
        }
    }
//    cout << "M = " << M << endl;
    for ( int x = 0; x < dim; x++ ) {
        for ( int y = 0; y < dim; y++ ) {
            E -= ( double ) SpinMatrix[y][x]*
                    ( SpinMatrix[PeriodicBoundary( y, dim, -1 )][x] +
                    SpinMatrix[y][PeriodicBoundary( x, dim, -1 )] );
        }
    }
//    cout << "E = " << E << endl;
}

void MetropolisSampling ( int dim, int MCcycles, double T ) {
    //  Random number
    random_device rd;
    mt19937_64 gen( rd() );
    //  Uniform "REAL" distribution for x \in [0,1]
    uniform_real_distribution<double>RandomNumberGenerator( 0, 1 );
    //  Initialize lattice spin values
    uniform_int_distribution<int>SpinDistribution( 0, dim-1 );

    double **SpinMatrix = new double*[dim];
    for ( int i = 0; i < dim; i++ ) SpinMatrix[i] = new double[dim];
    double *ExpectVal = new double[5];
    double E, M = 0;
    for ( int i = 0; i < 5; i++ ) ExpectVal[i] = 0;

    InitializeLattice ( dim, SpinMatrix, E, M );
    double *EnergyDifference = new double[17];

    for ( int i = 0; i < 17; i += 4 ) EnergyDifference[i] = exp( -( i-8 )/T );

    //  Start Monte Carlo cycle
    for ( int cycle = 1; cycle <= MCcycles; cycle++ ) {
        for ( int x = 0; x < dim; x++ ) {
            for ( int y = 0; y < dim; y++ ) {
                int ix = SpinDistribution( gen );
//                cout << "ix = " << ix << endl;
                int iy = SpinDistribution( gen );
//                cout << "iy = " << iy << endl;
                int DeltaE = 2*SpinMatrix[ix][iy]*
                        ( SpinMatrix[ix][PeriodicBoundary( iy, dim, -1 )] +
                        SpinMatrix[PeriodicBoundary( ix, dim, -1 )][iy] +
                        SpinMatrix[ix][PeriodicBoundary( iy, dim, 1 )] +
                        SpinMatrix[PeriodicBoundary( ix, dim, 1 )][iy]);
                if ( RandomNumberGenerator( gen ) <= EnergyDifference[DeltaE + 8] ) {
                    SpinMatrix[ix][iy] += -1.0; // Flip one spin => new configuration
                    M += (double) 2*SpinMatrix[ix][iy];
                    E += (double) DeltaE;
                }
            }
        }
        ExpectVal[0] += E;
        ExpectVal[1] += E*E;
        ExpectVal[2] += M;
        ExpectVal[3] += M*M;
        ExpectVal[4] += fabs( M );


    }
    double E_expectVal  = ExpectVal[0]/(( double ) MCcycles );
    double E2_expectVal = ExpectVal[1]/(( double ) MCcycles );
    double M_expectVal  = ExpectVal[2]/(( double ) MCcycles );
    double M2_expectVal = ExpectVal[3]/(( double ) MCcycles );
    double M_absolute   = ExpectVal[4]/(( double ) MCcycles );
    double E_variance = ( E2_expectVal - E_expectVal*E_expectVal)/dim/dim;
    double M_variance = ( M2_expectVal - M_absolute*M_absolute )/dim/dim;

    cout << "T = " << T << endl;
    cout << "Number of cycles = " << MCcycles << endl;
    cout << "<E> = " << E_expectVal << endl;
    cout << "|M| = " << M_absolute << endl;
    cout << "Susceptibility = " << M_variance << endl;
    cout << "C_V = " << E_variance << endl;
}


int main()
{
    int dim = 2;
    int MCcycles = 100;
    double T = 1.0;
    MetropolisSampling ( dim, MCcycles, T );
    return 0;
}

//void WriteResultstoFile( int NSpins, int MCcycles, double T, double *ExpectVal ) {
//    double E_expectVal  = ExpectVal[0]/(( double ) MCcycles );
//    double E2_expectVal = ExpectVal[1]/(( double ) MCcycles );
//    double M_expectVal  = ExpectVal[2]/(( double ) MCcycles );
//    double M2_expectVal = ExpectVal[3]/(( double ) MCcycles );
//    double M_absolute   = ExpectVal[4]/(( double ) MCcycles );
//    // all expectation values are per spin, divide by 1/NSpins/NSpins
//    double E_variance = ( E2_expectVal - E_expectVal*E_expectVal)/NSpins/NSpins;
//    double M_variance = ( M2_expectVal - M_absolute*M_absolute )/NSpins/NSpins;
//    ofstream ofile;
//    ofile << setw(15) << setprecision(8) << T;
//    ofile << setw(15) << setprecision(8) << E_expectVal/NSpins/NSpins;
//    ofile << setw(15) << setprecision(8) << E_variance/T/T;
//    ofile << setw(15) << setprecision(8) << M_expectVal/NSpins/NSpins;
//    ofile << setw(15) << setprecision(8) << M_variance/T/T;
//    ofile << setw(15) << setprecision(8) << M_absolute/NSpins/NSpins << endl;
//} // end output function
