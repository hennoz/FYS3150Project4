#include <random>
using namespace std;

inline int PeriodicBoundary( int i, int limit, int add ) {
    return ( i + limit + add ) % ( limit );
}
//  Initial E, M, and SpinMatrix
//  && => R-value reference
void InitializeLattice ( int dim, double **SpinMatrix, double &E, double &M, int ordered ) {
    //  Random number
    random_device rd;
    mt19937_64 gen( rd() );
    //  Uniform "REAL" distribution for x \in [0,1]
    uniform_real_distribution<double>RandomNumberGenerator( 0, 1 );
    if ( ordered == 1 ) {
        for ( int x = 0; x < dim; x++ ) {
            for ( int y = 0; y < dim; y++ ) {
                SpinMatrix[x][y] = 1.0;
                M += ( double ) SpinMatrix[x][y];
            }
        }
    } else {
        for ( int x = 0; x < dim; x++ ) {
            for ( int y = 0; y < dim; y++ ) {
                SpinMatrix[x][y] = (RandomNumberGenerator( gen ) > 0.5) ? 1: -1;
                M += ( double ) SpinMatrix[x][y];
            }
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
