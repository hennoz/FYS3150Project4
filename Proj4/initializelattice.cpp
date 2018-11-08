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
