#ifndef INITIALIZELATTICE_H
#define INITIALIZELATTICE_H

inline int PeriodicBoundary( int i, int limit, int add );
void InitializeLattice ( int dim, double **SpinMatrix, double &E, double &M, int ordered );

#endif // INITIALIZELATTICE_H
