#include <random>

double E_() {
    double Z = 4*cosh( 8 ) + 12;
    return -32/Z*sinh( 8 );
}
double E2_() {
    double Z = 4*cosh( 8 ) + 12;
    return 256/Z*cosh( 8 );
}
double M_() {
    return 0;
}
double M2_() {
    double Z = 4*cosh( 8 ) + 12;
    return 32/Z*(exp( 8 ) + 1 );
}
double xhi() {
    return M2_()/* - M_()*M_()*/;
}
double CV() {
    return E2_() - E_()*E_();
}
double absM() {
    double Z = 4*cosh( 8 ) + 12;
    return 8/Z*(exp( 8 ) + 2 );
}

