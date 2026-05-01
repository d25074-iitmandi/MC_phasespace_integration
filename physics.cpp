#include<"physics.h">
#include<cmath>

const double alpha = 1/137;
const double vac_perm = 8.85*e-12;
const double h = 6.626*e-34;
const int c = 299792458;
const double e2 = 2.0 * alpha * vac_perm * h * c;

double matrix_element(double costheta){
	return e2*e2*(1+costheta*costheta);
}

double analytic_cross_section(double s){
	return (4*M_PI*alpha*alpha)/(3*s);
}
