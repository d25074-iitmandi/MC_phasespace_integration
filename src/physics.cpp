#include "physics.h"
#include <cmath>

const double alpha = 1.0/137.0;
const double e2 = 4.0 * M_PI * alpha;

double matrix_element(double costheta){
	return e2*e2*0.5*(1.0+costheta*costheta);
}
 
double analytic_cross_section(double s){
	return (4.0*M_PI*alpha*alpha)/(3.0*s);
}
