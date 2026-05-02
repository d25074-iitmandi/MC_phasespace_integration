#include "physics.hpp"
#include <cmath>

const double alpha = 1.0/137.0;
const double e2 = 4.0 * M_PI * alpha;
const double mu_mass = 0.1056; //Muon mass in GeV 

// Matrix element and cross-section changed to accomodate com energy closer to muon mass
double matrix_element(double costheta, double s_sqrt){ //Matrix element changed for energy closer to muon mass
	double s = s_sqrt*s_sqrt;
	double beta = 4.0 * pow(mu_mass,2)/s;
	return e2 * e2 * 0.5*(((1.0 + beta) + (1.0 - beta)*costheta*costheta));
}
 
double analytic_cross_section(double s_sqrt){
	double s = s_sqrt*s_sqrt;
	double beta = 4.0 * pow(mu_mass,2)/s;
	return ((4.0*M_PI*alpha*alpha)/(3.0*s)) * sqrt((1.0 - beta))*((1.0 + 0.5*beta));
}
