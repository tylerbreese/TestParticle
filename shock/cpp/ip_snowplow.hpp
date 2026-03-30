#ifndef IP_SNOWPLOW_HPP
#define IP_SNOWPLOW_HPP
#include <armadillo>

// main function
void ip_snowplow(double Us, double V_SW, const arma::vec& T, double dt, arma::vec& U, arma::vec& Cs, arma::vec& vA, arma::vec& R);

// helpers
double get_density(double r);
double get_dMdt(double rho, double R, double U_relative);
double get_accel(double dMdt, double M, double U_relative);
#endif