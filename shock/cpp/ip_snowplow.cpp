#include "ip_snowplow.hpp"
#include "inputs.hpp"
#include <iostream>
#include <armadillo>
#include <cmath>
using namespace arma;
using namespace std;
// --- Helper Functions ---

double get_density(double r) {
    const double AU = 1.496e13;
    // const double n_SW = 5.0;
    // const double n_PUI = 0.3;
    const double n_SW = inputs::n_SW;
    const double n_PUI = inputs::n_PUI;
    const double m = 1.67e-24;
    double R = (inputs::R0) * AU;

    double rho_SW = (m * n_SW) * pow(r / R, -2.45);
    double rho_PUI = (m * n_PUI) * pow(r / R, -1.0);

    return rho_SW + rho_PUI;
}

double get_dMdt(double rho, double R, double U_relative) {
    return inputs::PI * rho * pow(R, 2) * U_relative;
}

double get_accel(double dMdt, double M, double U_relative) {
    return -(dMdt / M) * U_relative;
}

// Function to calculate Cs and vA for a vector of R
void calculate_acoustic(const vec& r, vec& Cs, vec& vA) {
    const double AU = 1.496e13;
    const double B0 = 5.0e-5;
    const double mp = 1.67e-24;
    const double T_const = 0.5 * 11604;
    const double kb = 1.381e-16;
    double R = inputs::R0 * AU;

    Cs.set_size(r.n_elem);
    vA.set_size(r.n_elem);

    for (uword i = 0; i < r.n_elem; ++i) {
        double B = B0 * pow(r(i) / R, -2.0);
        double rho = get_density(r(i));
        double n = rho / mp;
        double Pressure = n * (kb * T_const);

        Cs(i) = sqrt((5.0 / 3.0) * Pressure / rho);
        vA(i) = B / sqrt(4.0 * inputs::PI * rho);
    }
}

// --- Main Snowplow Function ---

void ip_snowplow(double Us, double V_SW, const vec& T, double dt, vec& U, vec& Cs, vec& vA, vec& R) {
    const double AU = 1.496e13;
    int n_steps = T.n_elem;

    R = zeros<vec>(n_steps);
    U = zeros<vec>(n_steps);
    vec M = zeros<vec>(n_steps);

    R(0) = inputs::R0 * AU;
    M(0) = 6e16; // mass so far in grams, hard to find obs. comp.
    U(0) = Us;

    for (int ii = 0; ii < n_steps - 1; ++ii) {
        if (U(ii) < 0.0) break;

        double Uold = U(ii);
        double Mold = M(ii);
        double Rold = R(ii);

        // Halfstep 1
        double rho1 = get_density(Rold);
        double dMdt1 = get_dMdt(rho1, Rold, Uold - V_SW);
        
        double Mhalf = Mold + (dt / 2.0) * dMdt1;
        double Uhalf = Uold + (dt / 2.0) * get_accel(dMdt1, Mold, Uold - V_SW);
        double Rhalf = Rold + (dt / 2.0) * Uhalf;

        // Fullstep 2
        double rho2 = get_density(Rhalf);
        double dMdt2 = get_dMdt(rho2, Rhalf, Uhalf - V_SW);
        
        double Mfull = Mhalf + (dt / 2.0) * dMdt2;
        double Ufull = Uhalf + (dt / 2.0) * get_accel(dMdt2, Mfull, Uhalf - V_SW);
        double Rfull = Rhalf + dt * Ufull;

        // Save step
        M(ii + 1) = Mfull;
        R(ii + 1) = Rfull;
        U(ii + 1) = Ufull;
    }

    calculate_acoustic(R, Cs, vA);
    
}