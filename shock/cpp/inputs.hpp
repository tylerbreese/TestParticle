#ifndef INPUTS_HPP
#define INPUTS_HPP

#include <armadillo>
#include <cmath>

namespace inputs {
    using namespace arma;

    // --- Physical Constants ---
    const double c  = 3.0E10;          // speed of light, cm/s (cgs)
    const double e  = 4.8E-10;         // elementary charge, cgs
    const double Z  = 1.0;             // atomic number
    const double A  = 1.0;             // mass number
    const double q  = 1.0 * e;         // ion charge 
    const double m  = A * 1.67E-24;    // mass in grams
    const double PI = 3.14159265358979323846;

    // --- Shock Inputs ---
    const double n_SW = 5.0; // cm^-3 SW thermal density
    const double n_PUI = 0.3; // cm^-3 PUI density
    const double R0 = 62.0; // in AU
    const double th  = 60.0 * (PI / 180.0); 
    const double del = 0.0 * (PI / 180.0);
    // --- Simulation Inputs ---
    const int sample_size = 50;
    const double U0  = 800.0e5;
    const double B0  = 0.05e-5;
    const double Usw = 400.0e5;
    const double cycles = 10000.0;
    const double var = 0.5;
    inline mat init_U0() {
        return {U0 * cos(del), 0.0, U0 * sin(del)};
    }
    inline mat init_B0() {
        return {B0 * cos(th), 0.0, B0 * sin(th)};
    }
    inline mat init_Usw() {
        return {Usw * cos(del), 0.0, Usw * sin(del)};
    }

    // Gyrofrequency
    //inline double get_B1() { return as_scalar(norm(init_B0(), 2)); }
    inline double get_B1() { return B0; }
    inline double get_Om() { return (abs(q) * get_B1()) / (m * c); }
    inline double get_Rg() { return as_scalar(norm(init_Usw(), 2)) / get_Om(); }
}

#endif


// old method just drop this into main.exe
    // // --- Inputs ---
    // const int sample_size = 10;
    // const double th  = 90.0 * (3.14/180.0); // needs to be radians
    // const double del = 0.0 * (3.14/180.0);
    // mat U0(1,3); //shock speed (U1)
    // U0.col(0) = 1000.0e5 * cos(del);
    // U0.col(1) = 0.0;
    // U0.col(2) = 1000.0e5 * sin(del);
    // mat B0(1,3); //magnetic field strength Gauss
    // B0.col(0) = 5.0e-5 * cos(th);
    // B0.col(1) = 0.0;
    // B0.col(2) = 5.0e-5 * sin(th);
    // mat Usw(1,3); // Solar Wind Speed
    // Usw.col(0) = 400.0e5 * cos(del);
    // Usw.col(1) = 0.0;
    // Usw.col(2) = 400.0e5 * sin(del);

    // const double B1 = as_scalar(vecnorm(B0,2,1)); // field strength
    // const double s2 = 0.5 * pow(B1,2); // wave variance
    // const double Om = q*B1 / (m*c); // gyrofreq
    // const double Rg = as_scalar(vecnorm(Usw,2,1)) / Om; // gyro radii

    // const double cycles = 1; //number of time segments