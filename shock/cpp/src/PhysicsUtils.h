#ifndef PHYSICSUTILS_H
#define PHYSICSUTILS_H

#include <Eigen/Dense>
#include <vector>
#include <tuple>

// --- Global Constants Struct (Replaces MATLAB's 'global') ---
struct Constants {
    // MATLAB used cgs/g units, so we stick to them here for consistency
    const double c = 3.0E10;   // speed of light, cm/s
    const double e = 4.8E-10;  // elementary charge, cgs
    const double Z = 1.0;      // atomic number
    const double A = 1.0;      // mass number
    const double q = Z * e;    // ion charge (cgs)
    const double m = A * 1.67E-24; // ion mass, g
};
extern const Constants CONSTS; // Declaration

// --- Type Aliases for Clarity ---
// A split particle set is X, V, and a flag vector
using ParticleSet = std::tuple<Eigen::MatrixXd, Eigen::MatrixXd, Eigen::VectorXi>;
// The main split structure is a vector of optional particle sets
using SplitContainer = std::vector<std::optional<ParticleSet>>;

// --- Function Declarations ---

// Replaces MATLAB's integrate_boris: returns {X_next, V_next}
std::pair<Eigen::MatrixXd, Eigen::MatrixXd> integrate_boris(
    const Eigen::MatrixXd& Xn, 
    const Eigen::MatrixXd& Vn, 
    const Eigen::MatrixXd& U, // Velocity field (N_particles x 3)
    const Eigen::MatrixXd& B, // Magnetic field (N_particles x 3)
    double dt
);

// Replaces MATLAB's particle_split
SplitContainer particle_split(
    const Eigen::MatrixXd& X_all, // X(:, :, n+1) 
    const Eigen::MatrixXd& V_all, // V(:, :, n+1)
    const Eigen::VectorXd& Vmag_all_start, // Vmag(:, 1)
    const Eigen::VectorXd& Vmag_all_current, // Vmag(:, n+1)
    const Eigen::VectorXd& cutoff, 
    SplitContainer& split_in_out // Passed by reference for update
);

// Replaces MATLAB's particle_displace
void particle_displace(
    SplitContainer& split_in_out, 
    double R_g
);

#endif // PHYSICSUTILS_H