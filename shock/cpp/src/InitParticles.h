#ifndef INITPARTICLES_H
#define INITPARTICLES_H

#include <Eigen/Dense>
#include <vector>
#include <cmath>
#include <numeric>
#include <algorithm>

using namespace Eigen;
using namespace std;

class InitParticles {
public:
    // Properties
    Eigen::Vector3d U0;         // Upstream velocity [1/sqrt(2), 0, 1/sqrt(2)]
    int sample_size;            // Number of particles (1000)
    double del;                 // Delta angle (0.0)
    double th;                  // Theta angle (pi/2)

    // Struct to match the output of init_shock (same as in ShockFields)
    struct ShockInitParams {
        double r; // Density ratio
        double a; // B field change param
        double b; // U field change param
    };

    // Constructor
    InitParticles(
        const Eigen::Vector3d& U0_ = Eigen::Vector3d(1.0 / std::sqrt(2.0), 0.0, 1.0 / std::sqrt(2.0)),
        double th_ = M_PI / 2.0,
        double del_ = 0.0,
        int sample_size_ = 1000
    );

    // Method Declarations
    
    // Equivalent to init_shock
    ShockInitParams init_shock(const Eigen::Vector3d& U0_in, double th_in, double del_in);

    // Equivalent to sampling
    // Returns an Eigen::MatrixXd where each row is a particle's 3D velocity V0
    Eigen::MatrixXd sampling(int sample_size_in, const Eigen::Vector3d& U0_in);
    
    // Equivalent to pdfrnd (Probability Distribution Function Random)
    // C++ doesn't have a built-in interp1, so this will require manual implementation or a library
    Eigen::VectorXd pdfrnd(const Eigen::VectorXd& x, const Eigen::VectorXd& px, int sample_size_in);

private:
    // Utility for MATLAB's vecnorm(U, 2, 1) for a single 3D vector
    double vecnorm_3D(const Eigen::Vector3d& V) const {
        return V.norm(); 
    }

    // Helper function for MATLAB's trapz (Trapezoidal numerical integration)
    double trapz_eigen(const Eigen::VectorXd& x, const Eigen::VectorXd& y) const;

    // Helper for linear interpolation (used inside pdfrnd)
    double interp1_linear(const Eigen::VectorXd& X, const Eigen::VectorXd& V, double xq, double fill_val) const;
};

#endif // INITPARTICLES_H
