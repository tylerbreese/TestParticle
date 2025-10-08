#ifndef SHOCKFIELDS_H
#define SHOCKFIELDS_H

#include <Eigen/Dense>
#include <vector>
#include <cmath>
#include <numeric>
#include <random>

using namespace Eigen;
using namespace std;
// Define a struct to hold the persistent turbulence data (equivalent to MATLAB's 'persistent turby')
struct TurbulenceData {
    VectorXd kx;
    VectorXd ky;
    VectorXd kz;
    VectorXd beta;
    VectorXd Anx;
    VectorXd Bnx;
    VectorXd Any;
    VectorXd Bny;
    VectorXd Anz;
    bool initialized = false;
};

class ShockFields {
public:
    // Properties (with public access for easy translation)
    Vector3d x0;       // Shock location [0, 0, 0]
    Vector3d B0;       // Upstream magnetic field vector
    Vector3d U0;       // Upstream velocity vector
    int N;                    // Number of modes (must be positive integer)
    double w_g;               // Gyrofrequency
    double s;                 // Variance s^2 / B^2
    VectorXd An;       // Turbulence amplitudes
    VectorXd kk;       // k-vectors
    Vector3d B;        // Magnetic field vector at current X
    Vector3d U;        // Velocity field vector at current X

    // Constructor
    ShockFields(
        const Vector3d& B0_ = Vector3d(1.0 / std::sqrt(2.0), 0.0, 1.0 / std::sqrt(2.0)),
        const Vector3d& U0_ = Vector3d(1.0 / std::sqrt(2.0), 0.0, 1.0 / std::sqrt(2.0)),
        const Vector3d& x0_ = Vector3d::Zero(),
        int N_ = 200
    );

    // Placeholder for P.init_shock since P is not defined in MATLAB code
    // In a full C++ project, P would be a separate class/struct
    struct ShockInitParams {
        double r; // Density ratio
        double a; // B field change param
        double b; // U field change param
    };
    ShockInitParams init_shock_placeholder(const Vector3d& U0_in, double th, double del);

    // Member methods
    std::pair<Vector3d, Vector3d> shock_field(
        const Vector3d& B0_in, const Vector3d& U0_in, const Vector3d& x0_in,
        const Vector3d& X, const VectorXd& An_in, const VectorXd& kk_in, 
        ShockInitParams P_in, double th, double del
    );

    std::pair<VectorXd, VectorXd> init_turby(
        const Vector3d& U_in, double w_g_in, double s_in, int N_in
    );

    std::pair<VectorXd, VectorXd> get_k_vec(
        double l_min, double l_max, int N_in
    );

    std::pair<Vector3d, Matrix3Xd> sum_turby(
        const Vector3d& X, int N_in, const VectorXd& kk_in, const VectorXd& An_in
    );

private:
    TurbulenceData turby_data; // Persistent data equivalent
    // Utility for MATLAB's vecnorm(U, 2, 1) for a single 3D vector
    double vecnorm_3D(const Vector3d& V) const {
        return V.norm(); 
    }
};

#endif // SHOCKFIELDS_H