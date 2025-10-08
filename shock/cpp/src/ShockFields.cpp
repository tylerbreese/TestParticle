#include "ShockFields.h"
#include <iostream>

using namespace Eigen;
using namespace std;

// --- Constructor ---
ShockFields::ShockFields(
    const Vector3d& B0_, const Vector3d& U0_, const Vector3d& x0_, int N_
) : x0(x0_), B0(B0_), U0(U0_), N(N_), w_g(1.0), s(0.5)
{
    // C++ constructor initializes properties; default values handled in header
    // Initialize N, w_g, s with defaults if not set by a later call, or use defaults provided.
    // The MATLAB validation {mustBeNumeric, ...} is replaced by C++ types and checks if needed.
}

// --- Placeholder for external shock initialization ---
ShockFields::ShockInitParams ShockFields::init_shock_placeholder(const Eigen::Vector3d& U0_in, double th, double del) {
    // This is a placeholder since the definition of P.init_shock is missing.
    // Replace with actual physics logic if available.
    return { 4.0, 1.5, 0.2 }; // Example values for r, a, b
}

// --- shock_field method ---
pair<Vector3d, Vector3d> ShockFields::shock_field(
    const Vector3d& B0_in, const Vector3d& U0_in, const Vector3d& x0_in,
    const Vector3d& X, const VectorXd& An_in, const VectorXd& kk_in,
    ShockInitParams P_in, double th, double del
) {
    // Note: The original MATLAB code passes B0, U0, x0, An, kk as arguments
    // but uses 'obj' properties inside other methods. Here, we use the arguments 
    // for local scope and assume the relevant physics object is P.

    ShockInitParams P = init_shock_placeholder(U0_in, th, del); // P.init_shock equivalent
    double r = P.r;
    double a = P.a;
    double b = P.b;

    // Call sum_turby to get the turbulence magnetic field perturbation (dB)
    pair<Vector3d, Matrix3Xd> turb_results = sum_turby(X, N, kk_in, An_in);
    Vector3d dB = turb_results.first;

    Vector3d U, B;

    // Shock taken in xz plane with downstream being the positive direction (X(1) < x0(1) is upstream)
    if (X(0) < x0_in(0)) { // Upstream
        U = U0_in;
        B = B0_in + dB;
    } else { // Downstream (After shock)
        // Note: MATLAB's '.*' for scalar/vector multiplication is implicit in Eigen/C++
        U(0) = U0_in(0) / r;
        U(1) = U0_in(1) + b * U0_in(0);
        U(2) = U0_in(2) + b * U0_in(0);

        B(0) = B0_in(0) + dB(0);
        B(1) = B0_in(1) + dB(1);
        B(2) = a * (B0_in(2) + dB(2));
    }

    return make_pair(B, U);
}

// --- init_turby method ---
pair<VectorXd, VectorXd> ShockFields::init_turby(
    const Vector3d& U_in, double w_g_in, double s_in, int N_in
) {
    // Calculate gyroradius r_gn = |U| / w_g(1)
    double r_gn = vecnorm_3D(U_in) / w_g_in; 

    // Constants
    const double L = 1.496e+13; // 1 AU in cm
    double l1 = 5.0 * L;
    double l2 = 0.5 * r_gn;

    // Get k-vectors
    pair<VectorXd, VectorXd> k_results = get_k_vec(l2, l1, N_in);
    VectorXd K = k_results.first;   // kk (k-vectors)
    VectorXd dK = k_results.second; // dkk (delta K)

    const double gamma = 11.0 / 3.0; // 3D turbulence
    
    // P0 = kk.^2 .* dkk (Element-wise multiplication)
    VectorXd P0 = K.array().square() * dK.array(); 

    // GkN = P0 ./ (1 + (kk.*L).^gamma) (Element-wise division)
    VectorXd GkN = P0.array() / (1.0 + (K.array() * L).pow(gamma));

    // Gk = sum(GkN, 2) -> sum(GkN) (Summing all elements as GkN is a 1D vector)
    double Gk = GkN.sum();

    // An = 2.0 * sqrt( s^2 ./ Gk) .* sqrt(GkN)
    VectorXd An = 2.0 * std::sqrt(s_in * s_in / Gk) * GkN.array().sqrt();
    
    return make_pair(An, K); // Return An and kk
}

// --- get_k_vec method ---
pair<VectorXd, VectorXd> ShockFields::get_k_vec(
    double l_min, double l_max, int N_in
) {
    double k_min = (2.0 * M_PI) / l_max;
    double k_max = (2.0 * M_PI) / l_min;

    // MATLAB's logspace(a, b, N+1) is equivalent to N+1 points logarithmically spaced
    // between 10^a and 10^b.
    double a = std::log10(k_min);
    double b = std::log10(k_max);

    // Create N+1 logarithmically spaced points
    VectorXd K_full(N_in + 1);
    for (int i = 0; i <= N_in; ++i) {
        double exponent = a + (b - a) * (double)i / N_in;
        K_full(i) = std::pow(10.0, exponent);
    }

    // dK = diff(K, 1, 2)
    VectorXd dK(N_in);
    for (int i = 0; i < N_in; ++i) {
        dK(i) = K_full(i + 1) - K_full(i);
    }

    // K = K(:, 2:end) (Remove the first element)
    VectorXd K = K_full.tail(N_in); 

    return make_pair(K, dK);
}

// --- sum_turby method ---
pair<Vector3d, Matrix3Xd> ShockFields::sum_turby(
    const Vector3d& X, int N_in, const VectorXd& kk_in, const VectorXd& An_in
) {
    if (!turby_data.initialized) {
        // Initialize persistent variables
        
        // Random number generation setup
        static std::default_random_engine generator;
        std::uniform_real_distribution<double> dist_0_2pi(0.0, 2.0 * M_PI);
        std::uniform_real_distribution<double> dist_n1_1(-1.0, 1.0);
        
        // beta = (2*pi).*rand(N,1)
        turby_data.beta.resize(N_in);
        for (int i = 0; i < N_in; ++i) {
            turby_data.beta(i) = dist_0_2pi(generator);
        }

        // u = 2.*rand(n,1) - 1 (Assuming 'n' is size of X(1), which is 1 for a single point)
        // Here, we assume the calculation is for a single point, so n=1, and u is a scalar.
        // HOWEVER, the MATLAB code structure suggests X could be a matrix of samples, 
        // and u/phi/alpha are n x 1. Since X is 3x1 here, we assume n=1 for simplicity.
        // If X were n x 3, we'd need to resize all related vectors.
        
        // Assuming n=1 (single sample) for X and all related vectors:
        double u = dist_n1_1(generator); 
        double phi = dist_0_2pi(generator);
        double alpha = dist_0_2pi(generator);

        // Calculate rotation matrix components
        double sin_u = std::sqrt(1.0 - u * u);
        
        // Wavelets (Anx, Any, Anz, Bnx, Bny) - These are N x 1 vectors
        // kk and An are N x 1 vectors, u/phi/alpha are scalars (n=1)
        turby_data.Anx = An_in.array() * std::cos(alpha) * std::cos(phi) * u;
        turby_data.Any = An_in.array() * std::cos(alpha) * std::sin(phi) * u;
        turby_data.Anz = -An_in.array() * std::cos(alpha) * sin_u;
        turby_data.Bnx = An_in.array() * std::sin(alpha) * std::sin(phi);
        turby_data.Bny = -An_in.array() * std::sin(alpha) * std::cos(phi);
        
        // kx, ky, kz - These are N x 1 vectors
        turby_data.kx = kk_in.array() * std::cos(phi) * sin_u;
        turby_data.ky = kk_in.array() * std::sin(phi) * sin_u;
        turby_data.kz = kk_in.array() * u;

        turby_data.initialized = true;
    }
    
    // arg = (turby.kx .* X(1) + turby.ky .* X(2) + turby.kz .* X(3) + turby.beta')
    // X is 3x1, kx/ky/kz are N x 1. This is an N x 1 vector.
    VectorXd arg = turby_data.kx.array() * X(0) + 
                   turby_data.ky.array() * X(1) + 
                   turby_data.kz.array() * X(2) + 
                   turby_data.beta.array(); // beta is already N x 1

    // dBn is N x 3 matrix in C++ (where N is number of modes, 3 is vector component)
    MatrixXd dBn_comp(N_in, 3); 

    // dBn(:,:,1) = turby.Anx .* cos(arg) + turby.Bnx .* sin(arg)
    dBn_comp.col(0) = turby_data.Anx.array() * arg.array().cos() + 
                      turby_data.Bnx.array() * arg.array().sin();

    // dBn(:,:,2) = turby.Any .* cos(arg) + turby.Bny .* sin(arg)
    dBn_comp.col(1) = turby_data.Any.array() * arg.array().cos() + 
                      turby_data.Bny.array() * arg.array().sin();
    
    // dBn(:,:,3) = turby.Anz .* cos(arg)
    // Note: The MATLAB code only has the cos term for the Z-component, check physics for Bnz * sin(arg)
    dBn_comp.col(2) = turby_data.Anz.array() * arg.array().cos(); 
    
    // dB = squeeze(sum(dBn,2)) -> Sum across the N-modes dimension, resulting in a 3x1 vector
    Vector3d dB;
    dB(0) = dBn_comp.col(0).sum();
    dB(1) = dBn_comp.col(1).sum();
    dB(2) = dBn_comp.col(2).sum();
    
    // Return dB as 3x1 and dBn_comp as 3xN (transposed for potential C++ use, or keep N*3)
    // Using a 3xN matrix to represent the MATLAB 3D array dBn (3 components, N modes)
    Matrix3Xd dBn_3xN = dBn_comp.transpose(); 

    return make_pair(dB, dBn_3xN);
}