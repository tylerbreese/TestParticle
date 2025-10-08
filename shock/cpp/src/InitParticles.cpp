#include "InitParticles.h"
#include <iostream>
#include <random>

using namespace Eigen;
using namespace std;

// --- Helper Functions ---

// Equivalent of MATLAB's trapz(x, y)
double InitParticles::trapz_eigen(const VectorXd& x, const VectorXd& y) const {
    if (x.size() != y.size() || x.size() < 2) {
        return 0.0;
    }
    double integral = 0.0;
    for (int i = 0; i < x.size() - 1; ++i) {
        // Trapezoid area = 0.5 * (y[i] + y[i+1]) * (x[i+1] - x[i])
        integral += 0.5 * (y(i) + y(i + 1)) * (x(i + 1) - x(i));
    }
    return integral;
}

// Simplified linear interpolation for pdfrnd (CDF inverse lookup)
double InitParticles::interp1_linear(const VectorXd& X, const VectorXd& V, double xq, double fill_val) const {
    if (xq <= X(0)) return fill_val; 
    if (xq >= X(X.size() - 1)) return V(V.size() - 1);

    // Binary search for the index of the interval [X(i), X(i+1)] containing xq
    int i = 0;
    while (X(i + 1) < xq) {
        i++;
    }

    // Linear interpolation formula: V = V[i] + (V[i+1] - V[i]) * (xq - X[i]) / (X[i+1] - X[i])
    double x0 = X(i);
    double x1 = X(i + 1);
    double y0 = V(i);
    double y1 = V(i + 1);

    if (x1 == x0) return y0; // Should not happen with typical CDFs

    return y0 + (y1 - y0) * (xq - x0) / (x1 - x0);
}

// --- Constructor ---
InitParticles::InitParticles(
    const Vector3d& U0_, double th_, double del_, int sample_size_
) : U0(U0_), th(th_), del(del_), sample_size(sample_size_) 
{
    // The C++ constructor handles initial parameter setting.
}

// --- init_shock method ---
InitParticles::ShockInitParams InitParticles::init_shock(
    const Vector3d& U0_in, double th_in, double del_in
) {
    const double Va = 30e5; // cm/s (Alfvén velocity)
    const double Vs = 70e5; // cm/s (Sound velocity)
    
    // Calculate Mach numbers (using vecnorm)
    double U_norm = vecnorm_3D(U0_in);
    double Ms = U_norm / Vs; // Sound Mach number
    double Ma = U_norm / Va; // Alfvén Mach number
    
    // Adiabatic index assumed to be gamma = 5/3
    double r_val = (5.0/3.0 + 1.0) / (5.0/3.0 - 1.0 + (2.0 / (Ms * Ms)));
    
    // Note: MATLAB's ^ operator is std::pow in C++
    double Ma2 = Ma * Ma;
    double cos_th = std::cos(th_in);
    double sin_th = std::sin(th_in);
    double cos_del = std::cos(del_in);
    double cos_del2 = cos_del * cos_del;
    double cos_th2 = cos_th * cos_th;

    double denominator = Ma2 * cos_del2 - r_val * cos_th2;

    double a_val = r_val * (Ma2 * cos_del2 - cos_th2) / denominator;
    double b_val = ((r_val - 1.0) * cos_th * sin_th) / denominator;
    
    return {r_val, a_val, b_val};
}

// --- pdfrnd method ---
VectorXd InitParticles::pdfrnd(
    const VectorXd& x, const VectorXd& px, int sample_size_in
) {
    // 1. Calculate CDF
    // cdf = cumsum(px)
    VectorXd cdf = px.array().cumsum();
    
    // 2. Normalize (Divide by last element for max=1.0)
    // cdf = cdf / max(cdf)
    double max_cdf = cdf(cdf.size() - 1);
    cdf = cdf.array() / max_cdf; 

    // 3. Generate random numbers [0, 1]
    // rnd = rand(sample_size, 1)
    static std::default_random_engine generator;
    std::uniform_real_distribution<double> dist_0_1(0.0, 1.0);
    VectorXd rnd(sample_size_in);
    for (int i = 0; i < sample_size_in; ++i) {
        rnd(i) = dist_0_1(generator);
    }
    
    // 4. Interpolate (Inverse CDF lookup)
    // X = interp1(cdf, x, rnd, 'linear', 0)
    VectorXd X(sample_size_in);
    // MATLAB's 'linear', 0 means linear interpolation, and points outside the range 
    // are filled with 0. The inverse CDF lookup usually means the fill value 
    // for rnd=0 is x(1) (or x(0) in C++), and for rnd=1 is x(end). 
    // Assuming the MATLAB code intended to start X at x(0) when cdf is low.
    for (int i = 0; i < sample_size_in; ++i) {
        X(i) = interp1_linear(cdf, x, rnd(i), x(0));
    }
    
    return X;
}

// --- sampling method ---
MatrixXd InitParticles::sampling(int sample_size_in, const Vector3d& U0_in) {
    // 1. Define velocity range
    // v = linspace(1000, vecnorm(U0,2,1), 100)
    double U_norm = vecnorm_3D(U0_in);
    int v_size = 100;
    VectorXd v = VectorXd::LinSpaced(v_size, 1000.0, U_norm);

    // 2. Sample random angles (isotropic)
    // u = obj.pdfrnd(linspace(-1,1,2001), ones(1,2001), sample_size)
    // phi = obj.pdfrnd(linspace(0,2*pi,2001), ones(1,2001), sample_size)
    
    int angle_space_size = 2001;
    VectorXd u_space = VectorXd::LinSpaced(angle_space_size, -1.0, 1.0);
    VectorXd phi_space = VectorXd::LinSpaced(angle_space_size, 0.0, 2.0 * M_PI);
    VectorXd one_vec = VectorXd::Ones(angle_space_size); // Uniform PDF

    VectorXd u = pdfrnd(u_space, one_vec, sample_size_in);
    VectorXd phi = pdfrnd(phi_space, one_vec, sample_size_in);

    // 3. Calculate velocity distribution (Power-law in momentum, f(p) ~ p^1)
    // F = v.^power1; 
    int power1 = 1;
    VectorXd F = v.array().pow(power1);

    // 4. Normalize the PDF: F = F ./ trapz(v,F)
    double integral = trapz_eigen(v, F);
    F = F.array() / integral;
    
    // 5. Sample plasma velocity magnitude
    // V_plasma = obj.pdfrnd(v,F,sample_size)
    VectorXd V_plasma = pdfrnd(v, F, sample_size_in);

    // 6. Convert spherical (V_plasma, u, phi) to Cartesian (vx, vy, vz)
    // v = V_plasma, u = cos(theta), phi = phi angle
    // Note: Eigen's .array() is essential for element-wise operations (.*)
    
    // sqrt(1 - u.^2)
    VectorXd sin_theta = (1.0 - u.array().square()).sqrt();
    
    // vx = V_plasma .* sqrt(1 - u.^2) .* cos(phi);
    VectorXd vx = V_plasma.array() * sin_theta.array() * phi.array().cos();
    
    // vy = V_plasma .* sqrt(1 - u.^2) .* sin(phi);
    VectorXd vy = V_plasma.array() * sin_theta.array() * phi.array().sin();
    
    // vz = V_plasma .* u;
    VectorXd vz = V_plasma.array() * u.array();

    // 7. Apply Bulk Velocity Shift
    // V0 = [vx(:) + U0(1), vy(:) + U0(2), vz(:) + U0(3)]
    MatrixXd V0(sample_size_in, 3);
    
    // The .array() + scalar operation performs broadcasting in Eigen
    V0.col(0) = vx.array() + U0_in(0); 
    V0.col(1) = vy.array() + U0_in(1);
    V0.col(2) = vz.array() + U0_in(2);

    return V0; // Returns an N x 3 matrix
}