#include <iostream>
#include <cmath>
#include <iomanip>
#include <Eigen/Dense>

// Include the class headers
#include "ShockFields.h" 
#include "InitParticles.h"

using namespace Eigen;
using namespace std;

// Function to print a 3D Eigen Vector
void print_vector(const std::string& name, const Eigen::Vector3d& vec) {
    cout << name << ": [" 
         << fixed << setprecision(4) << vec(0) << ", " 
         << vec(1) << ", " 
         << vec(2) << "]" << endl;
}

int main() {
    // --- 1. Define Initial Parameters ---
    // Using default values for simplicity, but you can customize them here.
    const int N_modes = 200;
    const int N_particles = 1000;
    const double theta_shock = M_PI / 2.0;
    const double delta_angle = 0.0;
    
    // Default Upstream B and U from ShockFields
    Vector3d B0_default(1.0 / std::sqrt(2.0), 0.0, 1.0 / std::sqrt(2.0));
    Vector3d U0_default(1.0 / std::sqrt(2.0), 0.0, 1.0 / std::sqrt(2.0));
    Vector3d x0_shock(0.0, 0.0, 0.0); // Shock at the origin

    // Particle location for field calculation
    Vector3d X_upstream(-1.0, 0.0, 0.0); // Sample point upstream of the shock
    Vector3d X_downstream(1.0, 0.0, 0.0); // Sample point downstream of the shock

    // --- 2. Instantiate Objects ---
    
    // Initialize ShockFields object (handles B and U calculation, and turbulence init)
    ShockFields SF_sim(B0_default, U0_default, x0_shock, N_modes);
    
    // Initialize InitParticles object (handles particle sampling and shock parameters)
    InitParticles IP_sim(U0_default, theta_shock, delta_angle, N_particles);

    cout << "--- Simulation Setup ---" << endl;
    cout << "Number of modes (N): " << SF_sim.N << endl;
    cout << "Number of particles: " << IP_sim.sample_size << endl;
    print_vector("Upstream Velocity (U0)", SF_sim.U0);
    print_vector("Upstream B-Field (B0)", SF_sim.B0);
    cout << "------------------------" << endl;

    // --- 3. Initialize Turbulence Spectrum ---
    
    cout << "1. Initializing Turbulence Spectrum..." << endl;
    
    // Call init_turby using properties from the SF object (dummy w_g and s used)
    std::pair<VectorXd, VectorXd> turb_init = SF_sim.init_turby(
        SF_sim.U0, SF_sim.w_g, SF_sim.s, SF_sim.N
    );
    
    SF_sim.An = turb_init.first;   // Store Amplitudes
    SF_sim.kk = turb_init.second;  // Store k-vectors

    cout << "   - Calculated " << SF_sim.An.size() << " turbulence amplitudes (An)." << endl;

    // --- 4. Initialize Shock Jump Parameters ---
    
    cout << "2. Calculating Shock Jump Parameters..." << endl;
    
    // Get the jump ratios (r, a, b) from the InitParticles class
    InitParticles::ShockInitParams P_shock = IP_sim.init_shock(
        IP_sim.U0, IP_sim.th, IP_sim.del
    );
    
    cout << "   - Density Compression Ratio (r): " << P_shock.r << endl;
    cout << "   - B-Field Parameter (a): " << P_shock.a << endl;
    cout << "   - Velocity Parameter (b): " << P_shock.b << endl;
    cout << "------------------------" << endl;

    // --- 5. Sample Initial Particle Velocities ---
    
    cout << "3. Sampling Initial Particle Velocities..." << endl;
    
    MatrixXd V0_particles = IP_sim.sampling(IP_sim.sample_size, IP_sim.U0);
    
    cout << "   - Sampled " << V0_particles.rows() << " particle 3D velocities." << endl;
    cout << "   - First particle V0: [" 
         << V0_particles(0, 0) << ", " 
         << V0_particles(0, 1) << ", " 
         << V0_particles(0, 2) << "]" << endl;
    cout << "------------------------" << endl;

    // --- 6. Calculate Fields at Specific Locations ---

    cout << "4. Calculating Fields..." << endl;
    
    // --- Location 1: Upstream ---
    
    std::pair<Vector3d, Vector3d> fields_up = SF_sim.shock_field(
        SF_sim.B0, SF_sim.U0, SF_sim.x0, 
        X_upstream, SF_sim.An, SF_sim.kk, P_shock, SF_sim.th, SF_sim.del
    );
    
    cout << "--- At Upstream Location (" << X_upstream(0) << ", " << X_upstream(1) << ", " << X_upstream(2) << ") ---" << endl;
    print_vector("Magnetic Field (B)", fields_up.first);
    print_vector("Velocity Field (U)", fields_up.second);
    
    // --- Location 2: Downstream ---
    
    std::pair<Vector3d, Vector3d> fields_down = SF_sim.shock_field(
        SF_sim.B0, SF_sim.U0, SF_sim.x0, 
        X_downstream, SF_sim.An, SF_sim.kk, P_shock, SF_sim.th, SF_sim.del
    );

    cout << "--- At Downstream Location (" << X_downstream(0) << ", " << X_downstream(1) << ", " << X_downstream(2) << ") ---" << endl;
    print_vector("Magnetic Field (B)", fields_down.first);
    print_vector("Velocity Field (U)", fields_down.second);
    
    cout << "------------------------" << endl;
    
    return 0;
}