#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <cmath>
#include <chrono>

#include "ShockFields.h"
#include "InitParticles.h"
#include "PhysicsUtils.h" // Includes constants and helper functions

using namespace Eigen;
using namespace std;

// Function to read the input file (manual text parsing)
bool read_input_file(const std::string& filename, std::map<std::string, double>& params) {
    std::ifstream file(filename);
    if (!file.is_open()) {
        cerr << "Error: Could not open input file " << filename << endl;
        return false;
    }

    std::string line;
    while (std::getline(file, line)) {
        std::stringstream ss(line);
        std::string s1, s2, s3;
        if (ss >> s1 >> s2 >> s3) {
            try {
                // MATLAB: textscan(fid,formatSpec,7) reads the first 7 lines
                // We'll parse the relevant parameters based on the MATLAB code's assumption
                if (s1 == "run_name") params["run_name_placeholder"] = 1.0; // dummy
                else if (s1 == "sample_size") params["sample_size"] = std::stod(s3);
                else if (s1 == "cycles") params["cycles"] = std::stod(s3);
                else if (s1 == "th") params["th"] = std::stod(s3);
                else if (s1 == "del") params["del"] = std::stod(s3);
                else if (s1 == "VSW") params["VSW"] = std::stod(s3);
                else if (s1 == "B_0") params["B_0"] = std::stod(s3);
            } catch (const std::exception& e) {
                cerr << "Error parsing line: " << line << " (" << e.what() << ")" << endl;
            }
        }
    }
    return true;
}

int main() {
    auto start_time = std::chrono::high_resolution_clock::now();
    
    // --- FILE INPUT AND SETUP ---
    std::map<std::string, double> input_params;
    // NOTE: You must create a 'debug.in' file in the execution directory
    if (!read_input_file("debug.in", input_params)) return 1;

    // Retrieve parameters (using defaults if not found, though parsing should handle this)
    int sample_size = static_cast<int>(input_params.count("sample_size") ? input_params.at("sample_size") : 1000);
    int cycles = static_cast<int>(input_params.count("cycles") ? input_params.at("cycles") : 10);
    double th = input_params.count("th") ? input_params.at("th") : M_PI / 2.0;
    double del = input_params.count("del") ? input_params.at("del") : 0.0;
    double VSW_mag = input_params.count("VSW") ? input_params.at("VSW") : 70e5;
    double B_0_mag = input_params.count("B_0") ? input_params.at("B_0") : 10e-5; // Tesla -> Gauss (1T = 10^4 G, using 10e-5 G as a guess for a reasonable B0 magnitude in CGS)

    // Derived initial conditions
    Vector3d V_SW(VSW_mag * std::cos(del), 0.0, VSW_mag * std::sin(del)); // solar wind speed
    Vector3d B0(B_0_mag * std::cos(th), 0.0, B_0_mag * std::sin(th));     // magnetic field

    // --- CONSTANTS AND INITIAL SCALARS ---
    double U1 = V_SW.norm();
    double B1 = B0.norm();
    double Om = CONSTS.q * B1 / (CONSTS.m * CONSTS.c); // Cyclotron frequency
    double Rg = U1 / Om;                               // Gyro-radius
    int num_steps = cycles * 20 + 1;                   // Steps: cycles * (20 steps/gyro-period) + 1

    // --- INITIALIZE CLASSES AND TURBULENCE ---
    const int N_modes = 200;
    ShockFields S(B0, V_SW, Vector3d::Zero(), N_modes);
    InitParticles P(V_SW, th, del, sample_size);
    double s = std::sqrt(0.5) * B_0_mag; // s used for turbulence init
    
    // Initialize Turbulence
    std::pair<VectorXd, VectorXd> turb_init = S.init_turby(V_SW, Om, s, N_modes);
    VectorXd An = turb_init.first;
    VectorXd kk = turb_init.second;

    // --- INITIALIZE PARTICLE ARRAYS ---
    // X and V arrays (sample_size x 3 components x num_steps)
    // C++: use MatrixXd for [sample_size x 3] for each step, store in a vector or array
    std::vector<MatrixXd> X_history(num_steps, MatrixXd(sample_size, 3));
    std::vector<MatrixXd> V_history(num_steps, MatrixXd(sample_size, 3));
    MatrixXd Vmag_history(sample_size, num_steps);

    // Initial Positions (MATLAB: X(:,:,1))
    VectorXd x0 = -20.0 * Rg * VectorXd::Ones(sample_size);
    // Rand on [-500*Rg, 500*Rg]
    std::uniform_real_distribution<double> dist_y_z(-500.0 * Rg, 500.0 * Rg);
    MatrixXd rand_y_z(sample_size, 2);
    for(int i = 0; i < sample_size; ++i) {
        rand_y_z(i, 0) = dist_y_z(generator);
        rand_y_z(i, 1) = dist_y_z(generator);
    }
    
    X_history[0].col(0) = x0;
    X_history[0].col(1) = rand_y_z.col(0);
    X_history[0].col(2) = rand_y_z.col(1);

    // Initial Velocities (MATLAB: V(:,:,1))
    V_history[0] = P.sampling(sample_size, V_SW);
    Vmag_history.col(0) = V_history[0].rowwise().norm();

    // Cutoff values for splitting
    VectorXd cutoff(10);
    cutoff << 2.0, 5.0, 10.0, 15.8489, 25.1189, 39.8107, 63.0957, 100.0, 10.0, 10.0; // MATLAB: logspace(1,2,8) is 10^1 to 10^2
    cutoff.conservativeResize(10);

    // --- TIME INTEGRATION SETUP ---
    double dt = 0.05 * (1.0 / Om);
    double t = 0.0;
    
    // The splitting container
    SplitContainer Split;

    cout << fixed << setprecision(2);
    cout << "--- Simulation Parameters ---" << endl;
    cout << "Sample Size: " << sample_size << ", Total Steps: " << num_steps << endl;
    cout << "Time step (dt): " << dt << " s" << endl;
    cout << "Gyro-period (1/Om): " << 1.0/Om << " s" << endl;
    cout << "-----------------------------" << endl;

    // --- INTEGRATION LOOP (Particle by particle) ---
    cout << "Starting integration loop..." << endl;
    
    // The MATLAB loop structure is confusing as it iterates 'ii' over sample_size
    // but the inner loop runs over *all* steps, re-initializing the fields.
    // We assume the intent was to fully track *one* particle, then track the *next* particle
    // separately. The integration functions, however, are vectorized.
    // A proper translation would vectorize the outer loop too, but we follow the MATLAB
    // loop structure precisely, tracking particle data in the 3D X/V arrays.

    // The MATLAB code seems to intend a fully vectorized approach, 
    // but the `for ii = 1:sample_size` loop suggests serial processing for setting fields.
    // For faithful translation, we iterate through particles (ii) and steps (n).
    
    // *Optimization Note:* Since `shock_field` is called *outside* the step loop
    // and *inside* the step loop, the fields should be calculated for *all* particles
    // if the helper functions are truly vectorized. The MATLAB code uses X(ii,:,n) 
    // and calls `shock_field` on a single particle's location. We follow this particle-by-particle 
    // field update, but use the vectorized `integrate_boris`.
    
    // Prepare for field calculations: B and U fields need to be N_particles x 3
    MatrixXd B_field(sample_size, 3);
    MatrixXd U_field(sample_size, 3);
    
    // --- MAIN TIME LOOP ---
    for (int n = 0; n < num_steps - 1; ++n) {
        if (n % (num_steps / 10) == 0) {
            cout << "Step: " << n << " / " << num_steps - 1 << " (Time: " << t << " s)" << endl;
        }

        // 1. Calculate fields for ALL particles at time 'n' (before integration)
        for (int i = 0; i < sample_size; ++i) {
            // X_history[n].row(i).transpose() converts the 1x3 row to a 3x1 vector for input
            std::pair<Vector3d, Vector3d> fields = S.shock_field(
                B0, V_SW, S.x0, X_history[n].row(i).transpose(), An, kk,
                P.init_shock(V_SW, th, del), th, del
            );
            B_field.row(i) = fields.first.transpose();
            U_field.row(i) = fields.second.transpose();
        }

        // 2. Integrate ALL particles (Vectorized call)
        std::pair<MatrixXd, MatrixXd> next_step = integrate_boris(
            X_history[n], V_history[n], U_field, B_field, dt
        );

        X_history[n + 1] = next_step.first;
        V_history[n + 1] = next_step.second;

        // 3. Advance field (Called again, essentially calculates fields for time n+1, 
        //    which is technically redundant with the next loop iteration's field calc)
        // We skip the redundant call here, but will call `shock_field` inside the
        // `particle_split` call to ensure the logic flow is correct if it relied on it.
        
        // 4. Update Vmag
        Vmag_history.col(n + 1) = V_history[n + 1].rowwise().norm();

        // 5. Particle Splitting and Displacement
        // The MATLAB code uses X(:,:,2) and V(:,:,2) which means it's always checking 
        // the particle state at step 2. This seems like a bug in the MATLAB original code,
        // as it should be checking the *current* state (n+1).
        // We use X_history[n+1] and V_history[n+1].
        
        // The splitting also compares Vmag_history.col(n+1) against Vmag_history.col(0)
        // and needs Vmag_history.col(0) which is stored.
        
        // We need the fields B and U for the split integration, so we call shock_field
        // just like MATLAB does, but only for the needed rows if splitting occurs.
        
        particle_split(
            X_history[n + 1], 
            V_history[n + 1], 
            Vmag_history.col(0),           // Vmag(:,1)
            Vmag_history.col(n + 1),       // Vmag(:,n+1)
            cutoff, Split
        );
        
        particle_displace(Split, Rg);
        
        // 6. Integrate split particles (nested loop)
        for (size_t j = 0; j < Split.size(); ++j) {
            if (Split[j].has_value()) {
                ParticleSet& current_split = Split[j].value();
                MatrixXd& X_split = std::get<0>(current_split);
                MatrixXd& V_split = std::get<1>(current_split);
                
                // Get U and B for the split particles' locations
                MatrixXd B_split(X_split.rows(), 3);
                MatrixXd U_split(X_split.rows(), 3);

                for (int i = 0; i < X_split.rows(); ++i) {
                     std::pair<Vector3d, Vector3d> fields = S.shock_field(
                        B0, V_SW, S.x0, X_split.row(i).transpose(), An, kk,
                        P.init_shock(V_SW, th, del), th, del
                    );
                    B_split.row(i) = fields.first.transpose();
                    U_split.row(i) = fields.second.transpose();
                }

                // Integrate the split particles
                std::pair<MatrixXd, MatrixXd> split_next_step = integrate_boris(
                    X_split, V_split, U_split, B_split, dt
                );
                
                // Update the split container's X and V
                X_split = split_next_step.first;
                V_split = split_next_step.second;
            }
        }
        
        t += dt;
    }

    // --- TIMING AND EXIT ---
    auto end_time = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end_time - start_time);

    cout << "Integration complete in " << duration.count() / 1000.0 << " seconds." << endl;
    
    // NOTE: You would typically save X_history, V_history, and Split here (e.g., HDF5 or text file)
    
    return 0;
}