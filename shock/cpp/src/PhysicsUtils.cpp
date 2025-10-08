#include "PhysicsUtils.h"
#include <random>
#include <cmath>

// --- Global Constants Definition ---
const Constants CONSTS;

// Global random number generator (for displacement)
static std::default_random_engine generator;
static std::uniform_real_distribution<double> dist_0_1(0.0, 1.0);

// --- integrate_boris method ---
std::pair<Eigen::MatrixXd, Eigen::MatrixXd> integrate_boris(
    const Eigen::MatrixXd& Xn, 
    const Eigen::MatrixXd& Vn, 
    const Eigen::MatrixXd& U,
    const Eigen::MatrixXd& B,
    double dt
) {
    const int sample_size = Xn.rows();
    Eigen::MatrixXd E(sample_size, 3);
    
    // E = (-1/c) * (U x B) -> Convective Electric Field
    // E(:,1) = (-1/c) .* (U(:,2).*B(:,3) - U(:,3).*B(:,2));
    E.col(0) = (-1.0/CONSTS.c) * (U.col(1).array() * B.col(2).array() - U.col(2).array() * B.col(1).array());
    E.col(1) = (-1.0/CONSTS.c) * (U.col(2).array() * B.col(0).array() - U.col(0).array() * B.col(2).array());
    E.col(2) = (-1.0/CONSTS.c) * (U.col(0).array() * B.col(1).array() - U.col(1).array() * B.col(0).array());
    
    // E = 0.0; // Uncomment to turn off E field (as per MATLAB comment)

    // BORIS algorithm implementation (uses element-wise operations)
    Eigen::MatrixXd u = Vn;
    Eigen::MatrixXd R = (CONSTS.q / CONSTS.m) * E * (dt / 2.0);
    Eigen::MatrixXd T = (CONSTS.q / (CONSTS.m * CONSTS.c)) * B * (dt / 2.0);
    
    // S = 2.*T ./ (1 + sum(T.*T,2)); (Sum across the row/vector for each particle)
    Eigen::VectorXd T_dot_T = T.rowwise().squaredNorm();
    Eigen::MatrixXd S = 2.0 * (T.array().colwise() / (1.0 + T_dot_T.array()));
    //Eigen::VectorXd S = 2.0 * (T.array().colwise() / (1.0 + T_dot_T.array()));
    u += R; // first half step
    
    // V' = u + (u x T)
    Eigen::MatrixXd Vp(sample_size, 3);
    Vp.col(0) = u.col(0).array() + ( u.col(1).array() * T.col(2).array() - u.col(2).array() * T.col(1).array() );
    Vp.col(1) = u.col(1).array() + ( u.col(2).array() * T.col(0).array() - u.col(0).array() * T.col(2).array() );
    Vp.col(2) = u.col(2).array() + ( u.col(0).array() * T.col(1).array() - u.col(1).array() * T.col(0).array() );

    // V'' = u + (V' x S)
    Eigen::MatrixXd Vq(sample_size, 3);
    Vq.col(0) = u.col(0).array() + ( Vp.col(1).array() * S.col(2).array() - Vp.col(2).array() * S.col(1).array() );
    Vq.col(1) = u.col(1).array() + ( Vp.col(2).array() * S.col(0).array() - Vp.col(0).array() * S.col(2).array() );
    Vq.col(2) = u.col(2).array() + ( Vp.col(0).array() * S.col(1).array() - Vp.col(1).array() * S.col(0).array() );

    Eigen::MatrixXd V = Vq + R; // next half step
    
    // Position advance
    Eigen::MatrixXd X = Xn + V * dt;
    
    return {X, V};
}

// --- particle_split method ---
SplitContainer particle_split(
    const Eigen::MatrixXd& X_all, 
    const Eigen::MatrixXd& V_all, 
    const Eigen::VectorXd& Vmag_all_start,
    const Eigen::VectorXd& Vmag_all_current,
    const Eigen::VectorXd& cutoff, 
    SplitContainer& split_in_out 
) {
    const int sample_size = X_all.rows();
    const int n_cutoff = cutoff.size();

    // 1. Calculate ratio and flag
    Eigen::VectorXd ratio = Vmag_all_current.array().square() / Vmag_all_start.array().square(); // En / E0
    Eigen::Matrix<bool, Eigen::Dynamic, Eigen::Dynamic> temp_flag_bool(sample_size, n_cutoff);
    for (int j = 0; j < n_cutoff; ++j) {
        temp_flag_bool.col(j) = (ratio.array() > cutoff(j));
    }

    // Initialize/Resize split_in_out if needed
    if (split_in_out.empty()) {
        split_in_out.resize(n_cutoff);
    }
    
    // Handle split accumulation (MATLAB logic is tricky: flag = flag | temp_flag)
    Eigen::Matrix<int, Eigen::Dynamic, Eigen::Dynamic> final_flag_int(sample_size, n_cutoff);
    for(int j = 0; j < n_cutoff; ++j) {
        if (split_in_out[j].has_value()) {
            // Get previous flag vector
            const Eigen::VectorXi& prev_flag = std::get<2>(split_in_out[j].value());
            // Logical OR: new flag is set if it was set before OR if it is set now
            final_flag_int.col(j) = prev_flag.array() | temp_flag_bool.col(j).cast<int>().array();
        } else {
            // No previous split, use current flag
            final_flag_int.col(j) = temp_flag_bool.col(j).cast<int>();
        }
    }
    
    // 2. Populate the SplitContainer
    for (int j = 0; j < n_cutoff; ++j) {
        // Find indices where flag is true
        std::vector<int> indices;
        for (int i = 0; i < sample_size; ++i) {
            if (final_flag_int(i, j) == 1) {
                indices.push_back(i);
            }
        }
        
        if (!indices.empty()) {
            int count = indices.size();
            Eigen::MatrixXd X_split(count, 3);
            Eigen::MatrixXd V_split(count, 3);
            Eigen::VectorXi flag_vector = final_flag_int.col(j);

            // Extract the rows
            for (int k = 0; k < count; ++k) {
                X_split.row(k) = X_all.row(indices[k]);
                V_split.row(k) = V_all.row(indices[k]);
            }
            // Update the split container element
            split_in_out[j] = ParticleSet{X_split, V_split, flag_vector};
        } else {
            split_in_out[j].reset(); // Clear if no particles are split
        }
    }
    
    return split_in_out;
}

// --- particle_displace method ---
void particle_displace(
    SplitContainer& split_in_out, 
    double R_g
) {
    if (split_in_out.empty()) return;

    for (auto& optional_set : split_in_out) {
        if (optional_set.has_value()) {
            ParticleSet& current_set = optional_set.value();
            Eigen::MatrixXd& X_split = std::get<0>(current_set);
            
            const int n_rows = X_split.rows();
            
            // Generate random displacements (rand(N,1) * R_g)
            for (int i = 0; i < n_rows; ++i) {
                X_split(i, 0) += dist_0_1(generator) * R_g;
                X_split(i, 1) += dist_0_1(generator) * R_g;
                X_split(i, 2) += dist_0_1(generator) * R_g;
            }
            // Note: Velocity V_split (std::get<1>) is not modified here.
        }
    }
}