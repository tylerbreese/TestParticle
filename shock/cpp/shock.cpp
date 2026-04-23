#include <iostream>
#include <iomanip>
#include <chrono>
#include <random> 
#include <cmath>
#include <string>
#include <cstdio>
#include <fstream>
#include <armadillo>
#include "inputs.hpp"
using namespace std;
using namespace arma;

// --- Function Definitions ---
vec pdfrnd(const vec& x, const vec& px, int sampleSize) {
    vec cdf = cumsum(px);
    cdf /= as_scalar(sum(px));
    vec rnd = randu<vec>(sampleSize);
    vec X;
    interp1(cdf, x, rnd, X);
    X.replace(datum::nan, 0.0);
    return X;
}
mat sampling(const int& sample_size, const mat& U0) {
    double Umin = 1000.0; // cm/s
    double Umax = as_scalar(vecnorm(U0,2,1)); 
    mat V_shock(sample_size, 3);

    vec v = linspace<vec>(Umin, Umax, 100);
    vec u   = pdfrnd(linspace<vec>(-1, 1, 2001), ones<vec>(2001), sample_size);
    vec phi = pdfrnd(linspace<vec>(0, 6.28, 2001), ones<vec>(2001), sample_size);

    double area = as_scalar(trapz(v, v));
    vec F = v / area;

    vec V_plasma = pdfrnd(v, F, sample_size);
    vec vx = V_plasma % sqrt(1.0 - pow(u, 2)) % cos(phi);
    vec vy = V_plasma % sqrt(1.0 - pow(u, 2)) % sin(phi);
    vec vz = V_plasma % u;

    V_shock.col(0) = vx + U0(0); 
    V_shock.col(1) = vy + U0(1);
    V_shock.col(2) = vz + U0(2);
    
    return V_shock; 
}
struct Boris {
    mat Xnew, Vnew;
};
Boris integrate_boris(mat& Xn, mat& Vn, mat& U, mat& B, double dt) {
    const double c = inputs::c;
    const double m = inputs::m;
    const double q = inputs::q;

    mat E(1,3); // electric field
    mat X(1,3);
    mat V(1,3);
    mat Vp(1,3);
    mat Vq(1,3);

    E.col(0) = (-1/c) * (U.col(1)*B.col(2) - U.col(2)*B.col(1));
    E.col(1) = (-1/c) * (U.col(2)*B.col(0) - U.col(0)*B.col(2));
    E.col(2) = (-1/c) * (U.col(0)*B.col(1) - U.col(1)*B.col(0)); 

    mat u = Vn;
    mat R = (q/m) * E * (dt/2);
    mat T = (q/(m*c)) * B * (dt/2);
    mat S = (2.0 * T) / (1.0 + pow(norm(T),2));

    u = u + R; //first step
    //twirling
    Vp.col(0) = u.col(0) + ( u.col(1)*T.col(2) - u.col(2)*T.col(1) );
    Vp.col(1) = u.col(1) + ( u.col(2)*T.col(0) - u.col(0)*T.col(2) );
    Vp.col(2) = u.col(2) + ( u.col(0)*T.col(1) - u.col(1)*T.col(0) );
    //twirling
    Vq.col(0) = u.col(0) + ( Vp.col(1)*S.col(2) - Vp.col(2)*S.col(1) );
    Vq.col(1) = u.col(1) + ( Vp.col(2)*S.col(0) - Vp.col(0)*S.col(2) );
    Vq.col(2) = u.col(2) + ( Vp.col(0)*S.col(1) - Vp.col(1)*S.col(0) );

    V = Vq + R; // next half step

    X.col(0) = Xn.col(0) + V.col(0) * dt;
    X.col(1) = Xn.col(1) + V.col(1) * dt;
    X.col(2) = Xn.col(2) + V.col(2) * dt;

    Boris result;
    result.Xnew = X;
    result.Vnew = V;
    return result;

}
struct InitShock {
    double r, a, b;
    mat U2;
};
InitShock init_shock(mat U1, double th, double del, double Va, double Vs) {
    // --- Set Manually ---
    // Va = 30e5; % cm/s
    // Vs = 70e5; % cm/s
    // --- Calc Mach Numbers ---
    double Ms  = norm(U1) / Vs;
    double Ma  = norm(U1) / Va;
    double Ms2 = pow(Ms,2);
    double Ma2 = pow(Ma,2);

    double g = 1.67; // adiabatic index g = 5/3 in SW g = 7/5 in air
    double r = (g + 1) / ( (g - 1) + (2/Ms2) );
    double a = r * (Ma2*cos(del)*cos(del) - cos(th)*cos(th)) / (Ma2*cos(del)*cos(del) - r*cos(th)*cos(th));
    double b = ((r-1)*cos(th)*sin(th)) / (Ma2*cos(del)*cos(del) - r*cos(th)*cos(th));

    mat U2(1,3); // %shocked field 
    U2.col(0) = U1.col(0) / r;
    U2.col(1) = U1.col(1);
    U2.col(2) = U1.col(2) + b * U1.col(0);

    InitShock result;
    result.r  = r;
    result.a  = a;
    result.b  = b;
    result.U2 = U2;
    return result;
}
struct Field {
    mat Unow, Bnow;
};
Field shock_field(mat X, double x0, mat U0, mat B0, double r, double a, double b) {

    mat U(1,3); // flow field
    mat B(1,3); //magnetic field
    double x = as_scalar(X.col(0));
    if ( (x-x0) > 0.0 ) {
        U.col(0) = U0.col(0) / r;
        U.col(1) = U0.col(1);
        U.col(2) = U0.col(2) + b * U0.col(0);
        B.col(0) = B0.col(0);
        B.col(1) = B0.col(1);
        B.col(2) = a * B0.col(2);
    }
    else {
        U = U0;
        B = B0;
    }
    Field result;
    result.Unow = U;
    result.Bnow = B;
    return result;
}

struct TurbyInit {
    mat An, Bn, K;
};
TurbyInit init_turby(double Lmin, double Lmax, int N, double s2, mat U0, mat vA) {
    double kmax = 2.0 * inputs::PI / Lmin;
    double kmin = 2.0 * inputs::PI / Lmax;
    double D = (log10(kmax) - log10(kmin)) / (N - 1);

    vec kk(N);
    kk(0) = kmin;
    for (int i = 1; i < N; ++i) {
        kk(i) = pow(10.0, D) * kk(i - 1);
    }
    double dk = pow(10.0, D) - 1.0;

    // PSD calculation
    double Lc = 0.167 * 1.496e13;
    vec G = (4.0 * inputs::PI * pow(kk, 3) * dk) / (1.0 + pow(kk * Lc, 11.0 / 3.0));
    vec Amp = sqrt(s2 * G / as_scalar(sum(G)));
    // Random angles
    vec alpha = 2.0 * inputs::PI * randu<vec>(N);
    vec beta  = 2.0 * inputs::PI * randu<vec>(N);
    vec phi   = 2.0 * inputs::PI * randu<vec>(N);
    vec theta = inputs::PI * randu<vec>(N);

    // Amplitudes and Wave Vectors
    vec Anx = Amp % cos(alpha) % cos(phi) % cos(theta);
    vec Any = Amp % cos(alpha) % sin(phi) % cos(theta);
    vec Anz = Amp % (-1.0 * cos(alpha)) % sin(phi); // Explicitly a vec
    vec Bnx = Amp % sin(alpha) % sin(phi);
    vec Bny = Amp % (-1.0 * sin(alpha)) % cos(phi);
    vec Bnz = zeros<vec>(N); // Make sure this is size N
    vec knx = kk % cos(phi) % sin(theta);
    vec kny = kk % sin(phi) % sin(theta);
    vec knz = kk % cos(theta);

    mat combined_V = U0 + vA;
    vec ww = knx * combined_V(0) + kny * combined_V(1) + knz * combined_V(2);

    TurbyInit result;
    result.An.set_size(N,3);
    result.An.col(0) = Anx; 
    result.An.col(1) = Any; 
    result.An.col(2) = Anz; 
    result.Bn.set_size(N,3);
    result.Bn.col(0) = Anx; 
    result.Bn.col(1) = Any; 
    result.Bn.col(2) = Anz; 
    result.K.set_size(N, 5); // idk why arma freaks out if not done this way
    result.K.col(0) = knx;
    result.K.col(1) = kny;
    result.K.col(2) = knz;
    result.K.col(3) = ww;
    result.K.col(4) = beta;
    
    return result;
}
mat sum_turby(mat X, double T, mat K, mat An, mat Bn) {

    vec arg = (K.col(0) * X(0) + K.col(1) * X(1) + K.col(2) * X(2)) + K.col(4);
    //vec arg = sign * (K.col(0) * X(0) + K.col(1) * X(1) + K.col(2) * X(2)) + K.col(3) * T + K.col(4);
    
    mat dB(1,3);
    dB.col(0) = sum(An.col(0) % cos(arg) + Bn.col(0) % sin(arg));
    dB.col(1) = sum(An.col(1) % cos(arg) + Bn.col(1) % sin(arg));
    dB.col(2) = sum(An.col(2) % cos(arg) + Bn.col(2) % sin(arg));
    
    return dB;
}

double boundary(mat U1, mat V, double r) {

    double U2 = as_scalar(norm(U1)) / r;
    double v  = as_scalar(norm(V));
    double arg = (v-U2) / (v+U2);
    double P = 1 - pow(arg,2);

    return P;
}

struct Displace {
    mat Xdis, Vdis;
};
Displace particle_displace(mat X, mat V, double Rg){
    vec u   = pdfrnd(linspace<vec>(-1.0, 1.0, 2001), ones<vec>(2001), 1);
    vec phi = pdfrnd(linspace<vec>(0.0, 6.28, 2001), ones<vec>(2001), 1);
    mat Xdis(1,3);
    mat Vdis(1,3);

    Xdis.col(0) = X.col(0) + randu() * Rg;
    Xdis.col(1) = X.col(1) + randu() * Rg;
    Xdis.col(2) = X.col(2) + randu() * Rg;
    Vdis.col(0) = V.col(0) * sqrt(1.0 - pow(u, 2)) % cos(phi);
    Vdis.col(1) = V.col(1) * sqrt(1.0 - pow(u, 2)) % sin(phi);
    Vdis.col(2) = V.col(2) * u;

    Displace result;
    result.Xdis = Xdis;
    result.Vdis = Vdis;
    return result;
};

string get_timestamp_filename(string prefix, string extension) {
    // 1. Get the current time from the system clock
    auto now = chrono::system_clock::now();
    auto in_time_t = chrono::system_clock::to_time_t(now);

    // 2. Format the time into a string buffer
    stringstream ss;
    // Format: YearMonthDay_HourMinuteSecond
    ss << prefix << put_time(localtime(&in_time_t), "%Y-%m-%d_%H-%M-%S") << extension;
    
    return ss.str();
}

// --- Main Program ---
int main() {
    arma_rng::set_seed_random(); // random numbers
    // // --- constants ---
    const double c = 3E10; //speed of light, m/s
    const double e = 4.8E-10; // elementary charge, cgs
    const double Z = inputs::Z; // atomic number
    const double A = inputs::A; // mass number
    const double q = inputs::q; // ion charge 
    const double m = inputs::m; // grams

    // // --- Inputs ---
    const int sample_size = inputs::sample_size;
    const double th  = inputs::th;
    const double del =  inputs::del;
    mat U0  = inputs::init_U0();
    mat B0  = inputs::init_B0();
    mat Usw = inputs::init_Usw();

    double Om = inputs::get_Om();
    double Rg  = inputs::get_Rg();
    const double B1 = inputs::get_B1(); // field strength
    const double s2 = inputs::var * pow(B1,2); // wave spectrum

    const double cycles = inputs::cycles; //number of time segments
    const double dt = 0.05 * (1/Om);
    const int num_steps =  cycles/(0.05)  + 1; //number of steps
    vec T = regspace(0, dt, cycles/Om);

    // // --- Initialize ---
    mat V0(sample_size,3);
    V0 = sampling(sample_size,Usw);
    //vec xinit = (2.0 * randu<vec>(sample_size) - 1.0) * (100.0 * Rg);
    vec xinit = ones<vec>(sample_size) * (100.0 * Rg); // drop in particles upstream in fixed plane'
    // vec yinit = zeros<vec>(sample_size);
    // vec zinit = zeros<vec>(sample_size);
    vec yinit = (2.0 * randu<vec>(sample_size) - 1.0) * (100.0 * Rg);
    vec zinit = (2.0 * randu<vec>(sample_size) - 1.0) * (100.0 * Rg);

    cube X = zeros<cube>(sample_size, 3, num_steps);
    cube V = zeros<cube>(sample_size,3,num_steps);
    
    double x0 = 0.0 * Rg; // shock location IMPORTANT
    X.slice(0).col(0) = x0 - xinit;
    X.slice(0).col(1) = yinit;
    X.slice(0).col(2) = zinit;
    V.slice(0).col(0) = V0.col(0);
    V.slice(0).col(1) = V0.col(1);
    V.slice(0).col(2) = V0.col(2);

    // --- Generate output files ---
    string simfile   = get_timestamp_filename("sim_data_", ".csv");
    string splitinit = get_timestamp_filename("split_init_", ".csv");
    string splitout  = get_timestamp_filename("split_data_", ".csv");

    ofstream f1(simfile);
    ofstream f2(splitinit);
    ofstream f3(splitout);
    
    f1 << "PID,x0,xf,E0,Ef" << endl;
    f2 << "PID,tstep,weight,x,y,z,vx,vy,vz" << endl;
    f3 << "PID,weight,x0,xf,E0,Ef" << endl;
    
    // --- Init Shock and Turby
    double vA =  30.0e5; // cm/s
    double Cs = 149.0e5; // cm/s
    mat U1 = U0 - Usw; // IP shock use this
    //mat U1 = Usw; // if doing TS use this
    mat V_A(1,3); //vector alfven 
    V_A.col(0) = vA * cos(th); // follows from B 
    V_A.col(1) = 0.0;
    V_A.col(2) = vA * sin(th);
    //double r, a, b; mat U2;  
    auto [r,a,b,U2] = init_shock(U1,th,del,vA,Cs);
    const double N    = 201;
    const double Lmin = 0.5*Rg;
    const double Lmax = 1.496e13; 
    auto [An,Bn,K]  = init_turby(Lmin,Lmax,N,s2,U1,V_A);
    mat Btrack(T.n_elem,4);
    // --- Integration ---
    for (int i = 0; i < sample_size; i++) {

        const vec cutoff = { 2.0, 3.0, 5.0, 10.0, 15.0, 20.0, 25.0, 50.0, 75.0, 100.0 }; // splitting thresholds
        const double En0 = 0.5 * m * (V(i,0,0)*V(i,0,0) + V(i,1,0)*V(i,1,0) + V(i,2,0)*V(i,2,0));
        uword thresholds_passed = 0; // track thresholds for splitting

        for (int j = 0; j < num_steps; j++){

            mat dB = sum_turby(X.slice(j).row(i),T(j),K,An,Bn);
            auto [Unow,Bnow] = shock_field(X.slice(j).row(i), x0, U1, B0+dB, r, a, b);
            Btrack(j,0) = as_scalar(Bnow.col(0));
            Btrack(j,1) = as_scalar(Bnow.col(1));
            Btrack(j,2) = as_scalar(Bnow.col(2));
            Btrack(j,3) = as_scalar(norm(Bnow));
            mat Xold = X.slice(j).row(i);
            mat Vold = V.slice(j).row(i);

            double x = as_scalar(X.slice(j).row(i).col(0));
            if ( (x-x0) > 250.0*Rg ) {
                double P_return = randu();
                double P_escape = boundary(U1,V.slice(j).row(i),r);
                if (P_return < P_escape){ 
                    cout << "Particle " << i << " has left the building \n";
                    cout << "Particle " << i << " lasted " << j << " steps \n";
                    rowvec last_X = X.slice(j).row(i);
                    rowvec last_V = V.slice(j).row(i);
                    // Fill all remaining slices from ii+1 to the end
                    for (int s = j + 1; s < X.n_slices; s++) {
                        X.slice(s).row(i) = last_X;
                        V.slice(s).row(i) = last_V;
                    }
                    break;
                }
                auto [Xdis, Vdis] = particle_displace(Xold, Vold, Rg);
                Vold = Vdis; // "scatter" at the boundary
            }

            auto [Xnew, Vnew] = integrate_boris(Xold, Vold, Unow, Bnow, dt);
            X.slice(j+1).row(i) = Xnew;
            V.slice(j+1).row(i) = Vnew;

            // --- Splitting --- 
            double Enf = 0.5 * m * (V(i,0,j)*V(i,0,j) + V(i,1,j)*V(i,1,j) + V(i,2,j)*V(i,2,j));
            double ratio = Enf / En0;
            // Determine how many total thresholds the particle is currently qualified for
            uword current_count = accu(ratio > cutoff);
            if (current_count > thresholds_passed) {
                uword new_levels = current_count - thresholds_passed;
                for (uword lev = 0; lev < new_levels; ++lev) {
                    for (uword ii = 1; ii <= current_count; ++ii) {
                        auto [Xspl, Vspl] = particle_displace(Xnew, Vnew, Rg);
                        f2 << i << ","       // Particle ID
                        << j << ","       // Time Step
                        << ii << ","      // Split Index
                        << Xspl(0) << "," << Xspl(1) << "," << Xspl(2) << "," 
                        << Vnew(0) << "," << Vnew(1) << "," << Vnew(2) << "\n";
                    }
                }
                thresholds_passed = current_count;
    }

        }
    }
    cout << "all particles tested" << endl;
    for (int k = 0; k < sample_size; k++) {
        f1 << k << "," 
                << X(k,0,0) << ","
                << X(k,0,num_steps-1) << ","
                << 0.5 * m * ( pow(V(k,0,0),2) + pow(V(k,1,0),2) + pow(V(k,2,0),2) ) << ","
                << 0.5 * m * ( pow(V(k,0,num_steps-1),2) + pow(V(k,1,num_steps-1),2) + pow(V(k,2,num_steps-1),2) ) << "\n";
    }
    cout << "Data saved to sim_data.csv" << endl;
    
    // --- Splitting Integration ---
    cout << "Begin splitting" << endl;
    mat data;
    bool success = data.load(splitinit, csv_ascii);
    int num_rows = data.n_rows;
    if (success) {
        cout << "Data loaded successfully!" << endl;
        for (int ii = 0; ii < num_rows; ii++) {

            mat Xspl = zeros<mat>(num_steps,3);
            mat Vspl = zeros<mat>(num_steps,3);
            Xspl.row(0) = data(ii,span(3,5));
            Vspl.row(0) = data(ii,span(6,8));
            int tstart = as_scalar(data(ii,1));
            const double En0 = 0.5 * m * (Vspl(tstart,0)*Vspl(tstart,0) + Vspl(tstart,1)*Vspl(tstart,1) + Vspl(tstart,2)*Vspl(tstart,2));
            for (int jj = tstart; jj < num_steps; jj++){

                mat Xold = Xspl.row(jj);
                mat Vold = Vspl.row(jj);         

                mat dB = sum_turby(Xold,T(jj),K,An,Bn);
                auto [Unow,Bnow] = shock_field(Xold, x0, U1, B0+dB, r, a, b);
                double x = as_scalar(Xold.col(0));
                if ( (x-x0) > 250.0*Rg ) {
                    double P_return = randu();
                    double P_escape = boundary(U1,Vold,r);
                    if (P_return < P_escape){ 
                        cout << "Particle " << ii << " has left the building \n";
                        cout << "Particle " << ii << " lasted " << jj-data(ii,1) << " steps \n";
                        rowvec last_X = Xspl.row(jj);
                        rowvec last_V = Vspl.row(jj);
                        // Fill all remaining slices from ii+1 to the end
                        for (int s = jj + 1; s < X.n_slices; s++) {
                            Xspl.row(s) = last_X;
                            Vspl.row(s) = last_V;
                        }
                    break;
                    }
                    auto [Xdis, Vdis] = particle_displace(Xold, Vold, Rg);
                    Vold = Vdis; // "scatter" at the boundary
                }
                auto [Xnew, Vnew] = integrate_boris(Xold, Vold, Unow, Bnow, dt);
                Xspl.row(jj+1).col(0) = Xnew(0,0);
                Xspl.row(jj+1).col(1) = Xnew(0,1);
                Xspl.row(jj+1).col(2) = Xnew(0,2);
                Vspl.row(jj+1).col(0) = Vnew(0,0);
                Vspl.row(jj+1).col(1) = Vnew(0,1);
                Vspl.row(jj+1).col(2) = Vnew(0,2);
            }
            //cout << "Shape: " << Vspl.n_rows << "x" << Vspl.n_cols << std::endl;
            //cout << "Numsteps " << num_steps << endl;
            double Enf = 0.5 * m * (Vspl(num_steps-1,0)*Vspl(num_steps-1,0) + Vspl(num_steps-1,1)*Vspl(num_steps-1,1) + Vspl(num_steps-1,2)*Vspl(num_steps-1,2));
            // f3 << "PID,weight,x0,xf,E0,Ef" << endl;
            f3 << data(ii,0) << "," << data(ii,2) << ","
                             << data(ii,3) << "," << Xspl(num_steps-1,0) << ","
                             << En0 << "," << Enf << "\n";
        }
    }
    cout << "Data saved to split_data.csv" << endl;
    f1.close();
    f2.close();
    f3.close();
    return 0;
}

// run cmd
//g++ -std=c++17 prototype.cpp ip_snowplow.cpp -o prototype.exe -larmadillo 


// old output format
// good for debugging
    // --- Output ---
    // ofstream ef("distribution_data.csv");
    // ef << "PID,x0,xf,E0,Ef" << endl;
    // for (int k = 0; k < sample_size; k++) {
    //     ef << k << "," 
    //             << X(k,0,0) << ","
    //             << X(k,0,num_steps-1) << ","
    //             << 0.5 * m * ( pow(V(k,0,0),2) + pow(V(k,1,0),2) + pow(V(k,2,0),2) ) << ","
    //             << 0.5 * m * ( pow(V(k,0,num_steps-1),2) + pow(V(k,1,num_steps-1),2) + pow(V(k,2,num_steps-1),2) ) << "\n";
    // }
    // int pick = randi(distr_param(0, sample_size-1));
    // ofstream xf("position_data.csv");
    // xf << "Step,X,Y,Z" << endl;
    // for (int i = 0; i < num_steps; i++) {
    //     xf << i << "," 
    //             << X(pick, 0, i) << "," 
    //             << X(pick, 1, i) << "," 
    //             << X(pick, 2, i) << "\n"; 
    // }
    // xf.close();
    // ofstream vf("velocity_data.csv");
    // vf << "Step,Vx,Vy,Vz" << endl;
    // for (int i = 0; i < num_steps; i++) {
    //     vf << i << "," 
    //             << V(pick, 0, i) << "," 
    //             << V(pick, 1, i) << "," 
    //             << V(pick, 2, i) << "\n"; 
    // }
    // vf.close();
    // ofstream bf("magfield_data.csv");
    // bf << "Step,Bx,By,Bz,Bm" << endl;
    // for (int i = 0; i < num_steps; i++) {
    //     bf << i << "," 
    //             << Btrack(i, 0) << "," 
    //             << Btrack(i, 1) << ","
    //             << Btrack(i, 2) << "," 
    //             << Btrack(i, 3) << "\n"; 
    // }
    // bf.close();    // --- Output ---
    // ofstream ef("distribution_data.csv");
    // ef << "PID,x0,xf,E0,Ef" << endl;
    // for (int k = 0; k < sample_size; k++) {
    //     ef << k << "," 
    //             << X(k,0,0) << ","
    //             << X(k,0,num_steps-1) << ","
    //             << 0.5 * m * ( pow(V(k,0,0),2) + pow(V(k,1,0),2) + pow(V(k,2,0),2) ) << ","
    //             << 0.5 * m * ( pow(V(k,0,num_steps-1),2) + pow(V(k,1,num_steps-1),2) + pow(V(k,2,num_steps-1),2) ) << "\n";
    // }
    // int pick = randi(distr_param(0, sample_size-1));
    // ofstream xf("position_data.csv");
    // xf << "Step,X,Y,Z" << endl;
    // for (int i = 0; i < num_steps; i++) {
    //     xf << i << "," 
    //             << X(pick, 0, i) << "," 
    //             << X(pick, 1, i) << "," 
    //             << X(pick, 2, i) << "\n"; 
    // }
    // xf.close();
    // ofstream vf("velocity_data.csv");
    // vf << "Step,Vx,Vy,Vz" << endl;
    // for (int i = 0; i < num_steps; i++) {
    //     vf << i << "," 
    //             << V(pick, 0, i) << "," 
    //             << V(pick, 1, i) << "," 
    //             << V(pick, 2, i) << "\n"; 
    // }
    // vf.close();
    // ofstream bf("magfield_data.csv");
    // bf << "Step,Bx,By,Bz,Bm" << endl;
    // for (int i = 0; i < num_steps; i++) {
    //     bf << i << "," 
    //             << Btrack(i, 0) << "," 
    //             << Btrack(i, 1) << ","
    //             << Btrack(i, 2) << "," 
    //             << Btrack(i, 3) << "\n"; 
    // }
    // bf.close();
    // // Print the first 10 results for shock speed and Mach numbers
    // cout << "Step | Shock Speed (U) | Sound Speed (Cs) | Alfven Speed (vA)" << endl;
    // for (uword i = 0; i < 10 && i < U.n_elem; ++i) {
    //     printf("%4d | %15.2f | %15.2f | %15.2f\n", 
    //             (int)i, U(i), Cs(i), vA(i));
    // }