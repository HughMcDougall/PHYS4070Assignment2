//
// Main runtime for PHYS4070 assignment 2 part 1.
//

#include <iostream>
#include <vector>
#include <cmath>
#include <fstream>
#include <cassert>

#include "vector_utils.hpp"
#include "params.hpp"
#include "forces_and_integrators.hpp"
#include "sysvec_utils.hpp"
#include "manybody_example.hpp"

using vec = std::vector<double>;
using namespace params;

using sysvec  = std::vector<std::vector<double>>;
using diff_func  = std::function<sysvec(sysvec)>;

int main(int argc, char *argv[]) {
    /// Main simulation runtime for PHYS4070 Assignment 2 part 1
    /// Inputs are:
    /// mode        int     1-4         1                           Determines the type of simulation to run
    /// output_dir  str                 results/sim_results.dat     Location to save outputs to
    /// tmax        double  >0          100                         Max sim time
    /// v_launch    double  >0          1.1376                      Launch velocity for modes 1-3
    /// thet        double  >0          -1.08                       Starting angle of moon in orbit, radians
    /// dt          double  >0, <tmax   0.01                        Timestep for RK4 integration
    /// sparseness  int     >0, <maxits 1                           Number of int steps per output line

    //------------------------------------------------------------------
    // SIMULATION PARAMETERS (Defaults)

    // Mode selection
    int mode = 3;
    bool use_nbody = true;
    bool do_kick = false;

    // Sim Times
    double tmax = 100;
    double dt = 0.01;

    //Launch tuning params
    double thet = 2 * M_PI * (-62.1) / 360;
    thet = -89.7 / pow(r_moon, 1.5);
    double v_launch = params::v_proj_launch;

    //Output params
    int sparseness = 1;
    std::string output_dir = "results/sim_results.dat";

    //------------------------------------------------------------------
    // INPUTS
    // order of inps: mode, outpir_dir, tmax, v_launch, thet, dt, sparseness
    {
        if (argc > 1) { mode = std::stoi(argv[1]); }
        if (argc > 2) { output_dir = argv[2]; }
        if (argc > 3) { tmax = std::stod(argv[3]); }
        if (argc > 4) { v_launch = std::stod(argv[4]); }
        if (argc > 5) { thet = std::stod(argv[5]); }
        if (argc > 6) { dt = std::stod(argv[6]); }
        if (argc > 7) { sparseness = std::stoi(argv[7]); }
    }

    //------------------------------------------------------------------
    // INTERPRET INPUTS

    assert(mode<=4 && mode>0                            && "Incorrect mode input. Mode must be 0-4");
    assert(v_launch>0                                   && "Cannot have negative launch velocity");
    assert(dt > 0 && tmax > dt                          && "Bad time inputs. Need tmax > dt > 0");
    assert(sparseness>0 && sparseness<floor(tmax/dt) && "Sparseness must be between 1 and max its");

    if(mode==1){
        use_nbody = false;
        do_kick   = false;
    } else if(mode==2){
        use_nbody = true;
        do_kick   = false;
    } else if (mode==3){
        use_nbody = true;
        do_kick   = true;
    } else if(mode==4) {
        use_nbody = true;
        do_kick   = false;
    }

    //------------------------------------------------------------------
    // INITIAL CONDITIONS

    // Projectile
    vec proj_pos = {R_planet, 0};
    vec proj_vel = {v_launch, 0};

    // Moon
    vec moon_pos = {r_moon*cos(thet), r_moon*sin(thet)};
    vec moon_vel = {v_moon*-sin(thet), v_moon*cos(thet)};

    // Construct system vector
    sysvec X = {proj_pos, moon_pos, proj_vel, moon_vel};

    // Force function dX/dt
    std::vector<double> M = {m_proj, m_moon};
    std::vector<bool> LIVE = {0, 1};

    if (mode==4){ M=manybody::M; LIVE = manybody::LIVE; X = manybody::X;}

    // Determine which dX/dt = f(X) to use
    std::function<sysvec(sysvec)> f;
    if (use_nbody){
        f = [M,LIVE](sysvec X) { return f_nbody(X,M,LIVE);};
    } else{
        f = [](sysvec X) { return f_planet(X);};
    }

    //------------------------------------------------------------------
    // INTEGRATION LOOP
    std::cout<<"Doing simulation with "<<floor(tmax / dt ) << " itterations, outputting every "<< sparseness<< "\n";
    std::cout<<"Starting integration...\n";

    double t = 0;
    std::ofstream output(output_dir);
    vec dr;
    vec kick;
    for(int i = 0; i<tmax/dt; i++){
        t+=dt;
        X = runge_step(X,f,dt);

        // Determine distance to moon (Only used if do_kick is active)
        dr = X[1]-X[0];
        if (do_kick and vnorm(dr)<2*R_moon){
            kick = orbital_kick(dr, X[2], X[3]);
            X[2]+= kick;
            do_kick = false;
            std::cout << "kick performed at t= "<<t<<"\n";
            std::cout << "delta-v of kick:\t"<<vnorm(kick)<<"\n";
            std::cout << "Angle of kick:\t"<<atan(kick[0]/kick[1])<<"\n";
        }

        if(i%sparseness==0){
            output<<t<<"\t";
            print_to_stream(X, output);
        }

    }
    output.close();
    std::cout<<"Done\n";

    //------------------------------------------------------------------
    // FINISHING UP
    return 0;
}