/// singlegrid.cpp

#include <iostream>
#include "grid.hpp"
#include "rand_utils.hpp"
#include "monte_carlo.hpp"



int main(int argc, char *argv[]) {
    /// Runtime to execute and print a single montecarlo run for an NxN grid in Assignment 2 Part 2
    /// Inputs are:
    /// Var         Type    Constraints Default Desc
    /// N           int     >0          8       Dimension of NxN Grid
    /// T           double  >=0         1       Temperature for Monte carlo simulation
    /// Nits        int     >0          1E4     Number of full sweeps of the grid / samples to draw
    /// Nburn       int     >=-0        1E3     Number of burn-in itterations in monte carlo
    /// flips_per   int     >0          0(NxN)  Number of bit flips per 'sample' If set to zero, will use flips_per=NxN
    /// seed        int                 0       Seed for random number generation

    //------------------------------------------------------------------
    // SIMULATION PARAMETERS (Defaults)
    int N=8;
    double T = 2.77*0.75;

    int Nits = (int)1E4;
    int Nburn = (int)1E3;
    int flips_per = 0;
    int seed = 0;

    //------------------------------------------------------------------
    // INPUTS
    {
        if (argc > 1) { N = std::stoi(argv[1]); }
        if (argc > 2) { T = std::stod(argv[2]); }
        if (argc > 3) { Nits = std::stoi(argv[3]); }
        if (argc > 4) { Nburn = std::stoi(argv[4]); }
        if (argc > 5) { flips_per = std::stoi(argv[5]); }
        if (argc > 6) { seed = std::stoi(argv[6]); }
    }

    //------------------------------------------------------------------
    // INTERPRET INPUTS
    assert(N>0              && "Grid size must be >0");
    assert(T>=0             && "Monte Carlo Temp must be >=0");
    assert(Nits >0          && "Number of samples must be >0");
    assert(Nburn>=0         && "Cannot have negative burn in itterations");
    assert(flips_per >=0    && "Number of flips per samples cannot be negative");

    // Make a grid
    grid::grid grid_instance(N, seed);
    grid_instance.reset(1);
    std::cout << "--------\n";
    std::cout << "Grid Generated. Initial State:\n";
    std::cout << "--------\n";
    grid_instance.print();
    std::cout << "--------\n";
    std::cout << "Initial Energy And Spin:\t" << grid_instance.energy() << "\t" << grid_instance.netspin() << "\n";

    // Perform MCMC run
    std::cout << "--------\n";
    std::cout << "T = "<< T <<"\n";
    std::cout << "Doing MCMC Run with "<< Nits <<" Samples and "<<Nburn<<" Burn-in \n";
    std::vector<std::vector<int>> outs = grid_monte(grid_instance, Nits, Nburn, T, flips_per, seed);

    //Do outputs
    std::cout << "--------\n";
    std::cout << "MCMC Done. Final State:\n";
    grid_instance.print();

    std::cout << "--------\n";
    std::cout << "Final Energy And Spin:\t" << grid_instance.energy() << "\t" << grid_instance.netspin() << "\n";
    grid_instance.calc_state();
    std::cout << "Final Energy And Spin:\t" << grid_instance.energy() << "\t" << grid_instance.netspin() << "\n";

    std::cout<<"Mean and std of Energy:\t"<< mean(outs[0]) <<"\t"<< sqrt_var(outs[0]) <<"\n";
    std::cout<<"Mean and std of Spin:\t"<< mean(outs[1]) <<"\t"<< sqrt_var(outs[1]) <<"\n";

    return 0;
}
