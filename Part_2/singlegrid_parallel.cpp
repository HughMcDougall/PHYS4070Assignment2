/// singlegrid.cpp

#include <iostream>
#include <omp.h>
#include "grid.hpp"
#include "rand_utils.hpp"
#include "monte_carlo.hpp"

int main(int argc, char *argv[]) {
    /// Runtime to execute and print a multi-threaded montecarlo run for an NxN grid. For Benchmarking
    /// Inputs are:
    /// Var         Type    Constraints Default Desc
    /// N           int     >0          8       Dimension of NxN Grid
    /// T           double  >=0         1       Temperature for Monte carlo simulation
    /// Nits        int     >0          1E4     Number of full sweeps of the grid / samples to draw
    /// seed        int                 0       Seed for random number generation in grid creation (redundant)

    //------------------------------------------------------------------
    // SIMULATION PARAMETERS (Defaults)
    int N=64;
    double T = 2.77*0.75;

    int Nits = (int)1E4;
    int Nburn = (int)1E3;
    int Nchains = 5;

    //------------------------------------------------------------------
    // INPUTS
    {
        if (argc > 1) { N = std::stoi(argv[1]); }
        if (argc > 2) { T = std::stod(argv[2]); }
        if (argc > 3) { Nits = std::stoi(argv[3]); }
        if (argc > 4) { Nburn = std::stoi(argv[4]); }
        if (argc > 5) { Nchains = std::stoi(argv[5]); }
    }

    //------------------------------------------------------------------
    // INTERPRET INPUTS
    assert(N>0              && "Grid size must be >0");
    assert(T>=0             && "Monte Carlo Temp must be >=0");
    assert(Nits >0          && "Number of samples must be >0");
    assert(Nburn>=0         && "Cannot have negative burn in itterations");

    //------------------------------------------------------------------
    // SET UP
    double time_start;
    double time_finish;
    std::vector<double> runtimes(Nchains);

    std::vector<std::vector<int>> chain;

    // Make a grid
    grid::grid grid_instance(N);
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
    std::cout << "Running multi-threaded MCMC runs with "<< Nchains <<" chains \n";

    for (int chain_no=0;chain_no<Nchains;chain_no++){
        time_start = omp_get_wtime();
        chain = monte_multithread(grid_instance, Nits, Nburn, T);
        time_finish = omp_get_wtime();
        runtimes[chain_no]=time_finish-time_start;
    }

    //Do outputs
    std::cout << "--------\n";
    std::cout << "MCMC Done. Final State:\n";
    grid_instance.print();

    std::cout << "--------\n";
    std::cout << "Final Energy And Spin:\t" << grid_instance.energy() << "\t" << grid_instance.netspin() << "\n";
    grid_instance.calc_state();

    std::cout<<"Mean and std of Energy:\t"<< mean(chain[0]) <<"\t"<< sqrt_var(chain[0]) <<"\n";
    std::cout<<"Mean and std of Spin:\t"<< mean(chain[1]) <<"\t"<< sqrt_var(chain[1]) <<"\n";

    std::cout << "--------\n";
    std::cout << "Time Elapsed (s):\n";
    for (int chain_no=0;chain_no<Nchains;chain_no++){
        std::cout << "\t"<< runtimes[chain_no] <<"\n";
    }



    return 0;
}
