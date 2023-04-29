/// temp_sweep.cpp

#include <iostream>
#include <string>
#include <fstream>

#include "grid.hpp"
#include "rand_utils.hpp"
#include "vector_utils.hpp"
#include "monte_carlo.hpp"


int main(int argc, char *argv[]) {
    /// Runs a temp sweep in the range T_min to T_max
    /// Inputs are:
    /// Var         Type    Constraints Default Desc
    /// N           int     >0          16                          Dimension of NxN Grid
    /// T_min       double  >=0         0.0                         Min Temp in Sweep
    /// T_max       double  >0          5.0                         Max Temp in sweep
    /// N_temps     int     >0          32                          Number of temps to test
    /// output_dir  str                 ./results/temp_sweep.dat    Filename of output
    /// Nits        int     >0          1E4                         Number of full sweeps of the grid / samples to draw
    /// Nburn       int     >=0         1E3                         Number of burn-in itterations in monte carlo
    /// Nchains       int   >0          1E3                         Number of independent monte carlo chains to run
    /// flips_per   int     >0          0(NxN)                      Number of bit flips per 'sample' If set to zero, will use flips_per=NxN
    /// seed        int                 0                           Seed for random number generation

    //------------------------------------------------------------------
    // SIMULATION PARAMETERS (Defaults)
    int    N = 16;

    double T_min = 0.0;
    double T_max = 5.0;
    int N_temps = 32;


    int Nits = (int)1E4;
    int Nburn = (int)1E3;
    int Nchains = 8;
    std::string output_dir = "./results/temp_sweep.dat";

    int flips_per = 0;
    int seed = 0;

    //------------------------------------------------------------------
    // INPUTS
    {
        if (argc > 1) { N = std::stoi(argv[1]); }

        if (argc > 2) { T_min = std::stod(argv[2]); }
        if (argc > 3) { T_max = std::stod(argv[3]); }
        if (argc > 4) { N_temps= std::stoi(argv[4]); }
        if (argc > 5) { output_dir= argv[5]; }

        if (argc > 6) { Nits = std::stoi(argv[6]); }
        if (argc > 7) { Nburn = std::stoi(argv[7]); }
        if (argc > 8) { Nchains = std::stoi(argv[8]); }

        if (argc > 9) { flips_per = std::stoi(argv[9]); }
        if (argc > 10) { seed = std::stoi(argv[10]); }
    }

    //------------------------------------------------------------------
    // INTERPRET INPUTS
    assert(N>0                      && "Grid size must be >0");
    assert(T_min<T_max && T_max>=0  && "Monte Carlo Temps must be >=0");
    assert(Nits >0                  && "Number of samples must be >0");
    assert(N_temps>2                && "Number of temps must be >2");
    assert(Nburn>=0                 && "Cannot have negative burn in itterations");
    assert(flips_per >=0            && "Number of flips per samples cannot be negative");

    //------------------------------------------------------------------
    // MONTE CARLO RUN
    std::cout << "Doing Temp Seep for Square L = " << N <<" grid with "<< N_temps<<" temps between T = "<< T_min<< " and T = "<< T_max<< "\n";
    std::cout << "Doing set up for MCMC run...\n";
    std::cout << "Doing MCMC Run with "<< Nits <<" Samples and "<<Nburn<<" Burn-in \n";
    std::cout << "Each sample pulled after "<<flips_per<<" bit flips\n";

    // Create grid of temps
    std::vector<double> TEMPS = make_rgrid(T_min,T_max,N_temps);
    // Make grids
    grid::grid grid_instance(N, seed);
    grid_instance.reset(rand_plusminus()); // Set to pure up or pure down to avoid local 50/50 equilibrium
    // Prepare output file
    std::ofstream output(output_dir);


    std::cout<< "Doing sweep over temps ...\n";

    double T;

    for (int k=0; k<Nchains;k++){

        if(Nchains>1){  std::cout<<"Chain number "<<k<<"\n"; grid_instance.reset(rand_plusminus());	}

        //Do a sweep over all temps
        for (int i=0;i<N_temps;i++){

            // Get temp
            T = TEMPS[i];

            // Perform MCMC run
            std::cout<<"Doing MCMC run "<< k*N_temps+i<<" of " << Nchains*N_temps<<" at temp T="<<T<<":\n";
            std::vector<std::vector<int>> chain = grid_monte(grid_instance, Nits, Nburn, T, flips_per, seed);

            // Save results
            std::cout<<"\tMCMC done, saving outputs...\n";

            output << T << "\t" << mean(chain[0]) << "\t" << var(chain[0]) <<"\t" << mean(chain[1]) << "\t" << var(chain[1]);
            if(Nchains>1){  output<<"\t"<<k;}
            output << "\n";

            std::cout<<"\tOutputs done.\n";
        }
    }

    //------------------------------------------------------------------
    // FINISH

    std::cout<<"Temp Sweep complete \n";
    output.close();
    return 0;
}
