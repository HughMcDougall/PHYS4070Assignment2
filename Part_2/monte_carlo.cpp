//
// Created by hughm on 29/03/2023.
//

#include "monte_carlo.hpp"
#include <vector>
#include <cmath>
#include "grid.hpp"
#include "rand_utils.hpp"
#include <omp.h>

//-----------------------------------------------------------
// Utilities
std::vector<int> operator*=(std::vector<int> & a, const std::vector<int> & b){
    assert(a.size()==b.size() && "Tried to multiply two vectors with different lengths");
    for (int i=0; i<a.size(); i++){
        a[i]*=b[i];
    }
    return(a);
}
std::vector<int> operator*(std::vector<int> a, const std::vector<int> & b){return a*=b;}

double mean(const std::vector<int> &x){
    /// Gets the mean of a vector of ints
    double out=0;
    for (int a:x){
        out+=a;
    }
    out = double(out) / (double)x.size();
    return(out);
}

double var(const std::vector<int> &x){
    /// Gets the variance of a vector of ints
    double out=mean(x*x);
    out-=mean(x)*mean(x);
    return(out);
}

double sqrt_var(const std::vector<int> &x){
    /// Gets the standard deviation of a vector of ints
    double out=0;
    out = var(x);
    out = pow(out,0.5);
    return(out);
}

//-----------------------------------------------------------
// Monte Carlo Runs
std::vector<std::vector<int>> grid_monte(grid::grid & targ_grid, int Nits, int Nburn, double T, int flips_per_it, int seed){
    /// Runs a single MCMC run on a grid. Returns as vec of vec of ints: {energy, spin}

    assert(T>=0 && "Temperature cannot be negative");
    assert(Nits>0 && "Must have non-zero number of itterations");
    assert(targ_grid.N()>1 && "Grid size must be 2 or greater");
    assert(flips_per_it>=0 && "Need positive flips per itteration");

    //-----------------
    // Setup

    //Vectors to store results in
    std::vector<int> out_energies(Nits);
    std::vector<int> out_spins(Nits);

    //Indices for output and no. itterations
    int total_its = Nits+Nburn;
    int write_index = 0;
    bool burnin_finished  = false;
    if (flips_per_it == 0){flips_per_it = targ_grid.size();}

    //Objects for monte-carlo step
    double roll;
    int prop_energy;
    int prop_i, prop_j;
    bool use_step;
    double chance;

    //-----------------
    //Main monte-carlo loop
    for (int itt=0; itt<=total_its; itt++){
        // Do a sequence of flips to update the grid to a new sample state
        for(int j=0; j<flips_per_it; j++){
            //Choose a cell to flip
            prop_i= randbetween(0,targ_grid.N());
            prop_j= randbetween(0,targ_grid.N());

            //Do monte carlo step
            prop_energy = targ_grid.echange_from_flip(prop_i,prop_j);   //Calculate E change

            //Determine if proposal is accepted
            use_step=false;
            if (prop_energy<=0){
                use_step=true;
            } else if (T==0){
                use_step=false;
            } else{
                roll=randd();
                chance = exp(-1.0*(double)prop_energy / (double)T  );
                if( chance > roll){ use_step=true; } else {use_step=false;}
            }

            //If accepted, apply proposal
            if(use_step){
                targ_grid.flip(prop_i,prop_j, false, prop_energy);
            }
        }

        //Check if burn-in is done
        if (not burnin_finished) { if (itt > Nburn) { burnin_finished = true; }}

        //If out of burn-in phase, save to outputs
        if (burnin_finished){
            out_energies[write_index] = targ_grid.energy();
            out_spins[write_index] = targ_grid.netspin();
            write_index++;
        }
    }

    //-----------------
    //Output
    std::vector<std::vector<int>> out = {out_energies,out_spins};
    return(out);
}

//-----------------------------------------------------------
// Monte-Carlo with openMP
std::vector<std::vector<int>> monte_multithread(grid::grid & targ_grid, int Nits, int Nburn, double T){
    /// Runs a single MCMC run on a grid. Returns as vec of vec of ints: {energy, spin}

    assert(T>=0 && "Temperature cannot be negative");
    assert(Nits>0 && "Must have non-zero number of itterations");
    assert(targ_grid.N()>1 && "Grid size must be 2 or greater");


    //-----------------
    // Setup

    //Vectors to store results in
    std::vector<int> out_energies(Nits);
    std::vector<int> out_spins(Nits);

    //Indices for output and no. itterations
    int total_its = Nits+Nburn;
    int write_index = 0;
    bool burnin_finished  = false;

    // MISSINGNO - add thread specific rng

    //-----------------
    //Main monte-carlo loop, parallelized
    #pragma omp parallel default(none) shared(std::cout, targ_grid, Nits, Nburn, T, out_energies, out_spins, total_its, write_index, burnin_finished)
    {

    int ntid = omp_get_thread_num();
    #pragma omp critical
    {
        std::cout << "Hello from thread : " << ntid << std::endl;
    }

    #pragma omp master
    {
        std::cout << "Big Hello from thread : " << ntid << std::endl;
        std::cout << "Number of threads = " << omp_get_num_threads() << std::endl;
    }


    // Thread specific MCMC step variables
    double roll;
    int prop_energy;
    bool use_step;
    double chance;

    // Thread accumulators for energy and spin
    int thread_delta_M;
    int thread_delta_E;

    // Loop over all itterations
    for (int itt=0; itt<=total_its; itt++){
        thread_delta_M = 0;
        thread_delta_E = 0;

        // For both parities (odd and even)
        for(int parity = 0; parity<2; parity++){

            //For every other row...
            # pragma omp for
            for(int i=parity; i<targ_grid.N(); i+=2){

                // For every entry in that row...
                for(int j=0; j<targ_grid.N(); j++){

                    //Do monte carlo step
                    prop_energy = targ_grid.echange_from_flip(i, j);   //Calculate E change

                    //Determine if proposal is accepted
                    use_step=false;
                    if (prop_energy<=0){
                        use_step=true;
                    } else if (T==0){
                        use_step=false;
                    } else{
                        roll=randd();
                        chance = exp(-1.0*(double)prop_energy / (double)T  );
                        if( chance > roll){ use_step=true; } else {use_step=false;}
                    }

                    //If accepted, apply proposal
                    if(use_step){
                        thread_delta_M+= -2 * targ_grid.at(i,j);
                        thread_delta_E+= prop_energy;
                        targ_grid.flip(i,j, false, 0, false);
                    }

                } // element
            }

            // Barrier, all threads must finish odds / evens before moving to next parity
            #pragma omp critical
            {
            std::cout << "Waiting at barrier for thread : " << ntid << std::endl;
            }
            #pragma omp barrier
            #pragma omp critical
            {
            std::cout << "Passed barrier for thread : " << ntid << std::endl;
            }

        } // parity

        // After full sweep, have each thread apply its total energy and spin changes.
        #pragma omp critical
        {
            targ_grid.force_update(targ_grid.netspin()+thread_delta_M, targ_grid.energy()+thread_delta_E);
        } // update

        // Check if burn-in is done. If out of burn-in phase, write results
        // This task performed only by master thread
        # pragma omp master
        {
            if (not burnin_finished) { if (itt > Nburn) { burnin_finished = true; }}

            if (burnin_finished){
                out_energies[write_index] = targ_grid.energy();
                out_spins[write_index] = targ_grid.netspin();
                write_index++;
            }
        } // master

    } // itteration loop
    } // pragma omp parallel


    //-----------------
    //Output
    std::vector<std::vector<int>> out = {out_energies,out_spins};
    return(out);
}