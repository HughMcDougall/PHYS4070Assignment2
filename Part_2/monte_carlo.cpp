//
// Created by hughm on 29/03/2023.
//

#include "monte_carlo.hpp"
#include <vector>
#include <cmath>
#include "grid.hpp"
#include "rand_utils.hpp"

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

    //Objects for monte-carlo step
    double roll;
    int prop_energy;
    int prop_i, prop_j;
    bool use_step;
    if (flips_per_it == 0){flips_per_it = targ_grid.size();}
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
