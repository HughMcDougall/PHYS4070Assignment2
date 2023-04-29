//
// This file contains all physically meaningful functions for PHYS4070 Assignment 2 part 1
//

#include "forces_and_integrators.hpp"

#include <vector>
#include <functional>
#include <cmath>
#include <cassert>
#include <iostream>

#include "params.hpp"
#include "vector_utils.hpp"
#include "sysvec_utils.hpp"

using vec = std::vector<double>;
using sysvec  = std::vector<std::vector<double>>;
using diff_func  = std::function<sysvec(sysvec)>;
using namespace params;

//================================================================
// FORCE CALCULATORS
/// planet only force update
sysvec f_planet(const sysvec &X){
    /// Returns the time derivative of state vector X under the influence of central planet gravity
    int dim = X[0].size();
    int N = X.size() / 2;

    for (int i=0; i<X.size(); i++){assert(X[i].size()==dim && "Tried to calculate forces in f_planet() with different dimensions in sysvec entries");}

    sysvec out = X * 0;

    // Fill out first 'N' elements with velocities
    for(int i=0; i<N; i++){
        out[i] = X[N+i];
    }

    // Fill out remaining elements with gravity acceleration terms
    for(int i=0; i<N; i++){
        out[N+i] = -1*m_planet * X[i] / pow(vnorm(X[i]),3);
    }

    return out;
}

/// n-body planet update
sysvec f_nbody(const sysvec &X, const std::vector<double> &M, const std::vector<bool> &LIVE){
    /// Returns the time derivative of state vector X under the influence of central planet gravity
    /// Takes system vector X, double vector of masses M and bool vector "LIVE"
    /// Must be collapsed to f=f(X) only before feeding to runge_int or euler_int

    int dim = X[0].size();
    int N = X.size() / 2;

    for (int i=0; i<X.size(); i++){assert(X[i].size()==dim && "Tried to calculate forces in f_nbody() with different dimensions in sysvec entries");}

    // Begin by getting update under gravity only. Impact of velocity already included in here.
    sysvec out = f_planet(X);

    // Create an empty force vector
    sysvec forces(N);
    std::vector<double> force(dim);
    std::vector<double> dr(dim);
    for (int i=0; i<N; i++){forces[i].resize(dim);}

    // Loop over all objects, but not repeating index pairs
    for (int i=0; i<N; i++){
        for (int j=0; j<i; j++){

            // Only calculate the force if objects are different and at least one of them is massive enough to care about
            if ((LIVE[i] or LIVE[j]) && i!=j){
                // Calculate force _from_ j acting on i
                dr = X[i] - X[j];
                force = -M[i]*M[j] * dr / pow(vnorm(dr),3);

                // Add resulting acceleration to both bodies. To prevent drift, only do so if the other body is massive
                if (LIVE[j]) { forces[i] += force / M[i]; }
                if (LIVE[i]) { forces[j] -= force / M[j]; }
            }
        }
    }

    //Write these accelerations to the output vector
    for (int i=0; i<N; i++){
        out[N+i]+=forces[i];
    }

    return out;
}

//================================================================
// INTEGRATORS
sysvec euler_step(sysvec &X, const diff_func &f, double dt){
    /// Updates a system vector using euler step
    /// For debug purposes only
    sysvec dx = f(X) * dt;
    return X + dx;
}

sysvec runge_step(sysvec &X, const diff_func &f, double dt){
    /// Updates a system vector using the fourth order runge kutta method

    sysvec k1 = f(X);
    sysvec k2 = f(X + dt/2 * k1);
    sysvec k3 = f(X + dt/2 * k2);
    sysvec k4 = f(X + dt   * k3);

    return X + (k1/6 + k2/3 + k3/3 + k4/6) * dt;
}

//================================================================
// ORBITAL KICK

vec orbital_kick(const vec & dr, const vec & v0, const vec & vmoon){
    /// Calculates the impulse required to achieve a stable circular orbit about the moon

    //magntiude of orbital velocity relative to moon
    double orb_vel = pow(m_moon / vnorm(dr),0.5);

    //Directions
    vec r_dir = dr/vnorm(dr);
    vec tangent_dir_1 = {r_dir[1], -1 * r_dir[0]};
    vec tangent_dir_2 = -1 * tangent_dir_1;

    //Calculates prograde and retrograde kicks
    vec kick_1 = vmoon + orb_vel * tangent_dir_1 - v0;
    vec kick_2 = vmoon + orb_vel * tangent_dir_2 - v0;

    //Return the smaller of the two
    if (vnorm(kick_1)<vnorm(kick_2)){
        return(kick_1);
    }
    else{
        return(kick_2);
    }

}