//
// Used to generate the system inputs for a many body example. Just a fun test case
//

#ifndef ASSIGNMENT_2_MANYBODY_EXAMPLE_HPP
#define ASSIGNMENT_2_MANYBODY_EXAMPLE_HPP

#include <vector>
#include "params.hpp"
#include <cmath>
#include "vector_utils.hpp"
#include <random>

using vec = std::vector<double>;
using namespace params;

using sysvec  = std::vector<std::vector<double>>;

namespace manybody{
//----------------------------------------
vec polar_to_rect(double r, double theta){
    vec out = {cos(theta), sin(theta)};
    out = out*r;
    return(out);
}

std::mt19937 mt_generator(0);
std::uniform_real_distribution<double> uniform_rn(0.0, 1.0);


double randd(){
    double out = uniform_rn(mt_generator);
    return out;
}


inline std::vector<double> _make_M(int N_asteroids){
    std::vector<double> out(N_asteroids+2);
    out[0]=m_moon;
    out[1]=m_moon;
    for(int i=2; i<N_asteroids+2; i++){out[i] = m_moon/1000;}
    return out;
}

inline std::vector<bool> _make_LIVE(int N_asteroids){
    std::vector<bool> out(N_asteroids+2);
    out[0]=1;
    out[1]=1;
    for(int i=2; i<N_asteroids+2; i++){out[i] = 0;}
    return out;

}

inline sysvec _make_X(int N_asteroids){
    /// Generates the positions and velocities of the moons and asteroids
    int N = N_asteroids+2;
    sysvec X((N_asteroids+2)*2);

    //Positions of the two moons
    vec moon_pos_1 = {r_moon, 0};
    vec moon_vel_1 = {0, v_moon};
    vec moon_pos_2 = {-r_moon, 0};
    vec moon_vel_2 = {0, -v_moon};

    X[0] = moon_pos_1;
    X[1] = moon_pos_2;
    X[N+0] = moon_vel_1;
    X[N+1] = moon_vel_2;

    //Positions of the 4 lagrange points
    vec pos_1 = polar_to_rect(r_moon, M_PI / 180 * (60));;
    vec vel_1 = polar_to_rect(v_moon, M_PI / 180 * (60 + 90));;

    vec pos_2 = polar_to_rect(r_moon, M_PI / 180 * (-60));
    vec vel_2 = polar_to_rect(v_moon, M_PI / 180 * (-60 + 90));

    vec pos_3 = polar_to_rect(r_moon, M_PI / 180 * (120));;
    vec vel_3 = polar_to_rect(v_moon, M_PI / 180 * (120 + 90));;

    vec pos_4 = polar_to_rect(r_moon, M_PI / 180 * (240));
    vec vel_4 = polar_to_rect(v_moon, M_PI / 180 * (240 + 90));

    // For randomization in loop
    vec dr;
    vec dV;
    double dr_theta;
    double dr_mag;

    double dv_theta;
    double dv_mag;

    double flip;

    for(int i=2; i<N; i++){

        // Generate random deviations from lagrange point positions and velocities
        dr_mag      = randd() * 1;
        dr_theta    = randd() * 2*M_PI;

        dv_mag      = randd() * v_moon/20;
        dv_theta    = randd() * 2*M_PI;

        dr = polar_to_rect(dr_mag,dr_theta);
        dV = polar_to_rect(dv_mag,dv_theta);

        // Used to determine which L point each asteroid ends up at
        flip  = randd();

        // Generate Positions
        if (flip<0.25){
            X[i]    = dr + pos_1;
            X[N+i]  = dV + vel_1;
        } else if (flip<0.5){
            X[i]    = dr + pos_2;
            X[N+i]  = dV + vel_2;
        }else if (flip<0.75){
            X[i]    = dr + pos_3;
            X[N+i]  = dV + vel_3;
        }else{
            X[i]    = dr + pos_4;
            X[N+i]  = dV + vel_4;
        }
    }

    return X;

}

//----------------------------------------
// Actually generate the system
int N_asteroids = 20;

std::vector<double>     M = _make_M(N_asteroids);
std::vector<bool>    LIVE = _make_LIVE(N_asteroids);
sysvec                  X = _make_X(N_asteroids);

}
#endif //ASSIGNMENT_2_MANYBODY_EXAMPLE_HPP
