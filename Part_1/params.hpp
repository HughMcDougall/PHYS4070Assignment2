//
// Parameters (mass, radius) and best estimates of required values
// Stored here to avoid doubling up in runtimes
// HM Apr 23
//
#pragma once
#ifndef ASSIGNMENT_2_PARAMS_H
#define ASSIGNMENT_2_PARAMS_H

#endif //ASSIGNMENT_2_PARAMS_H
#include <cmath>

namespace params{
    /// Masses
    inline double m_planet = 1;
    inline double m_moon = 0.1;
    inline double m_proj = 1E-24;

    /// radii
    inline double R_planet = 1;
    inline double R_moon = 0.25;
    inline double R_proj = 1E-6;

    /// Orbital properties, moon about planet
    inline double r_moon = 19*R_planet;
    inline double v_moon  = pow(r_moon,-0.5);

    inline double v_proj_launch = pow(2*(1/R_planet - 1/(r_moon - R_moon)),0.5);

    /// Orbital properties, projectile about moon
    inline double r_proj_moon = 2*R_moon;
    inline double v_proj_moon  = pow(m_moon/ r_proj_moon,0.5);

}