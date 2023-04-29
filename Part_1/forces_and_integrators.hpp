//
// Created by hughm on 19/04/2023.
//

#ifndef ASSIGNMENT_2_FORCES_AND_INTEGRATORS_HPP
#define ASSIGNMENT_2_FORCES_AND_INTEGRATORS_HPP

#include <vector>
#include <functional>

#include "params.hpp"
#include "vector_utils.hpp"
#include "sysvec_utils.hpp"

using vec = std::vector<double>;
using sysvec  = std::vector<std::vector<double>>;
using diff_func  = std::function<sysvec(sysvec)>;

sysvec f_planet(const sysvec &X);
sysvec f_nbody(const sysvec &X, const std::vector<double> &M, const std::vector<bool> &LIVE);
sysvec euler_step(sysvec &X, const diff_func &f, double dt);
sysvec runge_step(sysvec &X, const diff_func &f, double dt);

vec orbital_kick(const vec & dr, const vec & v0, const vec & vmoon);

#endif //ASSIGNMENT_2_FORCES_AND_INTEGRATORS_HPP
