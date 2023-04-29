//
// Created by hughm on 29/03/2023.
//

#include <vector>
#include "grid.hpp"

double mean(const std::vector<int> &x);
double var(const std::vector<int> &x);
double sqrt_var(const std::vector<int> &x);
std::vector<std::vector<int>> grid_monte(grid::grid & targ_grid, int Nits, int Nburn, double T=0, int flips_per_it=0, int seed =0);

#ifndef ASSIGNMENT_2_MONTE_CARLO_HPP
#define ASSIGNMENT_2_MONTE_CARLO_HPP

#endif //ASSIGNMENT_2_MONTE_CARLO_HPP
