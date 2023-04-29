//
// Created by hughm on 19/04/2023.
//

#ifndef ASSIGNMENT_2_SYSVEC_UTILS_HPP
#define ASSIGNMENT_2_SYSVEC_UTILS_HPP

#include <vector>
#include <functional>
#include <cassert>
#include <fstream>
#include <iostream>

#include "vector_utils.hpp"

using sysvec  = std::vector<std::vector<double>>;
using diff_func  = std::function<sysvec(sysvec)>;

//================================================================
// UTILITIES

void print_to_stream(const sysvec &X, std::ostream &output = std::cout);
sysvec vconcat(const sysvec  & V1, const sysvec & V2);

//Overload vector operations to make direct products easier
//V-V Multiplication
sysvec operator*=(sysvec & a, const sysvec & b);
sysvec operator*(sysvec a, const sysvec & b);

//V-V Division
sysvec operator/=(sysvec & a, const sysvec & b);
sysvec operator/(sysvec a, const sysvec & b);

//V-V Addition
sysvec operator+=(sysvec & a, const sysvec & b);
sysvec operator+(sysvec a, const sysvec & b);

//V-V Subtraction
sysvec operator-=(sysvec & a, const sysvec & b);
sysvec operator-(sysvec a, const sysvec & b);

//V-D Multiplication
sysvec operator*=(sysvec & v, const double & a);
sysvec operator*(sysvec v, const double & a);
sysvec operator*(const double & a, sysvec v);

//V-D Division
sysvec operator/=(sysvec & v, const double & a);
sysvec operator/(sysvec v, const double & a);
sysvec operator/(const double & a, sysvec v);

//V-D Addition
sysvec operator+=(sysvec & v, const double & a);
sysvec operator+(sysvec v, const double & a);
sysvec operator+(const double & a, sysvec v);

//V-D Subtraction
sysvec operator-=(sysvec & v, const double & a);
sysvec operator-(sysvec v, const double & a);
sysvec operator-(const double & a, sysvec v);

#endif //ASSIGNMENT_2_SYSVEC_UTILS_HPP

