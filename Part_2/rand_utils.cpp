//
// Created by hughm on 29/03/2023.
//

#include <random>
#include <cassert>
#include "rand_utils.hpp"

std::mt19937 mt_generator(0);
std::uniform_real_distribution<double> uniform_rn(0.0, 1.0);

double randd() {
    double out = (double)rand() / (double)RAND_MAX;
    return out;
}

int rand_plusminus(){
    return (rand() % 2)*2 -1;
}

int randbetween(int a, int b){
    assert( b>a && "in randbetween on interval [a,b), b must be > a");
    return (rand() % (b-a)) + a;
}