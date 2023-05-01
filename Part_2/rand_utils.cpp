//
// Contains functions for easily generating random numbers. Use this for non-parallelized MCMC runs,
// As it's easier to pass an SRAND seed for a particular run
//

#include <random>
#include <cassert>
#include "rand_utils.hpp"

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