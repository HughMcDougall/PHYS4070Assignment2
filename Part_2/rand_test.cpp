//
// Created by hughm on 23/04/2023.
//


#include <iostream>
#include "grid.hpp"
#include "rand_utils.hpp"
#include "monte_carlo.hpp"

int main(){
    double T = 2.27*8;
    int prop_i;
    int prop_j;
    int prop_energy;
    bool use_step;
    double marge;
    double roll;
    double chance;

    grid::grid targ_grid(8,0);
    targ_grid.reset(1);


    for (int i=0; i<2000;i++){

        prop_i= randbetween(0,targ_grid.N());
        prop_j= randbetween(0,targ_grid.N());

        //Do monte carlo step
        prop_energy = targ_grid.echange_from_flip(prop_i,prop_j);   //Calculate E change

        //Determine if proposal is accepted
        if (prop_energy<=0){
            use_step=true;
        } else if (T==0){
            use_step=false;
        } else{
            roll=randd();
            marge = -1.0*(double)prop_energy / double(T);
            chance = exp(marge);
            if( chance - roll > 0){ use_step=true;} else{use_step=false;}
        }

        if (i % (2000/20)==0){
            std::cout<<"At index\t"<<prop_i<<'\t'<<prop_j<<"\n";
            std::cout<<"Echange: "<<prop_energy<<", at T ="<<T<<"\n";
            std::cout<<"Chance:\t"<<chance<<"\tagainst\t"<<roll<<"\n";
            std::cout<<"using step: "<<use_step<<"\n";
        }


        //If accepted, apply proposal
        if(use_step){
            targ_grid.flip(prop_i,prop_j, false, prop_energy);
        }
    }

    targ_grid.print();




    return 0;
}