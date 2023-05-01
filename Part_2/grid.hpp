#ifndef ASSIGNMENT_2_GRID_H
#define ASSIGNMENT_2_GRID_H

#pragma once
#include <random>
#include <cassert>
#include <vector>
#include <iostream>

#include "rand_utils.hpp"


namespace grid{
    /// Object class for for 2D NxN grid of up/down spins
    /// Contains methods for editing states and calculating spin / energy

    //-----------------------------------------------------------
    //Utils

    inline int smartmod(int i,int N){
        /// variant of modulo operator that works on negative indices
        while (i<0){i+=N;}
        return i % N;
    }

    using unt = std::size_t;
    //-----------------------------------------------------------

    class grid{
        //-----------------------------------------------------------
        private:
            std::vector<int> _datavec;
            unt _N;
            unt _size;
            int _energy;
            int _netspin;

        //-----------------------------------------------------------
        public:
            //==============
            //Setup
            grid(int N, int seed=0): _N(N), _size(N*N){
                /// Constructor function. By default, grid is random.
                _energy = 0;
                _netspin= 0;
                _datavec.resize(_size);
                shuffle(seed);
            }

            void shuffle(int seed=0){
                ///Randomizes the grid with srand(seed)
                srand(seed);

                for (int i=0; i<_N; i++){
                    for (int j=0; j<_N; j++){
                        at(i,j)=rand_plusminus();
                    }
                }
                calc_state();
            }

            void reset(int mode){
                assert(mode == 1 or mode == -1 && "in grid.reset(), mode must be 1 (up) or -1 (down)");
                ///Randomizes the grid with srand(seed)

                for (int i=0; i<_N; i++){
                    for (int j=0; j<_N; j++){
                        at(i,j)=mode;
                    }
                }
                calc_state();
            }

            void calc_state(){
                /// Does a full evaluation of the systems hamiltonian and updates _energy.
                /// For use after setting up or randomizing a grid
                _energy = 0;
                _netspin= 0;

                int e_per_site = 0;

                for (int i=0; i<_N; i++){
                    for (int j=0; j<_N; j++){
                        e_per_site  = (at(i+1,j)+at(i-1,j+1));
                        e_per_site *= -at(i,j);
                        _energy    += e_per_site;

                        _netspin+= at(i,j);
                    }
                }
                /// _energy*=4;

            }

            void force_update(int M, int E){
                /// Forces an update to the grid's energy and spin.
                /// For use in multithreading to do bulk application of state changes. Not for normal use.
                _netspin=M;
                _energy=E;
            }

            //==============
            //Gets
            int& at(unt i, unt j)        {return _datavec[smartmod(i,_N) * _N + smartmod(j,_N)];}  //For assigning
            int  at(unt i, unt j) const  {return _datavec[smartmod(i,_N) * _N + smartmod(j,_N)];}  //Const version
            int &operator()(std::size_t i, std::size_t j) { return at(i, j); }         //Quick format of .at
            int operator()(std::size_t i, std::size_t j) const { return at(i, j); }    //Const Version

            unt N() const {return _N;}                                                 //Shield for _N
            unt size() const {return _size;}                                           //Shield for _size=_N*_N
            int energy() const{return _energy;}                                        //Shield for _energy
            int netspin() const{return _netspin;}                                      //Shield for _netspin
            int *data() {return _datavec.data();}                                      //Shortcut to data array

            //==============
            //State Changes
            void flip(unt i, unt j, bool calc_energy = true, int energy_change = 0, bool change_state=true){
                /// Flips grid location i,j and updates grid's energy and spin
                if (calc_energy)    {energy_change  =   echange_from_flip(i,j);}

                // Update grid properties. Skip E and M updates if flagged, to avoid double write in openMP
                at(i,j) *= -1;
                if (change_state){
                    _energy  += energy_change;
                    _netspin += at(i,j)*2; // double as we over-write the old spin
                }
            }

            int echange_from_flip(unt i, unt j){
                /// Returns the energy change that will result from flipping grid location i,j
                int out = 0;

                out+=at(i+1,j);
                out+=at(i,j+1);
                out+=at(i-1,j);
                out+=at(i,j-1);

                out*=-at(i,j);
                out*=-2;

                /// out*=4;

                return out;
            }

            //==============
            //Outputs
            void print(std::ostream& stream = std::cout){
                /// Outputs an ascii render  of the grid to 'stream'. Defaults to std::cout
                for (unt i=0; i < _N; i++){
                    for (unt j=0; j < _N; j++){
                        if (at(i,j) == 1){stream << 0 << " ";} else if (at(i,j) == -1){stream << "-" << " ";} else stream << "!" << " ";
                        }
                    std::cout << "\n";
                    }
                }
        //-----------------------------------------------------------
        }; //Class

    } //Namespace

#endif //ASSIGNMENT_2_GRID_H
