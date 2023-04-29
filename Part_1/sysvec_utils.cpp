//
// Utilities to make working with system vectors easier. Based on vector_utils.cpp from assignment 1
//

#include "sysvec_utils.hpp"

#include <vector>
#include <functional>
#include <cassert>
#include <fstream>

#include "vector_utils.hpp"

using sysvec  = std::vector<std::vector<double>>;
using diff_func  = std::function<sysvec(sysvec)>;


//================================================================
// UTILITIES

void print_to_stream(const sysvec &X, std::ostream &output){
    /// Prints all entries in a system vector to 'output', which defaults to std::cout
    int N = X.size();
    int dim = X[0].size();

    for(int i=0; i<N;i++){
        for(int j=0; j<dim;j++){
            output << X[i][j] <<"\t";
        }
    }
    output << "\n";
}
sysvec vconcat(const sysvec  & V1, const sysvec & V2){
    ///Concatenates two sysvecs
    int N1 = V1.size();
    int N2 = V2.size();
    sysvec out (N1+N2);

    for (int i=0; i<N1; i++){
        out[i]=V1[i];
    }
    for (int i=N1; i<N2; i++){
        out[i]=V2[i-N1];
    }

    return out;
}
//================================================================
// OVERLOADS
//Overload vector operations to make direct products easier

//V-V Multiplication
sysvec operator*=(sysvec & a, const sysvec & b){
    assert(a.size()==b.size() && "Tried to multiply two vectors with different lengths");
    for (int i=0; i<a.size(); i++){
        a[i]*=b[i];
    }
    return(a);
}
sysvec operator*(sysvec a, const sysvec & b){return a*=b;}

//V-V Division
sysvec operator/=(sysvec & a, const sysvec & b){
    assert(a.size()==b.size() && "Tried to divide two vectors with different lengths");
    for (int i=0; i<a.size(); i++){
        a[i]/=b[i];
    }
    return(a);
}
sysvec operator/(sysvec a, const sysvec & b){return a/=b;}

//V-V Addition
sysvec operator+=(sysvec & a, const sysvec & b){
    assert(a.size()==b.size() && "Tried to add two vectors with different lengths");
    for (int i=0; i<a.size(); i++){
        a[i]+=b[i];
    }
    return(a);
}
sysvec operator+(sysvec a, const sysvec & b){return a+=b;}

//V-V Subtraction
sysvec operator-=(sysvec & a, const sysvec & b){
    assert(a.size()==b.size() && "Tried to subtract two vectors with different lengths");
    for (int i=0; i<a.size(); i++){
        a[i]-=b[i];
    }
    return(a);
}
sysvec operator-(sysvec a, const sysvec & b){return a-=b;}

//V-D Multiplication
sysvec operator*=(sysvec & v, const double & a){
    for (int i=0; i<v.size(); i++){
        v[i]*=a;
    }
    return(v);
}
sysvec operator*(sysvec v, const double & a){return(v*=a);}
sysvec operator*(const double & a, sysvec v){return(v*a);}

//V-D Division
sysvec operator/=(sysvec & v, const double & a){
    for (int i=0; i<v.size(); i++){
        v[i]/=a;
    }
    return(v);
}
sysvec operator/(sysvec v, const double & a){return(v/=a);}
sysvec operator/(const double & a, sysvec v){return(v/a);}

//V-D Addition
sysvec operator+=(sysvec & v, const double & a){
    for (int i=0; i<v.size(); i++){
        v[i]+=a;
    }
    return(v);
}
sysvec operator+(sysvec v, const double & a){return(v+=a);}
sysvec operator+(const double & a, sysvec v){return(v+a);}

//V-D Subtraction
sysvec operator-=(sysvec & v, const double & a){return(v+=a*-1);}
sysvec operator-(sysvec v, const double & a){return(v-=a);}
sysvec operator-(const double & a, sysvec v){return(v-a);}