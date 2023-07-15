//====================================================
// Tuple.h
// Donny Chan, Billy Taj
// todo: refactor. rename vars


#ifndef TUPLE_H
#define TUPLE_H

#include <math.h>
#include <string>
#include <map>
#include <iostream>
using namespace std;
class Tuple{

    public:
        Tuple(int);
        Tuple(Tuple *);
        ~Tuple();

        void add(Tuple *);
        void scalarize(float);
        void copyValue(Tuple *);
        void negate();
        void display();
        void multByElement(Tuple *);
        void flush();
        void divByElement(Tuple *);
        int getNumElements();
        float sum();
        float * tuple;
        map<string, float> string_to_float_map; // this was designed to take the place of the float array.  
        
    private:
        int numElements;
};

#endif

