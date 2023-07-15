//=========================================================
// Tuple.cpp
// Donny Chan


#include "../inc/Tuple.h"
#include <stdio.h>

Tuple :: Tuple(int num){
    numElements = num;

    tuple = new float [numElements];

    flush();
}

Tuple :: Tuple(Tuple * t){
    numElements = t->getNumElements();
    tuple = new float[numElements];
    for(int i = 0 ; i < numElements ; i ++){
        tuple[i] = t->tuple[i];
    }
}

Tuple :: ~Tuple(){
    // float array now gone.  it's now a string map with a count
    delete [] tuple;
}

void Tuple :: add(Tuple * t){
    if(numElements == t->getNumElements()){
        for(int i = 0 ; i < numElements ; i++){
            tuple[i] += (t->tuple)[i];
        }
    }
}

float Tuple :: sum(){
    float value = 0;
    for(int i = 0 ; i < numElements ; i ++)
        value += tuple[i];
    return value;
}

void Tuple :: scalarize(float s){
    for(int i = 0 ; i < numElements ; i++){
        tuple[i] *= s;
    }
}

void Tuple :: multByElement(Tuple * t){
    if(numElements == t->getNumElements())
    for(int i = 0 ; i < numElements ; i ++)
        tuple[i] *= t->tuple[i];
}

void Tuple :: divByElement(Tuple * t){
    if(numElements == t->getNumElements())
    for(int i = 0 ; i < numElements ; i ++)
            tuple[i] = tuple[i]/(t->tuple[i]);
}


void Tuple :: copyValue(Tuple * t){
    if(numElements == t->getNumElements()){
        for(int i = 0 ; i < numElements ; i ++){
            tuple[i] = t->tuple[i];
        }
    }
}

void Tuple :: negate(){
    for(int i = 0 ; i < numElements ; i ++)
        tuple[i] = ((-1) * tuple[i]);
}

int Tuple :: getNumElements(){
    return numElements;
}

void Tuple :: display(){
	if (numElements > 0) {
	    for(int i = 0 ; i < numElements ; i ++){
	        cout << tuple[i] << "\t";
	    }
	    cout << endl;
	}
}

void Tuple :: flush(){
    for(int i = 0 ; i < numElements ; i ++){
        tuple[i] = 0;
    }
}
