//====================================================
// functionlib.h
// Donny Chan, Billy Taj

#ifndef FUNCTIONLIB_H
#define FUNCTIONLIB_H

#include <sstream>
#include <math.h>
#include "ParameterManager.h"
#include "RandomNG.h"

static string floatToStr(float i){
    string s;
    stringstream ss(stringstream::in | stringstream::out);
    ss << i;
    ss >> s;
    return s;
}

static int roundToInt(float val){
    if(val < 0)
        return int(-val + 0.5) * -1;
    return int(val + 0.5);
}

static void output(string s){
    printf("%s", s.c_str());
}

static bool prob(float p){
     if(p <= 0.0) {return false;}
     if(p >= 1.0) {return true;}
     float result = (RandomNG::randFloat(0.0,1.0));
     return result <= p;
}

static int getI(string key){
    ParameterManager * pm = ParameterManager::getParameterManager();
    return pm->get_int(key);
}

#endif
