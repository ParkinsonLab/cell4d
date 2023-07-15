//===============================================
// RandomNG.h
// Donny Chan, Billy Taj

#ifndef RANDOMNG_H
#define RANDOMNG_H

#include <algorithm>
#include <stdlib.h>
#include <cstdlib>
#include <stdio.h>
#include <ctime>
#include <iomanip>
#include <iostream>
#include <vector>
#include <math.h>
#include <unistd.h>
#include <map>
#include <climits>
#include <random>


#define MBIG 1000000000
#define MSEED 161803398
#define MZ 0
#define FAC (1.0/MBIG)

using namespace std;

class RandomNG { 
    public:

        static int inext,inextp;
        static long ma[56];
        static int iff;

        static int custom_discrete_distribution(vector<float>);
        static float custom_normal_distribution(float, float);
        static unsigned long randLong();
        static int randInt(int, int);
        static float randFloat(float, float);
        static float randNormal(float, float);
        static float randNegExp(float);
        static float randErlang(float, float);

        static vector<int> create_shuffled_index_vector(int vec_size);
        static vector<int> randIntVec(int min, int max, int size);

        static long getSeed();
        static void setSeed(long);

        static mt19937_64 & get_generator();


    protected:

        RandomNG();

    private:

        static float orand();

}; 

#endif
