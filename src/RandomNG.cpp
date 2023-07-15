//===================================================
// RandomNG.cpp
// Donny Chan, Billy Taj

#include "../inc/RandomNG.h"


using namespace std;

int RandomNG :: inext;
int RandomNG :: inextp;
long RandomNG :: ma[56];
int RandomNG :: iff = 0;
long seed;
// potentially reseed later
// long seed = 3777855842;
mt19937_64 gen;

RandomNG::RandomNG(){
    
}


// Performs mersenne twister to select 1 item from the provided map
int RandomNG :: custom_discrete_distribution(vector<float> input_distribution) { 
    discrete_distribution<> reaction_weights(input_distribution.begin(),input_distribution.end());
    return reaction_weights(gen);
}

// gaussian sampling for particle movement, usually mean of 0 and SD of root(2 * difc * timescale), which is rms step length of particles in 1D.
float RandomNG::custom_normal_distribution(float mean, float mobility) {
    normal_distribution<float> distribution(mean, mobility);
    return(distribution(gen));
}

// generate an unsigned long int using mersenne twister
unsigned long RandomNG::randLong() {
    uniform_int_distribution<unsigned long> distribution(1, ULONG_MAX);
    unsigned long new_long = distribution(gen);
    while(new_long == 0) {
        new_long = distribution(gen);
    }
    return new_long;
}

// Random integer value between min and max (min <= value <= max)
int RandomNG :: randInt(int min, int max){
    if(min > max)    swap(min, max);
    uniform_int_distribution<int> distribution(min, max);
    int output = distribution(gen);

    return output;
}

// generate vector of ints sampled inclusively between min and max
vector<int> RandomNG::randIntVec(int min, int max, int size) {
    if(min > max)    swap(min, max);
    uniform_int_distribution<int> distribution(min, max);
    vector<int> output_vec;
    for(int counter = 0; counter < size; counter++) {
        output_vec.push_back(distribution(gen));
    }
    return output_vec;
}

// Random float value between min and max (min < value <= max)
float RandomNG :: randFloat(float min, float max){
    if(min > max)    swap(min, max);
    uniform_real_distribution<float> distribution(min, max);
    float output = distribution(gen);

    return output;

}

// Negative Expontential
float RandomNG :: randNegExp(float val){
    return -val * log(orand()) ;
}

// Erlang
float RandomNG :: randErlang(float x, float s){

    int i, k;
    float z;

    z = x / s;
    k = (int) (z * z);
    for (i = 0, z = 1.0; i < k; i++)
        z *= orand();
    return -(x/k) * log(z);
}

// Normal Distribution with mean and standard deviation
float RandomNG :: randNormal(float mean,float stdev){
    float v1, v2, w, z1;

    do {
        v1 = 2.0 * orand() - 1.0;
        v2 = 2.0 * orand() - 1.0;
        w = v1 * v1 + v2 * v2;
    } while ( w >= 1.0 );

    w = sqrt( (-2.0 * log( w ) ) / w );
    z1 = v1 * w;

    return (mean + z1 * stdev);
}


// Get next Random float between 0.0 and 1.0
float RandomNG :: orand()
{

    long mj,mk;
    int i,ii,k;

    if(iff == 0){
        long idum = getSeed();

        iff=1;
        mj=MSEED-(idum < 0 ? -idum : idum);
        mj = abs(mj);
        mj %= MBIG;
        ma[55]=mj;
        mk=1;
        for(int i = 1 ; i <= 54 ; i++){
            int ii = (21 * i) % 55;
            ma[ii] = mk;
            mk = mj - mk;
            if (mk < MZ) mk += MBIG;
            mj = ma[ii];
        }
        for (int k = 1 ; k <= 4 ; k++)
            for(int i = 1 ; i <= 55 ; i++){
                ma[i] -= ma[1+(i+30) % 55];
                if (ma[i] < MZ) ma[i] += MBIG;
            }

        inext=0;
        inextp=31;
    }

    if(++inext == 56) inext=1;
    if(++inextp == 56) inextp=1;

    mj = ma[inext]-ma[inextp];

    if(mj < MZ) mj += MBIG;

    ma[inext]=mj;
    return mj*(1.0/MBIG);
}

void RandomNG :: setSeed(long input_seed) {
    if(!input_seed == 0) {
        seed = input_seed;
        gen.seed(RandomNG::getSeed());
    } else {
        seed = static_cast<unsigned>(time(0)) * getpid();
        gen.seed(RandomNG::getSeed());
    }
    cout << "starting Seed: " << RandomNG::getSeed() << endl;
}

// Set the seed for the random sequence
long RandomNG :: getSeed(){
    // printf("Seed\t%d\n", seed);
    return seed;
}

mt19937_64 & RandomNG::get_generator() {
    return gen;
}

vector<int> RandomNG::create_shuffled_index_vector(int vec_size) {
    vector<int> vec(vec_size);
    iota(begin(vec), end(vec), 0);
    shuffle(vec.begin(), vec.end(), gen);
    return(vec);
}
