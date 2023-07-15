//=======================================================
// Matrix.h -> 2021
// Donny Chan, Billy Taj
// Header for the Matrix class


#ifndef MATRIX_H
#define MATRIX_H

#include "Tuple.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <vector>
#include <string>
#include <sstream>
#include <fstream>

using namespace std;

class Matrix{

    public:
        Matrix(Tuple *);
        Matrix(int, int);
        ~Matrix();

        void add(Matrix *);
        void addAndClear(Matrix *);
        void scalarize(float);
        void addValueToRow(int, Tuple *);
        void copyValueToRow(int, Tuple *);
        void appendRow(Tuple *);
        float get(int, int);
        Tuple * collapse();

        Tuple *  getRow(int);

        int getNumRows();
        int getNumCols();
        void display();
        void output_matrix(string filename);
        string displayToString();

        vector <Tuple *> matrix;

    private:
        int numRows;
        int numCols;
};

#endif

