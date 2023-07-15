//====================================================
// Matrix.cpp
// Donny Chan

#include "../inc/Matrix.h"
#include <stdio.h>
extern vector<string> small_molecule_id_list;

Matrix :: Matrix(Tuple * t){

//printf("Matrix Created:\t%d\n",(long)this);
    matrix.push_back(t);

    numRows = matrix.size();
    numCols = t->getNumElements();
}

Matrix :: Matrix(int rows, int cols){
    for(int r = 0 ; r < rows ; r ++){
        Tuple * t = new Tuple(cols);
        matrix.push_back(t);
    }

    numRows = rows;
    numCols = cols;
}

Matrix :: ~Matrix(){
//printf("Deleting Matrix:\t%d\n",(long)this);
    for(int i = 0 ; i < numRows ; i ++){
        delete matrix[i];
    }
}

void Matrix :: add(Matrix * m){
    if(numRows == m->getNumRows() && numCols == m->getNumCols()){
        for(int row = 0 ; row < numRows ; row++){
            matrix[row]->add((m->matrix)[row]);
        }
    }
}

void Matrix :: addAndClear(Matrix * m){
    if(numRows == m->getNumRows() && numCols == m->getNumCols()){
        for(int row = 0 ; row < numRows ; row++){
            matrix[row]->add((m->matrix)[row]);
            (m->matrix)[row]->flush();
        }
    }
}

void Matrix :: scalarize(float s){
    for(int row = 0 ; row < numRows ; row++){
        matrix[row]->scalarize(s);
    }
}

void Matrix :: addValueToRow(int row, Tuple * addRow){
    matrix[row]->add(addRow);
}

void Matrix :: copyValueToRow(int row, Tuple * copyRow){
    matrix[row]->copyValue(copyRow);
}

Tuple * Matrix :: collapse(){
    Tuple * t = new Tuple(numCols);
    for(int col = 0 ; col < numCols ; col ++){
        float sum = 0;
        for(int row = 0 ; row < numRows ; row ++){
            sum += matrix[row]->tuple[col];
        }
        t->tuple[col] = sum;
    }
    return t;
}

Tuple * Matrix :: getRow(int rowNum){
    if(rowNum < numRows){
        return matrix[rowNum];
    }
}

int Matrix :: getNumRows(){
    return numRows;
}

int Matrix :: getNumCols(){
    return numCols;
}

void Matrix :: display(){
    for(int i = 0 ; i < numRows ; i ++){
        Tuple * t = matrix[i];
        for(int j = 0 ; j < numCols ; j ++){
            printf("%f\t",t->tuple[j]);
        }
        printf("\n");
    }
}

void Matrix :: output_matrix(string filename){
    ofstream small_mol_count(filename);
    string delimiter = "";
    small_mol_count << "small_mol_summary" << "\n";
    for (string small_mol : small_molecule_id_list) {
        small_mol_count << delimiter << small_mol;
        delimiter = "\t";
    }
    delimiter = "";
    for(int i = 0 ; i < numRows ; i ++){
        Tuple * t = matrix[i];
        for(int j = 0 ; j < numCols ; j ++){
            small_mol_count << delimiter << t->tuple[j];
            delimiter = "\t";
        }
        small_mol_count << "\n";
        delimiter = "";
    }
    small_mol_count << endl;
    small_mol_count.close();
}

string Matrix :: displayToString(){
    stringstream mat("");
    for(int i = 0 ; i < numRows ; i ++){
        Tuple * t = matrix[i];
        for(int j = 0 ; j < numCols ; j ++){
            mat << to_string((t->tuple[j])) << "\t";
        }
        mat << endl;
    }
    return mat.str();
}

float Matrix :: get(int rowNum, int colNum){
    if(rowNum < numRows && colNum < numCols){
        return matrix[rowNum]->tuple[colNum];
    }
    return 0;
}

void Matrix :: appendRow(Tuple * appRow){

    if(numCols == appRow->getNumElements()){
        numRows ++;
        matrix.push_back(appRow);
    }
}
