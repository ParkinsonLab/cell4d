//===================================================
// Vector3D.cpp
// Donny Chan, Billy Taj

#include "../inc/Vector3D.h"
#include <stdio.h>

// Default Constructor
Vector3D :: Vector3D(){
    x = y = z = 0;
}

// Parametized Constructor
Vector3D :: Vector3D(float xVal, float yVal, float zVal){
    x = xVal; y = yVal; z = zVal;
}

// Parametized Constructor
Vector3D :: Vector3D(const Vector3D & v){
    x = v.x ; y = v.y; z = v.z;
}

// Take the cross product and return the vector
Vector3D Vector3D :: cross(const Vector3D & v){
    float c_x = y * v.z - z * v.y;
    float c_y = z * v.x - x * v.z;
    float c_z = x * v.y - y * v.x;

    Vector3D crossed(c_x, c_y, c_z);

    return crossed;
}

// Take the dot product and return the value
float Vector3D :: dot(Vector3D v){
    float dot_product = x * v.x + y * v.y + z * v.z;
    return dot_product;
}

// Set the vector by value
void Vector3D :: set(float xVal, float yVal, float zVal){
    x = xVal; y = yVal; z = zVal;
}

// Set the vector by vector
void Vector3D :: set(const Vector3D & v){
    x = v.x; y = v.y; z = v.z;
}

// Reverse vector and return
void Vector3D :: reverse(){
    x = -1 * x;
    y = -1 * y;
    z = -1 * z;
}

// Return the difference vector of 2 points
void Vector3D :: setDiff(Point & a, Point & b){
    x = a.x - b.x;
    y = a.y - b.y;
    z = a.z - b.z;
}

// Calculate the magnitude of the vector
float Vector3D :: magnitude(){
    return sqrt(pow(x,2) + pow(y,2) + pow(z,2));
}

// Display the values of the vector
void Vector3D :: display(){
    printf("(%f,%f,%f)\n",x,y,z);
}

// Normalize the vector to have magnitude of one
void Vector3D :: normalize(){

    float mag = magnitude();
    x = x / mag;
    y = y / mag;
    z = z / mag;
}

// Scalarize the vector by a specific factor
void Vector3D :: scalarize(float s){

    x *= s;
    y *= s;
    z *= s;
}

Vector3D Vector3D::normalize(Vector3D input_vec) {
    Vector3D out_vec;
    float mag = sqrt(pow(input_vec.x,2) + pow(input_vec.y,2) + pow(input_vec.z,2));
    out_vec.x = input_vec.x / mag;
    out_vec.y = input_vec.y / mag;
    out_vec.z = input_vec.z / mag;
    return out_vec;
}
