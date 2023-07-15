//========================================================
// Point.cpp
// Donny Chan, Billy Taj

#include "../inc/Point.h"

Point :: Point(){}

Point :: Point(float _x, float _y, float _z){
    x = _x; y = _y, z = _z;
}

Point :: Point(int _x, int _y, int _z){
    x = _x; y = _y, z = _z;
}

Point::Point(std::vector<double> & vec_coordinates){
    x = vec_coordinates[0]; y = vec_coordinates[1], z = vec_coordinates[2];
}

Point :: Point(const Point & p){

    x = p.x;
    y = p.y;
    z = p.z;
}

Point :: ~Point(){}

void Point :: add(Point & p){
    x += p.x;
    y += p.y;
    z += p.z;
}

void Point :: set(const Point & p){
    x = p.x;
    y = p.y;
    z = p.z;
}

float Point :: distFromOrigin(){
    return sqrt(pow(x,2) + pow(y,2) + pow(z,2));
}

Point Point :: mid_point(Point & a, Point & b){
    Point new_point;
    new_point.x = (a.x + b.x)/2;
    new_point.y = (a.y + b.y)/2;
    new_point.z = (a.z + b.z)/2;
    return new_point;
}
