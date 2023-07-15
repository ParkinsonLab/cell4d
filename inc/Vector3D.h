//==================================================
// Vector3D.h
// Donny Chan, Billy Taj

#ifndef VECTOR3D_H
#define VECTOR3D_H

#include "Point.h"
#include <math.h>

class Vector3D{

    public:
        float x,y,z;

        Vector3D();
        Vector3D(float, float, float);
        Vector3D(const Vector3D &);

        Vector3D cross(const Vector3D &);
        float dot(Vector3D);

        void set(float,float,float);
        void set(const Vector3D &);
        void reverse();
        void display();
        void setDiff(Point &, Point &);
        float magnitude();
        void normalize();
        static Vector3D normalize(Vector3D);
        void scalarize(float);
};

#endif
