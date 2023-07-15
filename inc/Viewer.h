//=====================================================
// Viewer.h
// Donny Chan, Billy Taj
// Include OpenGL libraries
#include <GL/gl.h>
#include <GL/glut.h>
#include <X11/X.h>

#include <math.h>
#include <stdio.h>

#include "Point.h"
#include "Vector3D.h"

class Viewer{

    public:
        Viewer();
        void set(Vector3D eye, Vector3D look, Vector3D up);
        void roll(float angle);
        void pitch(float angle);
        void turn(float angle);
        void yaw(float angle);
        void rotate(float angle, int axis);
        void slide(float delU, float delV, float delN);
        void setShape(float vAng, float asp, float nearD, float farD);

    private:
        Vector3D eye;
        Vector3D u,v,n;
        double viewAngle, aspect, nearDist, farDist;
        void setModelViewMatrix();

};
