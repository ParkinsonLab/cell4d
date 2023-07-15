//===========================================================
// Viewer.cpp
// Donny Chan, Billy Taj


#include "../inc/Viewer.h"
#define PI 3.14159265

Viewer :: Viewer(){
    viewAngle = 0; aspect = 1; nearDist = 10; farDist = 100;
}

void Viewer :: set(Vector3D v_eye, Vector3D v_look, Vector3D v_up){
    // Eye Position
    eye.set(v_eye);
    // Construct n
    n.set(eye.x - v_look.x, eye.y - v_look.y, eye.z - v_look.z);
    // u = up X n
    u.set(v_up.cross(n));

    // normalize n, u (unit length)
    n.normalize();
    u.normalize();

    // v = n X u
    v.set(n.cross(u));

//    printf("(%f,%f,%f), (%f,%f,%f), (%f,%f,%f)\n",v_eye.x,v_eye.y,v_eye.z,v_look.x,v_look.y,v_look.z,v_up.x,v_up.y,v_up.z);

    // notify OpenGL
    setModelViewMatrix();

}

void Viewer :: rotate(float angle, int axis){

    int i, j, k;
    i = j = k = 0;
    if(axis == 0)    i = 1;
    if(axis == 1)    j = 1;
    if(axis == 2)    k = 1;

    glRotatef(angle,i,j,k);

    // Roll the viewer through angle degrees
//    float cs = cos(PI/180 * angle);
//    float sn = sin(PI/180 * angle);
//    Vector3D t(i,j,k);
//    n.set(n.x , cs*n.y - sn*n.z, sn*n.y + cs*n.z);
//    u.set(u.x , cs*u.y - sn*u.z, sn*u.y + cs*u.z);
//    v.set(v.x , cs*v.y - sn*v.z, sn*v.y + cs*v.z);
//    printf("(%f,%f,%f), (%f,%f,%f), (%f,%f,%f)\n",n.x,n.y,n.z,u.x,u.y,u.z,v.x,v.y,v.z);
    setModelViewMatrix();
}


void Viewer :: roll(float angle){
    // Roll the viewer through angle degrees
    float cs = cos(PI/180 * angle);
    float sn = sin(PI/180 * angle);
    Vector3D t = u;
    u.set(cs*t.x - sn*v.x, cs*t.y - sn*v.y, cs*t.z - sn*v.z);
    v.set(sn*t.x + cs*v.x, sn*t.y + cs*v.y, sn*t.z + cs*v.z);
    setModelViewMatrix();
}

void Viewer :: pitch(float angle){
    // Pitch the viewer through angle degrees
    float cs = cos(PI/180 * angle);
    float sn = sin(PI/180 * angle);
    Vector3D t = n;
    n.set(cs*t.x - sn*v.x, cs*t.y - sn*v.y, cs*t.z - sn*v.z);
    v.set(sn*t.x + cs*v.x, sn*t.y + cs*v.y, sn*t.z + cs*v.z);
    setModelViewMatrix();
}

void Viewer :: turn(float angle){
    // Pitch the viewer through angle degrees
    float cs = cos(PI/180 * angle);
    float sn = sin(PI/180 * angle);
    Vector3D tu = u;
    Vector3D tn = n;
    u.set(cs*tu.x + sn*tn.x, cs*tu.y + sn*tn.y, cs*tu.z + sn*tn.z);
    n.set(-sn*tu.x + cs*tn.x, -sn*tu.y + cs*tn.y, -sn*tu.z + cs*tn.z);
    setModelViewMatrix();
}

void Viewer :: slide(float delU, float delV, float delN){
//    printf("(%f,%f,%f)\n",delU,delV,delN);
    eye.x += delU * u.x + delV * v.x + delN * n.x;
    eye.y += delU * u.y + delV * v.y + delN * n.y;
    eye.z += delU * u.z + delV * v.z + delN * n.z;
    setModelViewMatrix();
}

void Viewer :: setShape(float vAng, float asp, float nearD, float farD){

    viewAngle = vAng;
    aspect = asp;
    nearDist = nearD;
    farDist = farD;

    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();

    gluPerspective(vAng, asp, nearD, farD);
}

void Viewer :: setModelViewMatrix(){

    // Set the light source property
    GLfloat light_intensity[] = {1.0f, 1.0f, 1.0f, 0.5f};
    GLfloat light_position[] = {n.x, n.y, n.z, 0.5f};

//    glLightfv(GL_LIGHT0, GL_POSITION, light_position);
//    glLightfv(GL_LIGHT0, GL_DIFFUSE, light_intensity);

    // load modelview matrix with existing viewer values;
    float m[16];

    // Eye vector
    Vector3D vec(eye.x, eye.y, eye.z);
//    printf("--> (%f,%f,%f)\n",eye.x, eye.y, eye.z);

    m[0] = u.x;    m[4] = u.y;    m[8] = u.z;    m[12] = vec.dot(u);
    m[1] = v.x;    m[5] = v.y;    m[9] = v.z;    m[13] = vec.dot(v);
    m[2] = n.x;    m[6] = n.y;    m[10] = n.z;    m[14] = vec.dot(n);
    m[3] = 0.0;    m[7] = 0.0;    m[11] = 0.0;    m[15] = 1.0;

    // Make Model View matrix current
    glMatrixMode(GL_MODELVIEW);

    // load OpenGL's modelview matrix
    glLoadMatrixf(m);
}
