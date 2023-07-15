//=================================================
// Point.h
// Donny Chan, Billy Taj

#ifndef POINT_H
#define POINT_H

#include <math.h>
#include <vector>

class Point{

    public:

        Point();
        Point(const Point &);
        Point(float, float, float);
        Point(int, int, int);
        Point(std::vector<double> &);

        ~Point();

        float x;
        float y;
        float z;

        void add(Point &);
        void set(const Point &);
        float distFromOrigin();
        static Point mid_point(Point &, Point &);
};

#endif
