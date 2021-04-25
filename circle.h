 
#ifndef CIRCLE_H
#define CIRCLE_H
//
//						 circle.h
//
/************************************************************************
            DECLARATION OF THE CLASS CIRCLE
************************************************************************/
// Class for Circle
// A circle has 7 fields:
//     a, b, r (of type reals), the circle parameters
//     s (of type reals), the estimate of sigma (standard deviation)
//     g (of type reals), the norm of the gradient of the objective function
//     i and j (of type int), the iteration counters (outer and inner, respectively)
#include <iostream>
#include <cmath>
#include <limits>
#include <iomanip>
#include <cstdlib>
using namespace std;


class Circle
{
public:

    // The fields of a Circle
    double a, b, r, s, g, Gx, Gy;
    int i, j;

    // constructors
    Circle();
    Circle(double aa, double bb, double rr);

    // routines
    void print(void);

    // no destructor we didn't allocate memory by hand.
};

class CircleData
{
public:

    int n;
    double *X;		//space is allocated in the constructors
    double *Y;		//space is allocated in the constructors
    double meanX, meanY;

    // constructors
    CircleData();
    CircleData(int N);
    CircleData(int N, double X[], double Y[]);

    // routines
    void means(void);
    void center(void);
    void scale(void);
    void print(void);

    // destructors
    ~CircleData();
};


//   Note: long double is an 80-bit format (more accurate, but more memory demanding and slower)

typedef long long integers;

//   next define some frequently used constants:

const double One=1.0,Two=2.0,Three=3.0,Four=4.0,Five=5.0,Six=6.0,Ten=10.0;
//const reals One=1.0L,Two=2.0L,Three=3.0L,Four=4.0L,Five=5.0L,Six=6.0L,Ten=10.0L;
const double Pi=3.141592653589793238462643383L;
const double REAL_MAX=numeric_limits<double>::max();
const double REAL_MIN=numeric_limits<double>::min();
const double REAL_EPSILON=numeric_limits<double>::epsilon();

//   next define some frequently used functions:

template<typename T>
inline T SQR(T t) { return t*t; };

Circle CircleFitByHyper (CircleData& data);
Circle CircleFitByKasa (CircleData& data);
Circle CircleFitByPratt (CircleData& data);
Circle CircleFitByTaubin (CircleData& data);
int CircleFitByLevenbergMarquardtFull (CircleData& data, Circle& circleIni, double LambdaIni, Circle& circle);

#endif // CIRCLE_H


