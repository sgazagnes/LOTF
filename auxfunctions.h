
#pragma once
#ifndef AUXFUNCTIONS_H
#define AUXFUNCTIONS_H

// Local headers


#define MIN(a,b)  ((a<=b) ? (a) : (b))
#define MAX(a,b)  ((a>=b) ? (a) : (b))

#define pdd pair<double, double> ;


typedef enum {
  DOWN = 1,
  SAME = 2,
  UP = 4
} Direction;

typedef enum {
  ONGOING = 1,
  TOMERGE = 2,
  FINISHED = 3
} Status;

struct point3D{
  // public:
point3D()
: m_x(0.000),
    m_y(0.000),
    m_z(0.000)
  {};
  
explicit point3D(float x, float y, float z)
  : m_x(x),
    m_y(y),
    m_z(z)
  {};
  
point3D(point3D const &ot)
  : m_x(ot.m_x),
    m_y(ot.m_y),
    m_z(ot.m_z)
  {};

  virtual ~point3D(){};

  point3D& operator=(point3D const &ot)
  {
    if(this != &ot) {
      this->m_x = ot.m_x;
      this->m_y = ot.m_y;
      this->m_z = ot.m_z;
    }
    return (*this);
  };

  // member parameters.
  float m_x;
  float m_y;
  float m_z;

private:
  bool operator==(point3D const &ot) const;
  bool operator>(point3D const &ot) const;
  bool operator<(point3D const &ot) const;
};

struct CurvatureParameters{
  // public:
  // Constructors
CurvatureParameters()
: m_a(0.0), m_b(0.0), m_ra(0.0), m_r(0.0), m_E(0.0)
  {};
  
  explicit CurvatureParameters(double a, double b, double ra, double r, double e)
    : m_a(a), m_b(b), m_ra(ra), m_r(r), m_E(e)
  {};
  // D-tor
  virtual ~CurvatureParameters(){};
  
  // Copy constructor
CurvatureParameters(CurvatureParameters const &ot)
: m_a(ot.m_a), m_b(ot.m_b), m_ra(ot.m_ra), m_r(ot.m_r), m_E(ot.m_E)
  {};
  // Assignment operator
  CurvatureParameters& operator=(CurvatureParameters const &ot)
  {
    if(this != &ot) {
      this->m_a = ot.m_a;
      this->m_b = ot.m_b;
      this->m_ra = ot.m_ra;
      this->m_r = ot.m_r;
      this->m_E = ot.m_E;
    }
    return (*this);
  };
  // data Members
  double m_a;// x coord of circle center
  double m_b;// y_coord of circle center
  double m_ra;//(radius)
  double m_r;// (1/Radius)
  double m_E;// Energy value(minimized)
  // protected:
private:
  // Operators to avoid mistakes
  bool operator==(CurvatureParameters const &ot) const;
  bool operator>(CurvatureParameters  const &ot) const;
  bool operator<(CurvatureParameters  const &ot) const;
};
//------------------- END class CurvatureParameters -------------


#include "CoordGrid.h"
#include "gridNode.h"
#include "pathCandidate.h"
#include "TNtuple.h"

bool sortbysec(const pair<int,unsigned short> &a,   const pair<int,unsigned short> &b);

double Cartesian_To_Polar(float const x, float const y,   std::pair<float,float>& polarOut, bool useSign = true);

TNtuple* GridToNtuple(std::vector < GridNode > const &Nodes, std::string const &name = "TestTuple");

double *polyFit(std::vector<double>  x, std::vector<double>  y, int n);

int returnDirection(double prev, double cur);

double returnAngle(double x1, double x2, double x3, double y1, double y2, double y3);

double returnCurvature(double x1, double x2, double x3, double y1, double y2, double y3);

double distanceBetweenTube(GridNode &tubeA, GridNode &tubeB);

double IntersectionPointTubes(GridNode &tubeA, GridNode &tubeB, GridNode &out, int xymid = 0);

bool LineLineIntersect( GridNode &tubeA, GridNode &tubeB, GridNode &tubeC, float &ixOut, float &iyOut, float &izOut); //Output
bool PointsLineIntersectLive( GridNode &tubeC, float x1, float x2, float y1, float y2); //Output 

bool PointsLineIntersectFinal( GridNode &tubeC, float x1, float x2, float y1, float y2); //Output 

void TubeIntersectionPointCoord(GridNode &tubeA, GridNode &tubeB);

void fit_circle(std::vector<point3D> const &pnts, CurvatureParameters &curvature);


#endif
