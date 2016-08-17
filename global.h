#include <iostream>
#include <vector>
#include <cmath>
#include "MersenneTwister.h"

#define PI 3.14159265

using namespace std;

/* Class def */
class Ball {

private:

public:
  Ball (double, double, double); // constructor
  Ball (double, double, int); //polar coordinate constructor
  double r[3],v[3],F[3];

  double m,R; //mass, radius
  int idx, pid; //index in vector, particle id as defined in initPID()
  bool isEdge;  //selects for boundary considerations
  int edgeType; //as defined in edgeType()

  double getForce(); //returns total F on Ball at given time
};

class Spring {

public:
  Spring (int, int); //constructor

  int n1,n2; //nodes: index the balls attached to this brings

  double L, k, m;  //length, spring constant, mass
  double eql; //equilibrium length (F = 0)

  double getLength();
  double getEnergy();
};


/* Global Declaration */
MTRand randi;

bool _msd = false;
int nBalls, nSprings, nTime; 

float tmax = 50; //time absolute stop
float dt = 0.001; //time step physics
int ts = 100;      //time step rendering

int _nSYS = 4;  //side length is Twice this #
int _nSpectrin = 5; //number of spectrin Springs between each actin
double _lActin = 0.6; //initial length between Actin
double _Contour = 2.5 * _lActin;
float _sigma = 0.1; // for now


float msd = 0;
double _k = 20; 
double _m = 1;

FILE *f1, *f2, *f3, *f4, *f5, *f6;
FILE *kx, *ky;

Ball* Particle = 0;
vector<Ball> v_balls;
vector<Spring> v_springs;

void init();
void hexInit();
void meshInit();
void spectrinInit(int);
void initPID();
void filesInit();
void filesClose();
void initKaran();
void moveParticle();

void physics();
void timeStep();
void doAnalysis();
void ForceSprings();
void updatePosition(Ball &);
void updateBrownian(Ball &);

void SurfaceForce();
double LJforce(double,double);
double zSurface(double, double);

float calcMSD(float**);
float calcMSDx(float*);

float distBall(Ball*, Ball*);
vector<float> normBall(Ball*, Ball*);
bool isEdge(int, int, int);
int edgeType(int, int, int);

pair<float,float> doStats(vector<float> &v);
void writeForceZ(FILE *f);

void writeBalls(FILE*);
void writeSprings(FILE*);
void measureEdge(FILE*);
void measureSpringEnergy();
void measureContour(FILE*);
void sampleForceZ();
void sampleForce3D();
void writeKaranXY(FILE*,FILE*);

/* function def */
Ball::Ball (double x, double y, double z) {

  for(int i=0; i<3; i++) {
    r[i]=0;
    v[i]=0;
    F[i]=0;
  }
  r[0] = x;
  r[1] = y;
  r[2] = z;

  idx = v_balls.size();
  pid = 1;
  isEdge = false;

  m = _m;
}

Ball::Ball (double a, double b, int hex) {

  double x,y; 
  if (!hex) {

    Ball(a, b, 0);//does this double init?
    cout << " hex logic Warning!" << endl;
  } else {
    x = a * cos(b);
    y = a * sin(b);

    r[0] = x;
    r[1] = y;
    r[2] = 0;
  }

  idx = v_balls.size();
  pid = 1;
  isEdge = false;

  m = _m;
}

/* returns magnitude of Force */
double Ball::getForce() {

  double f = 0;

  for(int j=0; j<3; j++) {
    f += F[j]*F[j];
  }
  f = sqrt(f);

  return f;
}

Spring::Spring (int b1, int b2) {

  /* connect balls to nodes at construction */
  n1 = b1;
  n2 = b2;

  /* sets equilibrium length 
     as some fraction f of initial length */
  float feq = 1;
  L = getLength();
  eql = feq * L; 
 
  m = _m;
  k = _k;
}

/* determines a spring's lengths from its two nodes */
double Spring::getLength(){

  double L = 0;
  double dr[3];
  
  Ball a = v_balls[n1];
  Ball b = v_balls[n2];

  for (int i=0; i<3; i++) {
    dr[i] = b.r[i] - a.r[i];

    L += dr[i]*dr[i];
  }
  L = sqrt(L);

  return L;
}

double Spring::getEnergy() {

  double E, L;
  L = getLength();
  E = k/2 * (L - eql) * (L - eql);

  return E;
}


bool isEdge(int i, int j, int N) {

  //returns true if corner
  bool foo = false;
  int beg = 0;
  int end = 2*N-1;

  if (i==end || j==beg) foo = true;
  if (i==beg || j==end) foo = true;

  if (i==beg || j==beg) foo = true;
  if (i==end || j==end) foo = true;

  return foo;
}

int edgeType(int i, int j, int N) {

  int beg = 0;
  int end = 2*N-1;

  if (i==end && j==beg) return 0;
  if (i==beg && j==end) return 0;
  if (i==beg && j==beg) return 0;
  if (i==end && j==end) return 0;

  //isn't corner
  if (j==beg || j==end) return 1; //horizontal
  if (i==beg || i==end) return 2; //vertical
  
  cout << "not edge?" << endl;
  getc(stdin);
  return -1;
}

double pow(double x, int n) {

  double ans = 1;
  for(int i=0; i<n; i++) {

    ans = ans*x;
  }
  return ans;
}
