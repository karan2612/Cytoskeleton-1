#include <iostream>
#include <vector>
#include <math.h>
#include "MersenneTwister.h"

#define PI 3.14159265

using namespace std;

/* Class def */
class Ball {

private:

public:
  Ball (double, double); // ctor
  Ball (double, double, int); //polar
  double dist2ball(Ball &);

  double x,y,z;
  double vx,vy;

  double T,U;

  double F;
  double Fx,Fy;
};

class Spring {

public:

  int n1,n2; //nodes index the balls that the spring is attached to
  double eql; //equilibrium length (F = 0)

  Spring (int, int); //constructor
  double Energy(vector<Ball> &);
  double Length();
};


/* General Declaration */
MTRand randi;
double k = 2; 
double L = 1; 
double m = 1;
double dt = 0.05;
float tmax = 100;
int tw = 10;

vector<Ball> ball_v;
vector<Spring> spring_v;
/*
Ball::Ball(double, double);
Ball::Ball(double, double, int);
double Ball::dist2ball(Ball &);
Spring::Spring(int, int);
double Spring::Length();
double Spring::Energy(vector<Ball> &);
*/
void hexInit();
vector<Spring> substringInit(int, vector<Ball> &, vector<Spring> &);
float distBall(Ball, Ball);
void updatePosition(Ball &);
void ForceSprings();
void physics();

void writePositions(FILE*);
void writeSprings(FILE*);


/* function def */
Ball::Ball (double x0, double y0) {

  x = x0;
  y = y0;
  z = 0;

  vx = 0;
  vy = 0;  

  F = 0; //to be retired
  Fx = 0;
  Fy = 0;
}

Ball::Ball (double r, double t, int hex) {

  if (!hex) {
    //cout << " but I should be ctor 1.. " << endl;
    Ball(r, t);

  } else {

    x = r * cos(t);
    y = r * sin(t);

    z = 0;
    vx = 0;
    vy = 0;  
    F = 0;
    Fx=0;
    Fy=0;
  }

}


double Ball::dist2ball(Ball &b) {

  double x,y,z;
  double dx,dy,dz;
  
  dx = b.x - x;
  dy = b.y - y;
  dz = b.z - z;

  double L;
  L = sqrt( dx*dx + dy*dy + dz*dz );

  return L;
}
//b1->dist2ball(b2);


Spring::Spring (int b1, int b2) {

  //eql = L;
  /* connect balls to nodes at construction */
  n1 = b1;
  n2 = b2;

}

double Spring::Length(){

  /* determines a spring's lengths from its two nodes */

  double dx, dy, L;
  
  Ball a = ball_v[n1];
  Ball b = ball_v[n2];

  dx = b.x - a.x;
  dy = b.y - a.y;

  L = sqrt( dx*dx + dy*dy );

  return L;
}

//unsued 
double Spring::Energy(vector<Ball> &v) {

  /* U = 1/2 k x*x */
  double dx = v[n1].x - v[n2].x;
  double dy = v[n1].y - v[n2].y;

  double mag = sqrt(dx*dx + dy*dy);

  double E = k/2 * (mag - eql) * (mag - eql);

  return E;
}



void hexInit() {

  Ball a0(0.001,0,0);
  ball_v.push_back(a0);

  double sd = PI/3; //'Sixty Degrees'

  /* Can I automate this in a 2line for loop? */
  Ball a1(1, 0*sd, 1);
  Ball a2(1, 1*sd, 1);
  Ball a3(1, 2*sd, 1);
  Ball a4(1, 3*sd, 1);
  Ball a5(1, 4*sd, 1);
  Ball a6(1, 5*sd, 1);

  ball_v.push_back(a1);
  ball_v.push_back(a2);
  ball_v.push_back(a3);
  ball_v.push_back(a4);
  ball_v.push_back(a5);
  ball_v.push_back(a6);

  /* now for springs */
  Spring s01 = Spring(0,1);
  Spring s02 = Spring(0,2);
  Spring s03 = Spring(0,3);
  Spring s04 = Spring(0,4);
  Spring s05 = Spring(0,5);
  Spring s06 = Spring(0,6);

  Spring s12 = Spring(1,2);
  Spring s23 = Spring(2,3);
  Spring s34 = Spring(3,4);
  Spring s45 = Spring(4,5);
  Spring s56 = Spring(5,6);
  Spring s61 = Spring(6,1);

  spring_v.push_back(s01);
  spring_v.push_back(s02);
  spring_v.push_back(s03);
  spring_v.push_back(s04);
  spring_v.push_back(s05);
  spring_v.push_back(s06);
  spring_v.push_back(s12);
  spring_v.push_back(s23);
  spring_v.push_back(s34);
  spring_v.push_back(s45);
  spring_v.push_back(s56);
  spring_v.push_back(s61);

  for (int i=0; i<12; i++) {
    spring_v[i].eql = 0.95;
  } 

}


/* what is the cost of not passing by ref? */
float distBall(Ball a, Ball b) {

  float x,y,z;
  x = b.x - a.x;
  y = b.y - a.y;
  z = b.z - a.z;

  float L;
  L = sqrt( x*x + y*y + z*z );

  return L;
}


vector<Spring> substringInit(int N, vector<Ball> &vb, vector<Spring> &vs) {

  /* subdivides each spring into N smaller springs*/
  vector<Spring> foo;
  return foo;
}






// For deletion 
void initString() {

  /* Sets Initial Conditions
     we start w 1D spring for 1D forces*/

  Ball a1(0,0);
  Ball b0(1,0);
  Ball a2(3,0); 

  ball_v.push_back(a1);
  ball_v.push_back(b0);
  ball_v.push_back(a2);

  /* need to build springs */

  Spring s1 = Spring(0,1);
  Spring s2 = Spring(1,2);

  s1.eql = (a2.x - a1.x) / 2;
  s2.eql = (a2.x - a1.x) / 2;

  spring_v.push_back(s1);
  spring_v.push_back(s2);

}


//also for deletion
void show() {
  cout << "showing..." << endl;

  int n = ball_v.size();
  cout << n << endl;

  Ball* foo;
  for (int i=0; i<n; i++) {
    foo = &ball_v.at(i);

    cout << foo->x << " "
	 << foo->y << endl;
  }
  getc(stdin);
}

