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

  double Fx,Fy;

  int pid, idx;
};

class Spring {

public:

  int n1,n2; //nodes index the balls that the spring is attached to
  double L;
  double eql; //equilibrium length (F = 0)

  Spring (int, int); //constructor
  double Energy(vector<Ball> &);
  double Length();
};


/* General Declaration */
MTRand randi;
double k = 2; 
double m = 1;
double dt = 0.01;
float tmax = 100;
float msd = 0;
bool _msd = false;
FILE *f1, *f2, *f3;


vector<Ball> v_balls;
vector<Spring> v_springs;

void hexInit();
void meshInit();
//vector<Spring> substringInit(int, vector<Ball> &, vector<Spring> &);
vector<Spring> springstringInit(vector<Spring> vs);

float distBall(Ball, Ball);
void updatePosition(Ball &);
void updateBrownianPosition(Ball &b);
void ForceSprings();
void physics();

void writePositions(FILE*);
void writeSprings(FILE*);

float calcMSD(float*);
void initMSD();
void MSD();

/* function def */
Ball::Ball (double x0, double y0) {

  x = x0;
  y = y0;
  z = 0;

  vx = 0;
  vy = 0;  

  Fx = 0;
  Fy = 0;

  idx = v_balls.size();
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

    Fx=0;
    Fy=0;

    idx = v_balls.size();
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

  /* connect balls to nodes at construction */
  n1 = b1;
  n2 = b2;

  /* sets equilibrium length 
     as some fraction f of initial length */
  float f = 0.95;
  L = Length();

  eql = f * L;
  
}

double Spring::Length(){

  /* determines a spring's lengths from its two nodes */

  double dx, dy, L;
  
  Ball a = v_balls[n1];
  Ball b = v_balls[n2];

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


void meshInit() {

  int N=3;
  double L=1.2;  

  double a = -L*(N - 1./2);
  double x0=a,y0=a;

  double h = sqrt(3)*L/2;

  /* Build Balls */
  for (int y=0; y<2*N; y++) {

    if (y%2 == 0) {
      //even
      for (int x=0; x<2*N; x++) {
	Ball b(x0 + x*L, y0 + y*h); 
	v_balls.push_back(b);
      }

    } else {
      //odd
      for (int x=0; x<2*N; x++) {
	Ball b(x0+L/2 + x*L, y0 + y*h);
	v_balls.push_back(b);
      }
    }
  }

  /* Build Springs */
  int j;
  for (int y=0; y<2*N; y++) {
    for (int x=0; x<2*N; x++) {

      j = 2*N*y + x;

      if (x != 2*N-1) {
	Spring s(j,j+1);
	v_springs.push_back(s);
      }

      if (y != 2*N-1) {
	Spring s1(j, j+2*N);
	v_springs.push_back(s1);

	if (y%2 == 0) { 
	  /* y even */
	  if (x != 0) {
	    Spring s0(j, j+2*N - 1);
	    v_springs.push_back(s0);
	  }
	} else { 
	  /* y odd */
	  if (x != 2*N-1) {
	    Spring s0(j, j+2*N + 1);
	    v_springs.push_back(s0);
	  }
	}
      }// y != 2*N-1

    }
  }

}
void hexInit() {

  Ball a0(0.001,0,0);
  v_balls.push_back(a0);

  double sd = PI/3; //'Sixty Degrees'

  /* Can I automate this in a 2line for loop? */
  Ball a1(1, 0*sd, 1);
  Ball a2(1, 1*sd, 1);
  Ball a3(1, 2*sd, 1);
  Ball a4(1, 3*sd, 1);
  Ball a5(1, 4*sd, 1);
  Ball a6(1, 5*sd, 1);

  v_balls.push_back(a1);
  v_balls.push_back(a2);
  v_balls.push_back(a3);
  v_balls.push_back(a4);
  v_balls.push_back(a5);
  v_balls.push_back(a6);

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

  v_springs.push_back(s01);
  v_springs.push_back(s02);
  v_springs.push_back(s03);
  v_springs.push_back(s04);
  v_springs.push_back(s05);
  v_springs.push_back(s06);
  v_springs.push_back(s12);
  v_springs.push_back(s23);
  v_springs.push_back(s34);
  v_springs.push_back(s45);
  v_springs.push_back(s56);
  v_springs.push_back(s61);

  // now redundant
  for (int i=0; i<12; i++) {
    v_springs[i].eql = 0.95;
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

vector<Spring> springstringInit(vector<Spring> vs) {
  //makes a copy
  int N = vs.size();
  int n = 3; //subdivisions
  Spring *spr;
  Ball *b1, *b2;
  double dx, dy;
  
  vector<Spring> newSpring_v;
  //for each spring
  for (int j=0; j<N; j++) {
    spr = &vs.at(j);

    int j1 = spr->n1; 
    int j2 = spr->n2; 

    b1 = &v_balls.at(j1);
    b2 = &v_balls.at(j2);

    dx = (b2->x - b1->x) / n;
    dy = (b2->y - b1->y) / n;

    //create N-1 new balls
    for (int i=1; i<n-1; i++) {
      Ball b(b1->x + i*dx, b1->y + i*dy);
      v_balls.push_back(b);
    }

    //connect w N new springs (+2 old balls)
    int k = v_balls.size() - 1;

    newSpring_v.push_back( Spring(j2, k) ); 
    for (int i=1; i<n; i++) {
      newSpring_v.push_back( Spring(k-i, k-i-1) );
    }
    newSpring_v.push_back( Spring(k-n, j1) );

  }

  return newSpring_v;
}







// For deletion 
void initString() {

  /* Sets Initial Conditions
     we start w 1D spring for 1D forces*/

  Ball a1(0,0);
  Ball b0(1,0);
  Ball a2(3,0); 

  v_balls.push_back(a1);
  v_balls.push_back(b0);
  v_balls.push_back(a2);

  /* need to build springs */

  Spring s1 = Spring(0,1);
  Spring s2 = Spring(1,2);

  s1.eql = (a2.x - a1.x) / 2;
  s2.eql = (a2.x - a1.x) / 2;

  v_springs.push_back(s1);
  v_springs.push_back(s2);

}


//also for deletion
void show() {
  cout << "showing..." << endl;

  int n = v_balls.size();
  cout << "Total Balls: " << n << endl;

  Ball* foo;
  for (int i=0; i<n; i++) {
    foo = &v_balls.at(i);

    cout << i << " "
	 << foo->idx << endl;
  }
  getc(stdin);
}

