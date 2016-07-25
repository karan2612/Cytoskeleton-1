#include <iostream>
#include <vector>
#include <math.h>

using namespace std;

/* Class def */
class Ball {

private:

public:
  Ball (double, double); // ctor

  double x,y,z;
  double vx,vy;

  double T,U;

  double F;
};

class Spring {

public:

  int n1,n2; //nodes index the balls that the spring is attached to

  Spring (int, int); //constructor
  double xeq; //equilibrium length (F = 0)
  double Energy(vector<Ball> &v);
  double Length();
};


/* adjustable parameters, globally defined */
double k = 2; 
double L = 1; 
double m = 1;
double dt = 0.1;
float tmax=100;

vector<Ball> ball_v;
vector<Spring> spring_v;



/* function def */
Ball::Ball (double x0, double y0) {

  x = x0;
  y = y0;
  z = 0;

  vx = 0;
  vy = 0;  

  F = 0;
}

Spring::Spring (int b1, int b2) {

  //xeq = L;
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

  if (0) {
    cout << dx << ", "
	 << dy << ", "
	 << L  << endl;
  }

  return L;

}

double Spring::Energy(vector<Ball> &v) {

  /* U = 1/2 k x*x */
  double dx = v[n1].x - v[n2].x;
  double dy = v[n1].y - v[n2].y;

  double mag = sqrt(dx*dx + dy*dy);

  double E = k/2 * (mag - xeq) * (mag - xeq);

  if (0) {
  cout << "Energy: dx dy mag xeq E\n "
       << dx << ", "  
       << dy << ", " 
       << mag << ", " 
       << xeq << ", " 
       << E << endl;
  }
  return E;
}
//maybe it would be better to store this as a force


void init() {

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

  s1.xeq = (a2.x - a1.x) / 2;
  s2.xeq = (a2.x - a1.x) / 2;

  spring_v.push_back(s1);
  spring_v.push_back(s2);

}


/* this test script only moves one spring, 
   so I have something to practice rendering with*/
void updatePosition(Ball &b) {

  //compute a
  double ax = b.F / m;
  
  //update v
  b.vx += ax * dt;

  //update x
  b.x += b.vx * dt;

  if (0) {
    cout << ax << ", "
	 << b.vx << ", "
	 << b.x  
	 << endl;
  }
}


void physics() {

  /* Everything here should be called once per time step*/

  for (size_t j=0; j<ball_v.size(); j++) {
    ball_v[j].F = 0;
  }

  /* loop through springs to compute force on all balls */
  Spring* spr;
  double L, X, F;

  for (size_t j=0; j<spring_v.size(); j++) {

    spr = &spring_v[j];

    L = spr->Length();
    X = spr->xeq;
    /* 
       either both are attracted
       or both repelled
       No matter what there will be opposite sign
    */

    F = -k * (L - X);

    if (0) {
    cout << L << ", "
	 << X << ", "
	 << F 
	 << endl;
    }

    ball_v[spr->n1].F -= F/2;
    ball_v[spr->n2].F += F/2;
    /* this /2 does not belong
       if on edge...
     */
  }

  /* add brownian motion to force bucket? */

  //loop over all moving particles, now there is only 1.
  updatePosition(ball_v[1]); 

}

void writePositions(FILE* f) {

  int n = ball_v.size();
  for(int j=0; j<n; j++) {
    fprintf(f, "%f" , ball_v[j].x);
    fprintf(f, " %f", ball_v[j].y);
    fprintf(f, " %f", ball_v[j].z);
    fprintf(f, "\n");
  }

  fprintf(f, "\n"); //new timestep
}

void draw() {

  /* 
     rather than write and read from .txt
     what if I ported directly into an array?
  */


  int n = ball_v.size();

  double x[3*n];

  for(int i=0; i<n; i++) {

    x[3*i + 0] = ball_v[i].x;
    x[3*i + 1] = ball_v[i].y;
    x[3*i + 2] = ball_v[i].z;

  }

  /*
  for (int i = 0; i < 3*n; i++) 
    cout << x[i] << ", ";
  cout << endl;
  */


}


int main() {

  init();

  FILE* fout;
  fout = fopen("newPositions.txt", "w");
  int T = (int)tmax/dt;
  int n = ball_v.size();
  fprintf(fout, "%i\n" , T);
  fprintf(fout, "%i\n\n" , n);

  float t=0;
  while (t<tmax) {
    physics();
    writePositions(fout);
    //draw();

    t += dt;
  }

  fclose(fout);
  return 0;
}
