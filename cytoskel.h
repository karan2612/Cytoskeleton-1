#include <iostream>
#include <vector>
//#include <math.h>
#include <cmath>
#include "MersenneTwister.h"

#define PI 3.14159265

using namespace std;

/* Class def */
class Ball {

private:

public:
  Ball (double, double); // ctor
  Ball (double, double, int); //polar
  double r[3],v[3],F[3];

  double T,U;
  int pid, idx;
  bool isCorner;
  bool isEdge;
};

class Spring {

public:

  int n1,n2; //nodes index the balls that the spring is attached to
  double L;
  double eql; //equilibrium length (F = 0)

  Spring (int, int); //constructor
  //  double Energy(vector<Ball> &);
  double Length();
};


/* Global Declaration */
MTRand randi;

bool _msd = false;
int nBalls, nSprings, nTime; 

float tmax = 5; //time absolute stop
float dt = 0.005; //time step physics
int ts = 8;      //time step rendering

int _nSYS = 4;
float _LENGTH = 0.9;
float _nSpectrin = 4;

float msd = 0;
double k = 2; 
double m = 1;

FILE *f1, *f2, *f3;
FILE *kx, *ky;

vector<Ball> v_balls;
vector<Spring> v_springs;

void init();
void hexInit();
void meshInit();
vector<Spring> subInit(vector<Spring> vs);
void filesInit();
void filesClose();
void initKaran();

void physics();
void timeStep();
void ForceSprings();
void BoundarySprings(); //?
void updatePosition(Ball &);
void updateBrownianPosition(Ball &b);

float calcMSD(float**);
float calcMSDx(float*);

float distBall(Ball*, Ball*);
bool isCorner(int, int, int);
bool isEdge(int, int, int);

void writeBalls(FILE*);
void writeSprings(FILE*);
void writeKaranXY(FILE*,FILE*);

/* function def */
Ball::Ball (double x0, double y0) {

  for(int i=0; i<3; i++) {
    r[i]=0;
    v[i]=0;
    F[i]=0;
  }
  r[0] = x0;
  r[1] = y0;

  idx = v_balls.size();
  pid = 1;
  isCorner = false;
  isEdge = false;
}

Ball::Ball (double r, double t, int hex) {

  if (!hex) {
    Ball(r, t);

  } else {

    float x,y; //does this double init?
    x = r * cos(t);
    y = r * sin(t);
    Ball(x, y);
  }

}


Spring::Spring (int b1, int b2) {

  /* connect balls to nodes at construction */
  n1 = b1;
  n2 = b2;

  /* sets equilibrium length 
     as some fraction f of initial length */
  float f = 0.98;
  L = Length();

  eql = f * L;  
}

/* determines a spring's lengths from its two nodes */
double Spring::Length(){

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


void meshInit() {

  int N=_nSYS;
  double L = _LENGTH;

  double a = -L*(N - 1./2);
  double x0=a,y0=a;

  double h = sqrt(3)*L/2;

  /* Build Balls */
  Ball *b;
  for (int y=0; y<2*N; y++) {
    for (int x=0; x<2*N; x++) {

      if (y%2 == 0) //even
	{ b = new Ball(x0 + x*L, y0 + y*h); }    
      else //odd
	{ b = new Ball(x0+L/2 + x*L, y0 + y*h); }

      if (isCorner(x,y,N)) 
	{b->isCorner = true; }
      if (isEdge(x,y,N)) 
	{b->isEdge = true; b->pid=2;}

      v_balls.push_back(*b);
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
float distBall(Ball *b1, Ball *b2) {

  float L = 0;
  float r1[3],r2[3],dr[3];
  for (int i=0; i<3; i++) {

    r1[i] = b1->r[i];
    r2[i] = b2->r[i];

    dr[i] = r2[i] - r1[i];

    L += dr[i]*dr[i];
  }
  L = sqrt(L);
  
  return L;
}

vector<Spring> subInit(vector<Spring> vs) {

  //makes a copy
  int N = vs.size();
  int n = _nSpectrin; //new springs :: 1 should be identity
  Spring *spr;
  Ball *b1, *b2;
  //  double dx, dy;
  double r1[3],r2[3],dr[3];
  double x,y;

  vector<Spring> newSpring_v;
  //for each spring
  for (int j=0; j<N; j++) {
    spr = &vs.at(j);

    int k = v_balls.size(); 
    int j1 = spr->n1; 
    int j2 = spr->n2; 

    b1 = &v_balls.at(j1);
    b2 = &v_balls.at(j2);

    for (int i=0; i<3; i++) {
      r1[i] = b1->r[i];
      r2[i] = b2->r[i];
      dr[i] = (r2[i] - r1[i]) / n;
    }

    //create N-1 new balls
    for (int i=1; i<n; i++) {
      x = r1[0] + i*dr[0];
      y = r1[1] + i*dr[1];
      Ball b(x,y);
      b.pid = 0;
      v_balls.push_back(b);
    }

    //connect w N new springs (+2 old balls)
    if (n == 1) return vs;

    newSpring_v.push_back( Spring(j1,k) ); 
    for (int i=1; i<n-1; i++) {
      newSpring_v.push_back( Spring(k, k+1) );
      k++;
    }
    newSpring_v.push_back( Spring(k, j2) );

  }

  return newSpring_v;
}


void initKaran() {

  /* temp */
  kx = fopen("karanX.txt", "w");
  fprintf(kx, "%i\n" , nTime);
  fprintf(kx, "%i\n\n" , nBalls);

  ky = fopen("karanY.txt", "w");
  fprintf(ky, "%i\n" , nTime);
  fprintf(ky, "%i\n\n" , nBalls);
}

void fileInit() {

  cout << " initializing file.." << endl;
  f1 = fopen("balls.dat", "w");
  f2 = fopen("springs.dat", "w");

  nBalls = v_balls.size();
  nSprings = v_springs.size();

  nTime = (int) tmax/dt;
  nTime = nTime/ts + 1;
  cout << "predicted time steps: " << nTime << endl;

  fprintf(f1, "%i\n" , nTime);
  fprintf(f1, "%i\n\n" , nBalls);
  fprintf(f2, "%i\n\n" , nSprings);

  //msdInit()
  if (_msd) f3 = fopen("outMSD.txt", "w");
}

void filesClose() {
  fclose(f1);
  fclose(f2);
  fclose(f3);
  fclose(kx);
  fclose(ky);
  cout << "files written." << endl;
}

void randInit() {

  cout << " initializing rand.." << endl;
  unsigned int saat = (unsigned int)time(0); 
  randi.seed(saat);
}

/* writes Initial Positions and PID */
void initBalls() {

  int N = nBalls;
  float r0[N][3];
  for(int j=0; j<N; j++) { 

    //write pid (1-anchor, 0-spectrin)
    fprintf(f1, " %i", v_balls[j].pid);

    //write initial positions
    for(int i=0; i<3; i++) {
      r0[j][i] = v_balls.at(j).r[i];
    }
    /* does't go anywhere yet! */

  }
  fprintf(f1, "\n\n");

}


void writeKaranXY(FILE* fx, FILE* fy) {

  size_t N = v_balls.size();
  float x=0;
  float y=0;
  for (int j=0; j<N; j++) {
    x = v_balls[j].r[0];
    fprintf(fx,"%f ",x);

    y = v_balls[j].r[1];
    fprintf(fy,"%f ",y);
  }
  fprintf(fx,"\n");
  fprintf(fy,"\n");
}


/* new development: no spring, just brown*/
float calcMSD(float **r0) {

  size_t N = v_balls.size(); 

  float rt[N][3];
  float sum=0;  

  for (size_t j=0; j<N; j++) {
    /* Zeros Forces */
    for (int i=0; i<3; i++) {
      v_balls[j].F[i] = 0;
    }

    /* update physics (no springs) */
    updateBrownianPosition(v_balls[j]); 

    /* compute MSD */
    for (int i=0; i<3; i++) {
      rt[j][i] = v_balls.at(j).r[i];
      sum += (rt[j][i] - r0[j][i])
   	    *(rt[j][i] - r0[j][i]);
    }
  }

  sum = sum/N; //norm
  return sum;
}

float calcMSDx(float *x0) {

  size_t N = v_balls.size(); 

  float xt[N];
  float sum=0;  

  for (size_t j=0; j<N; j++) {
    /* Zeros Forces */
    for (int i=0; i<3; i++) {
      v_balls[j].F[i] = 0;
    }

    /* update physics (no springs) */
    updateBrownianPosition(v_balls[j]); 

    /* compute MSD */
    xt[j] = v_balls.at(j).r[0];
    sum += (xt[j] - x0[j])
          *(xt[j] - x0[j]);
  }
  //cout << sum << " / " << N;
  sum = sum / N; //norm
  //cout << " = " << sum << endl;
  return sum;
}

/*Vestigial, this will be axed soon
  but I think the Algorithm may be reincarnated in L-J Force
*/
void findEdges() {

  Ball *b1, *b2;
  float L = 0;
  //  float a = _LENGTH/_nSpectrin;
  float a = 1.01;
  cout << a << endl;
  
  int count = 0;
  int N = v_balls.size();
  for (int i=0; i<N; i++) {
    for (int j=0; j<N; j++) {

      if (i==j) continue;
      b1 = &v_balls[i];
      b2 = &v_balls[j];
      
      L = distBall(b1, b2);

      L = abs(L - a*_nSpectrin);
      if (L < 0.2) {
	count++;

	if (0) {
	cout << i << " "
	     << j << " "
	     << L << endl;}

      }
    }
  }
  cout << "\n" << count << endl;
  getc(stdin);
}

bool isCorner(int i, int j, int N) {

  //returns true if corner
  bool foo = false;
  int beg = 0;
  int end = 2*N-1;

  if (i==end && j==beg) foo = true;
  if (i==beg && j==end) foo = true;

  if (i==beg && j==beg) foo = true;
  if (i==end && j==end) foo = true;

  return foo;
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

/* for misc purposes */
void show() {
  cout << "showing..." << endl;

  int n = v_balls.size();
  cout << "Total Balls: " << n << endl;
  cout << "(i, idx, pid)" << endl;

  Ball* foo;
  for (int i=0; i<n; i++) {
    foo = &v_balls.at(i);

    cout << i << " "
	 << foo->idx << " "
	 << foo->pid << endl;
  }
  getc(stdin);
}

