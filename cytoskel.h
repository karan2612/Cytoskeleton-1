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

float tmax = 30; //time absolute stop
float dt = 0.005; //time step physics
int ts = 30;      //time step rendering

int _nSYS = 3;  //side length is Twice this #
int _nSpectrin = 7; //number of spectrin Springs between each actin
double _lActin = 0.9; //initial length between Actin
double _Contour = 2.5 * _lActin;

float msd = 0;
double _k = 10; 
double _m = 1;

Ball* Particle = 0;
vector<Ball> v_balls;
vector<Spring> v_springs;

FILE *f1, *f2, *f3, *f4;
FILE *kx, *ky;

void init();
void hexInit();
void meshInit();
void spectrinInit();
void initPID();
void filesInit();
void filesClose();
void initKaran();

void physics();
void timeStep();
void doAnalysis();
void ForceSprings();
void updatePosition(Ball &);
void updateBrownianPosition(Ball &);

void SurfaceForce();
double LJforce(double,double);
double zSurface(double, double);

float calcMSD(float**);
float calcMSDx(float*);

float distBall(Ball*, Ball*);
vector<float> normBall(Ball*, Ball*);
bool isEdge(int, int, int);
int edgeType(int, int, int);

void writeBalls(FILE*);
void writeSprings(FILE*);
void measureEdge(FILE*);
void measureEnergy();
void measureContour(FILE*);
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

void meshInit() {

  int N=_nSYS;
  double L = _lActin;

  double a = -L*(N - 1./2);
  double x0=a,y0=a;
  double z0 = 0;

  double h = sqrt(3)*L/2;

  /* Build Balls */
  Ball *b;
  for (int y=0; y<2*N; y++) {
    for (int x=0; x<2*N; x++) {

      if (y%2 == 0) //even
	{ b = new Ball(x0 + x*L, y0 + y*h, z0); }    
      else //odd
	{ b = new Ball(x0+L/2 + x*L, y0 + y*h, z0); }

      if (isEdge(x,y,N)) {
	b->isEdge = true; 
	b->pid=2;

	b->edgeType = edgeType(x,y,N);
      }

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

  cout << "   mesh complete" << endl;
}

void hexInit() {

  int i = 1;
  Ball a0(0.001, 0 ,i);
  v_balls.push_back(a0);

  double sd = PI/3; //'Sixty Degrees'

  /* Can I automate this in a 2line for loop? */
  Ball a1(1, 0*sd, i);
  Ball a2(1, 1*sd, i);
  Ball a3(1, 2*sd, i);
  Ball a4(1, 3*sd, i);
  Ball a5(1, 4*sd, i);
  Ball a6(1, 5*sd, i);

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

}

void spectrinInit() {

  //makes a copy
  //vector<Spring> vs = v_springs;
  int N = v_springs.size();
  int n = _nSpectrin; //new springs :: 1 should be identity

  Spring *spr;
  Ball *b1, *b2;

  double r1[2],r2[2],dr[2]; //lives in 2D plane
  double x,y,z;

  double A = _lActin;
  double C = _Contour;
  double h = sqrt(C*C + A*A) / 2; //creates dagger strings

  vector<Spring> newSpring_v;
  //for each spring
  int j1,j2;
  for (int j=0; j<N; j++) {
    spr = &v_springs.at(j);

    j1 = spr->n1; 
    j2 = spr->n2; 

    b1 = &v_balls.at(j1);
    b2 = &v_balls.at(j2);

    for (int i=0; i<2; i++) {
      r1[i] = b1->r[i];
      r2[i] = b2->r[i];
      dr[i] = (r2[i] - r1[i]) / n;
    }

    //create n-1 new balls
    int k = v_balls.size(); //used later
    double a = A/2;
    for (int i=1; i<n; i++) {
      x = r1[0] + i*dr[0];
      y = r1[1] + i*dr[1];
      z = -h + h/a * abs( a - (i*A/n) );

      Ball b(x,y,z);
      b.pid = 0;
      v_balls.push_back(b);
    }

    //connect w N new springs (+2 old balls)
    if (n == 1) return;
    newSpring_v.push_back( Spring(j1,k) ); 
    for (int i=1; i<n-1; i++) 
    {
      newSpring_v.push_back( Spring(k, k+1) );
      k++;
    }
    newSpring_v.push_back( Spring(k, j2) );

  }

  //haben Sie memory leak??
  v_springs = newSpring_v;

  cout << "   spectrin spring complete " 
       << n << endl;
}


void initKaran() {

  /* MISSING: x0 init */
  kx = fopen("karanX.txt", "w");
  fprintf(kx, "%i\n" , nTime);
  fprintf(kx, "%i\n\n" , nBalls);

  ky = fopen("karanY.txt", "w");
  fprintf(ky, "%i\n" , nTime);
  fprintf(ky, "%i\n\n" , nBalls);
}

void toyInit() {
  Ball a1(0,0,0);
  Ball a2(2,0,0);
  v_balls.push_back(a1);
  v_balls.push_back(a2);

  Spring s1 = Spring(0,1);
  v_springs.push_back(s1);

}


void fileInit() {

  cout << " initializing file.." << endl;
  f1 = fopen("balls.dat", "w");
  f2 = fopen("springs.dat", "w");
  f3 = fopen("force.txt", "w");
  f4 = fopen("contour.txt", "w");

  /* Writes nTime, nBalls, and nSprings */
  nBalls = v_balls.size();
  nSprings = v_springs.size();

  nTime = (int) tmax/dt;
  nTime = nTime/ts + 1;
  cout << "   predicted time steps: " << nTime << endl;

  fprintf(f1, "%i\n" , nTime);
  fprintf(f1, "%i\n\n" , nBalls);
  fprintf(f2, "%i\n\n" , nSprings);

  /* Defines pid in Render file */
  initPID();

  //msdInit()
  //  if (_msd) f3 = fopen("outMSD.txt", "w");
}

void filesClose() {
  fclose(f1);
  fclose(f2);
  fclose(f3);
  fclose(f4);
  fclose(kx);
  fclose(ky);
  cout << "files written." << endl;
}

void randInit() {

  cout << " initializing rand.." << endl;
  unsigned int saat = (unsigned int)time(0); 
  randi.seed(saat);
}

/* writes PID */
void initPID() {

  int N = nBalls;
  for(int j=0; j<N; j++) { 

    //write pid (1-anchor, 0-spectrin)
    fprintf(f1, " %i", v_balls[j].pid);
  }
  fprintf(f1, "\n\n");

}

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

//returns norm from origin of a to origin of b
vector<float> normBall(Ball *a, Ball *b) {

  float u[3],v[3],n[3];
  for (int i=0;i<3;i++) {
    u[i] = a->r[i];
    v[i] = b->r[i];

    n[i] = v[i] - u[i];
  }

  float m=0; //magnitutde
  for(int j=0; j<3; j++) {
    m += n[j]*n[j];
  }
  m = sqrt(m);
  /* I could make a subroutine,
     but that would pass an array and return a double */

  vector<float> norm = vector<float>(3);
  for(int j=0; j<3; j++) {
    norm[j] = n[j]/m;
  }

  return norm;
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
  //  float a = _lActin/_nSpectrin;
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

/*
void newLJ (Ball *b1, Ball *b2) {
  
  a1 = b1->radius;
  a2 = b2->radius;

  x1[] = b1->r[];
  x2[] = b2->r[];
  r[] = x2[] - x1[];

  R = mag(r);
  n[] = r[] / R;

  a0 = _rSpectrin;

  sigma = 2 * a0;
  //sigma = a1 + a2;

  f = 1.12246;
  m = f * sigma;
  e = _e; //?? 
  // R = R - abs(a2 - a1); ??
  p = m/R;

  F = 12 * e / m * pow(p,7) * ( pow(p,6) - 1 );

*/
    /* I want my force to return ONLY the magnitude
       that I may handle the direction elsewhere
       as dictated by the Case */

double LJforce(double r, double sigma) {

  /* r_min : F = 0, truncate potential
     sigma : U = 0, is sum of radii (a+b) */

  double F, e, m, p;
  e = 0.2; //I choose to absorb 12 into e
  //m = 0.45; 
  m = sigma*1.1225;

  if (r > m) return 0;
  if (r < 0.05) {
    cout << "Warning: LJ Singularity from small r" << endl;
  }

  /* F = 12 e/r [p12 - p6]  
      where p = r_min/r */
  p = m/r;
  p = pow(p,6);
  F =  - e /r * p * (p - 1);

  return F;
}


double zSurface(double x, double y) {
  /* returns Membrane z(x,y) */

  //temp flat
  return 0.1;
}

void initParticle() {

  double z = 2;
  Particle = new Ball(0,0,z);  //need to add radius!

}
// also need shorter name

void ParticleInteraction() {
  //cout << " considering Particle interactions.." << endl;

  if (!Particle) {
    cout << "particle not initialized!" << endl;
    getc(stdin);
  }
 
  double radius = 2.15; //hard code for now
  double r;
  vector<float> nm;
  double F;

  int N = nBalls;
  int count = 0;
  Ball *b;
  for (int j=0; j<N; j++) {
    
    b = & v_balls.at(j);
    r = distBall(b, Particle);

    if ( r < radius ) {
      nm = normBall(b,Particle);
      F = -5;
      // F = LJForce(r);
      for (int i=0; i<2; i++) {
	b->F[i] += F * nm[i];
      }

      count++;
    }
  }   

  //  cout << count << endl;

  //how can I be sure that none are INSIDE the Particle?

}


void measureEdge(FILE* f) {
  int N = nBalls;
  Ball *b;
  double F;
  for (size_t j=0; j<N; j++) {
    b = &v_balls.at(j);
    if (b->edgeType == 1) {
      //cout << b->getForce() << endl;
      F = b->getForce();
      fprintf(f, "%f ", F);
      //hmmm there's probably a c++ way to do this
      //      f << F;
    }
  }

  //  getc(stdin);
  //  cout << endl;
  fprintf(f, "\n");
}

void measureEnergy() {

  int N = nSprings;

  Spring *s;
  float E = 0;
  for (int j=0; j<N; j++) {

    s = &v_springs[j];
    E += s->getEnergy();
  }

  //should I normalize this?
  //  cout << E << endl;
}

void measureContour(FILE *f) {

  /* the springs are created N at a time
     so I can exploit this in reading the segments */

  int N = nSprings;
  int n = _nSpectrin;

  Spring *s;
  float C=0; //contour
  for(int j=0; j<N; j++) {

    s = &v_springs.at(j);
    C += s->getLength();

    if (j == 0) continue;
    if (j%n == 0) {
      //cout << C << endl;
      fprintf(f, "%f ", C);
      C = 0;
    }
  }
  //cout << endl;
  fprintf(f,"\n");
}

/* for misc purposes */
void show() {
  cout << "showing..." << endl;

  int n = v_balls.size();

  Ball* b;
  for (int i=0; i<n; i++) {
    b = &v_balls.at(i);

    float m; //magnitutde
    for(int j=0; j<3; j++) {
      m += b->F[i] * b->F[i];
    }
    m = sqrt(m);    
    cout << m << endl;
  }

  //  getc(stdin);
}

