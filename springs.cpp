#include <iostream>
#include <vector>
#include <math.h>

using namespace std;

/* Class def */
class Ball {

private:

public:
  Ball (double, double); // ctor

  double x,y;
  double vx,vy;

  double T,U;

};

class Spring {

  int n1,n2; //nodes index the balls that the spring is attached to
  double xeq; //equilibrium length (F = 0)

public:
  Spring (int, int); //constructor
  double Energy(vector<Ball> &v);
};


/* adjustable parameters, globally defined */
double k = 2; 
double L = 1; 
double m = 1;

vector<Ball> balls;

/* function def */
Ball::Ball (double x0, double y0) {

  x = x0;
  y = y0;

  vx = 0;
  vy = 0;  

}

Spring::Spring (int b1, int b2) {

  xeq = L;
  /* connect balls to nodes at construction */
  n1 = b1;
  n2 = b2;

}

double Spring::Energy(vector<Ball> &v) {

  /* U = 1/2 k x*x */
  double dx = v[n1].x - v[n2].x;
  double dy = v[n1].y - v[n2].y;

  double mag = sqrt(dx*dx + dy*dy);

  double E = k/2 * (mag - xeq) * (mag - xeq);

  cout << "Energy: dx dy mag xeq E\n "
       << dx << ", "  
       << dy << ", " 
       << mag << ", " 
       << xeq << ", " 
       << E << endl;

  return E;
}
//maybe it would be better to store this as a force


void init() {

  /* Sets Initial Conditions
     we start w 1D spring for 1D forces*/

  Ball foo(0,0);
  Ball bar(1.5,0); 

  balls.push_back(foo);
  balls.push_back(bar);

}


/* this test script only moves one spring, 
   so I have something to practice rendering with*/
void physicsTest(Ball &b, float dt) {

  double xo = balls[0].x;

  //compute a: foo on bar

  //cout << b.x << endl;
  double dx = (b.x - xo) - L;
  double ax = -(k/m)*dx;
  
  //update v
  b.vx += ax * dt;

  //update x
  b.x += b.vx * dt;
}

void render(Ball &b, float t) {

  if (t - floor(t) < 0.1) {
    cout << b.x << endl;
  }

  /* this might actually be entirely separate
     eg, I could have this entire program output a script of positions(t)
     and have That be an input to a separate rendering program

     this is advantageous if the physics takes a long time to churn.
  */
}


int main() {

  init();

  cout << "hello world" << endl;
  cout << balls[1].y << endl;

  Spring* soloSpring = new Spring(0,1);
  double Energy = soloSpring->Energy(balls);

  cout << Energy << endl;

  FILE* fout;
  fout = fopen("positions.txt", "w");

  float t=0, dt=0.1;
  while (t<100) {
    physicsTest(balls[1], dt);
    render(balls[1], t);

    fprintf(fout, "%f\n", balls[1].x);
    t += dt;
  }

  fclose(fout);
  return 0;
}
