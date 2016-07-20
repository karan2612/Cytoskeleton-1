
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
