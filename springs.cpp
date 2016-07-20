Last login: Tue Jul 19 15:24:54 on ttys001
Tonys-Air:~ q$ python
Python 2.7.12 |Continuum Analytics, Inc.| (default, Jul  2 2016, 17:43:17) 
[GCC 4.2.1 (Based on Apple Inc. build 5658) (LLVM build 2336.11.00)] on darwin
Type "help", "copyright", "credits" or "license" for more information.
Anaconda is brought to you by Continuum Analytics.
Please check out: http://continuum.io/thanks and https://anaconda.org
>>> 250.34+1+57.5+43.94
352.78000000000003
>>> 352.78/4
88.195
>>> 
Tonys-Air:~ q$ cd Projects/Biophysics/Cytoskeleton/
Tonys-Air:Cytoskeleton q$ ls
Makefile	florian.C~	modFlorian.C~	render.h~
a.out		gr.h		notes		test.html
first		gr3.h		positions.txt	tony.C
florian.C	modFlorian.C	render.h	tony.C~
Tonys-Air:Cytoskeleton q$ emacs tony.C





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
