#include "gr.h"
#include "gr3.h"

#include "stdio.h"

#include <iostream>
#include <fstream>
#include <vector>

using namespace std;

int main(void) {

  double num;
  vector<double> x(1024);

  ifstream infile("newPositions.txt", ios::in);

  int T,N,X;
  infile >> T;
  infile >> N;
  X = 3*N;
  float in[T][X];
  for(int t=0; t<T; t++) {
    for(int x=0; x<X; x++) {
      infile >> num;
      in[t][x] = num;
    }
  }

  float *positions;
  float colors[X];  
  float radii[N];

  for(int i=0; i<N; i++) {
    colors[3*i + 0] = 1;
    colors[3*i + 1] = 0;
    colors[3*i + 2] = 0;

    radii[i] = 0.25;
  }

  /*
    float cyl_colors[] = {1,1,1, 1,1,1};
    float cyl_radius[] = {0.1, 0.1};
    float cyl_dir[] = {1,0,0, 1,0,0};
    float lengths[] = {in[t][3] - in[t][0],
		       in[t][6] - in[t][3]};
    //this is savage. what if non adjacent?
   */

  ifstream inf2("newSprings.txt", ios::in);
  int S,Y,V;
  inf2 >> S;
  Y = 7*S;
  V = 3*S; //bad names
  float in2[T][Y];
  for(int t=0; t<T; t++) {
    for(int y=0; y<Y; y++) {
      inf2 >> num;
      in2[t][y] = num;
    }
  }

  float inX[T][V], inD[T][V], inL[T][S];
  for(int t=0; t<T; t++) {
    for(int s=0; s<S; s++) {
   
      inX[t][3*s + 0] = in2[t][7*s + 0];
      inX[t][3*s + 1] = in2[t][7*s + 1];
      inX[t][3*s + 2] = in2[t][7*s + 2];

      inD[t][3*s + 0] = in2[t][7*s + 3];
      inD[t][3*s + 1] = in2[t][7*s + 4];
      inD[t][3*s + 2] = in2[t][7*s + 5];

      inL[t][s] = in2[t][7*s + 6];
    }
  } //maybe there is a more elegant way to do this, or is this very elegant?

  float *cyl_pos;
  float *cyl_dir;
  float *cyl_len;

  float cyl_col[V];
  float cyl_rad[S];
  for(int j=0; j<S; j++) {
    cyl_rad[j] = 0.1;

    cyl_col[3*j + 0] = 1;
    cyl_col[3*j + 1] = 1;
    cyl_col[3*j + 2] = 1;
  }


  for(int t=0; t<T; t++) {

    positions = in[t];
    cyl_pos = inX[t];
    cyl_dir = inD[t];
    cyl_len = inL[t];
    
    gr3_clear();
    gr3_drawspheremesh(N, positions, colors, radii);
    gr3_drawcylindermesh(S, cyl_pos, cyl_dir, cyl_col, cyl_rad, cyl_len);

    gr_clearws();
    gr_setviewport(0, 1, 0, 1);
    gr3_drawimage(0, 1, 0, 1, 500, 500, GR3_DRAWABLE_GKS);
    gr_updatews();    

  }

  getc(stdin);
  return 0;
}
