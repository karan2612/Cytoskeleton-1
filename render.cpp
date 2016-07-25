#include "gr.h"
#include "gr3.h"

#include "stdio.h"

#include <iostream>
#include <fstream>
#include <vector>

using namespace std;

void read(){
}

int main(void) {

  double num;
  vector<double> x(1024);

  ifstream infile("newPositions.txt", ios::in);
  string str;

  int T,X;
  infile >> T;
  infile >> X;
  X = 3*X;
  float in[T][X];
  for(int t=0; t<T; t++) {
    for(int x=0; x<X; x++) {
      infile >> num;
      in[t][x] = num;
    }
  }

  float colors[] = {1,0,0, 1,1,1, 1,0,0};
  float radii[] = {0.25, 0.25, 0.25};
  float directions[] = {1, 0, 0};


  float cyl_colors[] = {1,1,1, 1,1,1};
  float cyl_radius[] = {0.1, 0.1};
  float cyl_dir[] = {1,0,0, 1,0,0};


  float* positions;
  for(int t=0; t<T; t++) {
    float lengths[] = {in[t][3] - in[t][0],
		       in[t][6] - in[t][3]};
    //this is savage. what if non adjacent?

    positions = in[t];
    
    gr3_clear();
    gr3_drawspheremesh(3, positions, colors, radii);
    gr3_drawcylindermesh(2, positions, cyl_dir,
			 cyl_colors, cyl_radius, lengths);

    gr_clearws();
    gr_setviewport(0, 1, 0, 1);
    gr3_drawimage(0, 1, 0, 1, 500, 500, GR3_DRAWABLE_GKS);
    gr_updatews();    

  }

  getc(stdin);
  return 0;
}
