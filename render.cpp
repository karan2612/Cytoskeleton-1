#include "gr.h"
#include "gr3.h"
#include "stdio.h"

#include <iostream>
#include <fstream>
#include <vector>

using namespace std;

int main(void) {
  
  ifstream infile("positions.txt", ios::in);
  vector<double> x;

  double num  = 0.0;
  while ( infile >> num ) {
    x.push_back(num);
  }

  float positions[] = {-1, 0, 0, 1, 0, 0};

  for(size_t i =0; i < x.size(); i++) {

    float colors[] = {1, 1, 1, 1, 0, 0};
    float radii[] = {0.25, 0.25};
    float directions[] = {1, 0, 0};
    float lengths[] = {positions[3] - positions[0]};
    float cyl_radius[] = {0.1};

    //    cout << x[i] << endl;

    positions[3] = x[i];

    gr3_clear();
    gr3_drawspheremesh(2, positions, colors, radii);
    gr3_drawcylindermesh(1, positions, directions, colors, cyl_radius, lengths);

    gr_clearws();
    gr_setviewport(0, 1, 0, 1);
    gr3_drawimage(0, 1, 0, 1, 500, 500, GR3_DRAWABLE_GKS);
    gr_updatews();    
  }
  getc(stdin);
  //gr3_export("test.html", 1000, 1000);                                                                                 
  return 0;
}
