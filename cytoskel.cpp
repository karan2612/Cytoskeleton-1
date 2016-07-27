#include "cytoskel.h"

int main() {

  cout << "hello world!" << endl;

  /* Begin Init */
  cout << " initializing system.." << endl;
  //hexInit();
  meshInit();

  cout << " initializing file.." << endl;
  f1 = fopen("newBalls.txt", "w");
  f2 = fopen("newSprings.txt", "w");

  int T, B, S; 
  T = (int)tmax/dt;
  B = ball_v.size();
  S = spring_v.size();

  fprintf(f1, "%i\n" , T);
  fprintf(f1, "%i\n\n" , B);
  fprintf(f2, "%i\n\n" , S);

  cout << " initializing rand.." << endl;
  unsigned int saat = (unsigned int)time(0); //ok?
  randi.seed(saat);

  // init msd 
  int N = ball_v.size();
  float x0[N];
  for(int i=0; i<N; i++) {
    x0[i] = ball_v.at(i).x;
  }
  if (_msd) f3 = fopen("outMSD.txt", "w");
  /* End Init */
  
  cout << " beginning physics.." << endl;
  float t=0; int ts=0; 
  while (t<tmax) {

    t += dt;
    ts++;

    if (_msd) {
      msd = calcMSD(x0);
      fprintf(f3, "%f\n", msd);
      continue;
    }
    physics();
    //getc(stdin);

    writePositions(f1);
    writeSprings(f2);
  }

  cout << " finished physics.." << endl;
  fclose(f1);
  fclose(f2);

  cout << "files written." << endl;
  return 0;
}

/* new development: no spring, just brown*/
float calcMSD(float *x0) {

  size_t N = ball_v.size(); 

  float xt[N];
  float sum=0;  

  for (size_t j=0; j<N; j++) {
    /* Zeros Forces */
    ball_v[j].Fx = 0;
    ball_v[j].Fy = 0;

    /* update physics (no springs) */
    updateBrownianPosition(ball_v[j]); 

    /* compute MSD */
    xt[j] = ball_v.at(j).x;
    sum += (xt[j] - x0[j]) * (xt[j] - x0[j]);
  }

  return sum;
}

void initMSD() {

}


/* Everything here should be called once per time step*/
void physics() {

  size_t N = ball_v.size(); 

  /* Zeros Forces */
  for (size_t j=0; j<N; j++) {
    ball_v[j].F  = 0;
    ball_v[j].Fx = 0;
    ball_v[j].Fy = 0;
  }

  /* add springs */
  ForceSprings();

  /* add brownian motion */
  if (0) {

  double D = 0.1;

  for (size_t j=0; j<N; j++) {
    ball_v[j].Fx += sqrt(2*D/dt)*randi.randNorm(0,1);  //dividing by sqrt(dt) bc a_x needs to mult sqrt(dt)
    ball_v[j].Fy += sqrt(2*D/dt)*randi.randNorm(0,1);
  }

  }

  /* loop particles */

  for (size_t j=0; j<N; j++) {
    //updatePosition(ball_v[j]); 
    updateBrownianPosition(ball_v[j]); 
  }

}


/* loop all springs to compute force on all balls */
void ForceSprings() {

  Spring *spr;
  Ball *b1, *b2;

  int j1,j2; //vector index for the balls
  double L, X, F;
  double x1, x2, y1, y2;

  size_t N = spring_v.size();
  for (size_t j=0; j<N; j++) {

    spr = &spring_v[j];

    L = spr->Length();
    X = spr->eql;
    F = -k * (L - X);

    if (0) {
      cout << "Spring " << j << endl;
      cout << "(L, X, F): "
	   << L << ", "
	   << X << ", "
	   << F 
	   << endl;
    }

    j1 = spr->n1;
    j2 = spr->n2;
    b1 = &ball_v[j1];
    b2 = &ball_v[j2];

    x1 = b1->x;
    x2 = b2->x;

    b1->Fx += F * (x1-x2) / L; 
    b2->Fx -= F * (x1-x2) / L;

    y1 = b1->y;
    y2 = b2->y;

    b1->Fy += F * (y1-y2) / L; 
    b2->Fy -= F * (y1-y2) / L;

  }

}

void updateBrownianPosition(Ball &b) {

  double D = 0.1;

  //update x: see wiki
  b.x += (b.Fx/m)*dt + sqrt(2*D*dt)*randi.randNorm(0,1); 
  b.y += (b.Fy/m)*dt + sqrt(2*D*dt)*randi.randNorm(0,1); 

}
void updatePosition(Ball &b) {

  //compute a
  double ax = b.Fx / m;
  double ay = b.Fy / m;
  
  //update v
  b.vx += ax * dt;
  b.vy += ay * dt;

  //update x
  b.x += b.vx * dt;
  b.y += b.vy * dt;

  if (0) {
    cout << b.x << ", "
	 << b.vx << ", "
	 << ax << endl;

    cout << b.y << ", "
	 << b.vy << ", "
	 << ay << endl;
  }
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


void writeSprings(FILE* f) {

  Spring *spr;
  Ball *b1, *b2;
  int n = spring_v.size();
  float dx,dy,dz;
  float L;
  float nx,ny,nz;

  /* fetch information about position, direction, L */
  for(int j=0; j<n; j++) {

    spr = &spring_v.at(j);
    b1 = &ball_v.at(spr->n1);
    b2 = &ball_v.at(spr->n2);

    dx = b2->x - b1->x;
    dy = b2->y - b1->y;
    dz = b2->z - b1->z;

    L = sqrt( dx*dx + dy*dy + dz*dz );    

    nx = dx/L;
    ny = dy/L;
    nz = dz/L;
    
    fprintf(f, "%f" , b1->x);
    fprintf(f, " %f", b1->y);
    fprintf(f, " %f", b1->z);

    fprintf(f, " %f", dx);
    fprintf(f, " %f", dy);
    fprintf(f, " %f", dz);

    fprintf(f, " %f\n", L);
  }
  fprintf(f, "\n");

}
