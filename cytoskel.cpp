#include "cytoskel.h"

int main() {

  cout << "hello world!" << endl;

  init();
  physics();

  filesClose();
  return 0;
}

void init() {

  /* Build Cytoskeleton System */
  cout << " initializing system.." << endl;
  meshInit();

  cout << v_balls.size() << " "<< v_springs.size() << endl;
  v_springs = subInit(v_springs);
  cout << v_balls.size() << " "<< v_springs.size() << endl;

  /* Set Files and Rand*/
  fileInit();
  randInit();
  initKaran();

}

void physics() {

  /* for all balls (post init) */
  initBalls();

  /* Defined in Header:
     tmax - (float) absolute stop time
     dt - (float) time step physics
     ts - (int) time step output
   */

  cout << " beginning physics.." << endl;
  int t=0, t_count=0;

  float T=0; 
  while (T<tmax) {

    timeStep();

    if (t % ts == 0) 
    {
      t_count++;
      if (_msd) 
      {
	//msd = calcMSDx(x0);
	fprintf(f3, "%f\n", msd);
	//continue;
      }

      writeBalls(f1);
      writeSprings(f2);
      //  writeKaranXY(kx,ky);
    }

    t++;
    T += dt;
  }

  cout << "actual time steps: " << t_count << endl;
  cout << " finished physics.." << endl;
}


/* Everything here should be called once per time step*/
void timeStep() {

  size_t N = nBalls;

  /* Zeros Forces */
  for (size_t j=0; j<N; j++) {
    for(int i=0; i<3; i++) {
      v_balls[j].F[i] = 0;
    }
  }

  /* add springs */
  ForceSprings();
  BoundarySprings();

  /* loop particles */
  int isSpectrin;
  for (size_t j=0; j<N; j++) {

    //    if (v_balls[j].isCorner) continue; 
    if (v_balls[j].isEdge) continue; 

    //updatePosition(v_balls[j]); 
    isSpectrin = v_balls[j].pid;
    updateBrownianPosition(v_balls[j]); 
  }
}

void BoundarySprings() {

  //identify boundary balls

  //for each find boundary partners

  //calc periodic distance

  //compute force
}

/* loop all springs to compute force on all balls */
void ForceSprings() {

  Spring *spr;
  Ball *b1, *b2;

  int j1,j2; //vector index for the balls
  double L, X, F;
  double u[3], v[3], dr[3];

  size_t N = v_springs.size();
  for (size_t j=0; j<N; j++) {

    spr = &v_springs[j];

    L = spr->Length();
    X = spr->eql;
    F = -k * (L - X);

    j1 = spr->n1;
    j2 = spr->n2;
    b1 = &v_balls[j1];
    b2 = &v_balls[j2];

    for(int i=0; i<3; i++) {
      u[i] = b1->r[i];
      v[i] = b2->r[i];

      dr[i] = u[i]-v[i];
      b1->F[i] += F * dr[i] / L;
      b2->F[i] -= F * dr[i] / L;
    }

  }
}

void updateBrownianPosition(Ball &b) {

  double D = 0.01;

  for (int i=0; i<3; i++) {
    b.r[i] += b.F[i]/m * dt;
    b.r[i] += sqrt(2*D*dt)*randi.randNorm(0,1); 
  }

}

void updatePosition(Ball &b) {

  double a[3];
  for (int i=0; i<3; i++) {
    a[i] = b.F[i]/m;

    b.v[i] += a[i] * dt;
    b.r[i] += b.v[i] * dt;
  }
}


void writeBalls(FILE* f) {

  int n = v_balls.size();
  for(int j=0; j<n; j++) {
    for(int i=0; i<3; i++) {
      fprintf(f, "%f ", v_balls[j].r[i]);   
    }
    fprintf(f, "\n");
  }

  fprintf(f, "\n"); //new timestep
}


void writeSprings(FILE* f) {

  Spring *spr;
  Ball *b1, *b2;
  int n = v_springs.size();

  float L;
  float dx,dy,dz;
  float nx,ny,nz;
  float r1[3],r2[3],dr[3],nr[3];

  /* fetch information about position, direction, L */
  for(int j=0; j<n; j++) {

    spr = &v_springs.at(j);
    b1 = &v_balls.at(spr->n1);
    b2 = &v_balls.at(spr->n2);

    L = 0;
    for (int i=0; i<3; i++) {
      
      r1[i] = b1->r[i];
      r2[i] = b2->r[i];
      dr[i] = r2[i] - r1[i];
      
      L += dr[i]*dr[i];
    }
    L = sqrt(L);

    for (int i=0; i<3; i++) {
      fprintf(f, "%f " , r1[i]);
    }
    
    for (int i=0; i<3; i++) {
      nr[i] = dr[i]/L; //necessary?
      fprintf(f, "%f ", nr[i]);
    }

    fprintf(f, " %f\n", L);
  }

  fprintf(f, "\n");
}



