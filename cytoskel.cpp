#include "cytoskel.h"

int main() {

  cout << "hello world!" << endl;

  init();
  //  toyInit();

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
}

void physics() {

  /* Defined in Header:
     tmax - (float) absolute stop time
     dt - (float) time step physics
     ts - (int) time step output */

  cout << " beginning physics.." << endl;
  int t=0, t_count=0;

  float T=0; 
  while (T<tmax) {

    timeStep();

    if (t % ts == 0) 
    { t_count++;
      
      /* make measurements */
      writeBalls(f1);
      writeSprings(f2);
      measureEdge(f3);

      //  if (_msd) fprintf(f3, "%f\n", msd);
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

  /* Tally Forces */
  ForceSprings();
  SurfaceForce();
  //ParticleInteraction();

  /* Update Particles */
  for (size_t j=0; j<N; j++) {

    if (v_balls[j].isEdge) continue; 
    //updatePosition(v_balls[j]); 
    updateBrownianPosition(v_balls[j]); 
  }

}


/* loop all springs to compute force on all balls */
void ForceSprings() {

  Spring *spr;
  Ball *b1, *b2;

  int j1,j2; //vector index for the balls
  double L, X, K, F;
  double u[3], v[3], dr[3];

  size_t N = v_springs.size();
  for (size_t j=0; j<N; j++) {

    spr = &v_springs[j];

    L = spr->getLength();
    X = spr->eql;
    K = spr->k;
    F = - K * (L - X);

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

/* repels spectrin from surface, tethers actin to surface*/
void SurfaceForce() {
  int N = nBalls;
  Ball *b;

  int isActin;  
  for (size_t j=0; j<N; j++) {
    b = & v_balls[j];
    isActin = b->pid; 

    if (isActin) {
      b->r[2] = zSurface(b->r[0],
			 b->r[1]); //better way?
      continue;
    } 

    //else: is spectrin
    double x,y,z,z0;  
    double r,F;
    x = b->r[0];
    y = b->r[1];
    z = b->r[2]; //should be negative
    z0= zSurface(x,y);    

    if (z > z0) {
      cout << "Warning: z > z0!" << endl;
    }

    r = abs(z - z0);
    F = LJforce(r);
    b->F[2] += F; //sign?
  }
}


void updateBrownianPosition(Ball &b) {

  double D = 0.01;
  double m = b.m;

  for (int i=0; i<3; i++) {
    b.r[i] += b.F[i]/m * dt;
    b.r[i] += sqrt(2*D*dt)*randi.randNorm(0,1); 
  }

}

void updatePosition(Ball &b) {

  double a[3];
  double m = 1; //temp fix
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



