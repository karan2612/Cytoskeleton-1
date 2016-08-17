#include "global.h"
#include "init.h"
#include "physics.h"
#include "analysis.h"

int main() {

  cout << "hello world!" << endl;

  init();

  physics();
  writeForce3D();

  filesClose();
  return 0;
}

void init() {

  /* Build Cytoskeleton System */
  cout << " initializing system.." << endl;
  meshInit();
  spectrinInit(_nSpectrin);
  initParticle();

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
    t++;
    T += dt;

    if (T < 500) continue;
    if (t % ts == 0) { 

      t_count++;
      doAnalysis();

      cerr << t_count << "\r";
    }
  }
  cout << "   actual time steps: " << t_count << endl;
  cout << " finished physics.." << endl;
}


/* Everything here should be called once per time step*/
void timeStep() {

  size_t N = nBalls;

  /* Zeros Forces */
  for(int i=0; i<3; i++) {
    for (size_t j=0; j<N; j++) {
      v_balls[j].F[i] = 0;
    }
    Particle->F[i] = 0;
  }


  /* Tally Forces */
  ForceSprings();
  SurfaceForce();
  ParticleInteraction();

  /* Update Particles */
  for (size_t j=0; j<N; j++) {

    if (v_balls[j].isEdge) continue; 
    //updatePosition(v_balls[j]); 
    updateBrownian(v_balls[j]); 
  }
  
  sampleForceZ();
}


void sampleForceZ() {

  float z = Particle->F[2];
  statForce.push_back(z);
}

void moveParticle() {

  Particle->r[2] -= 0.013;

}
void doAnalysis() {
 
  /* make measurements */
  // measureEdge(f3);
  //measureSpringEnergy();
  //measureContour(f4);
  sampleForce3D();
  moveParticle();
  //  cout << Particle->r[2] << endl;

  //  if (_msd) fprintf(f3, "%f\n", msd);
  //  writeKaranXY(kx,ky);

  /* for rendering */
  writeBalls(f1);
  writeSprings(f2);
}


 
