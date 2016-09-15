#include "global.h"
#include "init.h"
#include "physics.h"
#include "analysis.h"

#include "inputModule.h"

int main() {

  cout << "hello world!" << endl;

  init();

  physics();

  filesClose();
  return 0;
}

void init() {

  readInput();
  //getc(stdin);

  /* Build Cytoskeleton System */
  cout << " initializing system.." << endl;
  meshInit();
  spectrinInit(_nSpectrin);
  initParticle();

  /* Set Files and Rand*/
  randInit();
  fileInit();
}

void physics() {

  /* Defined in Header:
     tmax - (float) absolute stop time
     _dt - (float) time step physics
     _tSamp - (int) time step output */

  cout << " beginning physics.." << endl;

  float T=0; //elapsed time
  int t=0, t_count=0;
  //  while (T<tmax) {
  while (t < _tMax) {

    timeStep(t);
    t++;
    T += _dt;

    //init equilibriation
    if (t < _tEqGlobal) continue;

    if (t % _tSamp == 0) { 

      t_count++;
      doAnalysis();

      cerr << t_count << "\r";
    }
  }

  cout << "   rendered time steps: " << t_count << endl;
  cout << "   final system time : " << T << endl;
  
  integrateWrappingEnergy();
  //  writeForce3D();

  cout << " finished physics.." << endl;
}


/* Everything here called once per Physics time step */
void timeStep(int t) {

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
  
  if ( (t % _tSamp) > _tEqLocal) {
    sampleForceZ();
  }

}

void doAnalysis() {
 
  /* make measurements */
  sampleForce3D();
  writeForceZ(f6); //mean Fz discovered in here!

  Particle->r[2] -= _dz; //moves particle

  /* for rendering */
  writeBalls(f1);
  writeSprings(f2);
}
 
