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
  randInit();
  fileInit();
}

void physics() {

  /* Defined in Header:
     tmax - (float) absolute stop time
     dt - (float) time step physics
     ts - (int) time step output */

  cout << " beginning physics.." << endl;

  float T=0; 
  int t=0, t_count=0;
  //  while (T<tmax) {
  while (t < nSteps) {

    timeStep(t);
    t++;
    T += dt;

    if (t < 500) continue;
    if (t % ts == 0) { 

      t_count++;
      doAnalysis();

      cerr << t_count << "\r";
    }
  }

  cout << "   actual time steps: " << t_count << endl;
  cout << "   final system time : " << T << endl;
  
  integrateWrappingEnergy();
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
  
  if (t % ts < 10) return;
  sampleForceZ();

}

void doAnalysis() {
 
  /* make measurements */
  sampleForce3D();
  writeForceZ(f6); //mean Fz discovered in here!
  moveParticle();

  /* for rendering */
  writeBalls(f1);
  writeSprings(f2);
}


void initParticle() {

  /* Let wrapping fraction c (0,2)
     then z = c * R (if z = 0, c=0, at psi=0)
   */
  double z,c,R;
  double buffer;
  
  c = -0.2; //wrapping fraction
  R = 1.1 * _lActin; //input radius
  z = R * (1-c);

  //  Particle = new Ball(0,0,z);
  Particle = new Ball(0,0,z);
  Particle->R = R;

}

/* somehow need to equilibriate system before this kicks in.
   maybe it would be useful to make global the ellapsed time and step time, 
   for easier access by skip logic
*/
void moveParticle() {

  Particle->r[2] -= _dz;
  //Particle->r[2] += _dz; //moving up

}

 
