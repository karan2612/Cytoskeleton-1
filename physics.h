/*
  This file contains the following Physics
    distBall(Ball*, Ball*)
    normBall(Ball*, Ball*)
    LJforce(r, sig)
    zSurface(x, y)

    ParticleInteraction()
 */

float distBall(Ball *b1, Ball *b2) {

  float L = 0;
  float r1[3],r2[3],dr[3];
  for (int i=0; i<3; i++) {

    r1[i] = b1->r[i];
    r2[i] = b2->r[i];

    dr[i] = r2[i] - r1[i];

    L += dr[i]*dr[i];
  }
  L = sqrt(L);
  
  return L;
}

//returns norm from origin of a to origin of b
vector<float> normBall(Ball *a, Ball *b) {

  float u[3],v[3],n[3];
  for (int i=0;i<3;i++) {
    u[i] = a->r[i];
    v[i] = b->r[i];

    n[i] = v[i] - u[i];
  }

  float m=0; //magnitutde
  for(int j=0; j<3; j++) {
    m += n[j]*n[j];
  }
  m = sqrt(m);
  /* I could make a subroutine,
     but that would pass an array and return a double */

  vector<float> norm = vector<float>(3);
  for(int j=0; j<3; j++) {
    norm[j] = n[j]/m;
  }

  return norm;
}


double LJforce(double r, double sigma) {

  /* r_min : F = 0, truncate potential
     sigma : U = 0, is sum of radii (a+b) */

  double F, e, m, p;
  e = 0.2; //I choose to absorb 12 into e
  //m = 0.45; 
  m = sigma*1.1225;

  if (r > m) return 0;
  if (r < 0.05) {
    cout << "Warning: LJ Singularity from small r" << endl;
  }

  /* F = 12 e/r [p12 - p6]  
      where p = r_min/r */
  p = m/r;
  p = pow(p,6);
  F =  - e /r * p * (p - 1);

  return F;
}


double zSurface(double x, double y) {
  /* returns Membrane z(x,y) */

  //temp flat
  return 0.1;
}


void ParticleInteraction() {
  //cout << " considering Particle interactions.." << endl;

  if (!Particle) {
    cout << "particle not initialized!" << endl;
    getc(stdin);
  }
 
  double radius = 2.15; //hard code for now
  double r;
  vector<float> nm;
  double F;

  int N = nBalls;
  int count = 0;
  Ball *b;
  for (int j=0; j<N; j++) {
    
    b = & v_balls.at(j);
    r = distBall(b, Particle);

    if ( r < radius ) {
      nm = normBall(b,Particle);
      F = -5;
      // F = LJForce(r);
      for (int i=0; i<2; i++) {
	b->F[i] += F * nm[i];
      }

      count++;
    }
  }   

  //  cout << count << endl;

  //how can I be sure that none are INSIDE the Particle?

}
