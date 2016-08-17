/*
  This file contains the following Physics
    distBall(Ball*, Ball*)
    normBall(Ball*, Ball*)
    magForce(Ball*)

    forceSprings()
    surfaceForce()
    updateBrownianPosition(Ball &)
    updatePosition(all &)
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

double magForce(Ball *b) {

  double F=0, f=0;
  for (int i=0; i<3; i++) {
    f = b->F[i];
    F += f*f;
  }
  return sqrt(F);
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
    F = LJforce(r,0.1); //hardcode surface sigma
    b->F[2] += -F; //sign?
  }
}


void updateBrownian(Ball &b) {

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

double LJforce(double r, double d) {

  /* r_min : F = 0, truncate potential
     sigma : U = 0, is const for all, 
     but it should be the sum of radii d = (a+b) */

  double F, e, m, p;
  e = 0.01;
  //m = 0.45; 
  m = 1.1225 * _sigma;
  r = r - d + _sigma;

  if (r > m) return 0; //yes, U(r') = 0

  /* F = 12 e/r [p12 - p6]  
      where p = r_min/r */
  p = m/r;
  if (1/p < 0.05) 
    cout << "Warning: LJ Singularity from small r" << endl;

  p = pow(p,6);
  F =  e /r * p * (p - 1);

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
 
  double radius = Particle->R;
  double r,d;
  vector<float> nm;
  double F;

  int N = nBalls;
  int count = 0;
  Ball *b;
  //  cout << "***time step***" << endl;
  for (int j=0; j<N; j++) {
    
    b = & v_balls.at(j);
    r = distBall(b, Particle);
    nm = normBall(b,Particle); // points b->P

    if ( r < radius ) 
    {

      F = -5; //F- rides the nm from P->b
      for (int i=0; i<3; i++) {
	b->F[i] += F * nm[i];
      }

    } else 
    {
      if (r > 1.5 * radius) {
	//somehow count skips?
	continue;
      }
      count++;

      d = (_sigma/2) + radius; // R_part + R_interactant
      F = LJforce(r,d);
      //      cout << F << endl;
      for (int i=0; i<3; i++) {
	b->F[i]        -= F *  nm[i];
	Particle->F[i] += F *  nm[i];

       	//cout << nm[i] << endl;
      }

      //investigate normals
      float rt = sqrt( nm[0]*nm[0] + nm[1]*nm[1] );
      float zt = nm[2];
      float tan= zt/rt;
      //      cout << tan << endl;
      //      cout << endl;      
    }

  }   

  //end of interactant loop

}
