/*
  This file contains the following Initialization functions
    meshInit()
    springInit()
    elemInit()
    edgeInit()
    spectrinInit()

    fileInit()
    fileClose()
    initPID()
    randInit()
    initParticle()

  and the retired functions as well
    hexInit()
    toyInit()
 */

void meshInit() {

  int N=_nSys;
  double L = _lActin;
  double h = sqrt(3)*L/2;

  double a = -(N - 1/2.);
  double x0 = L*a;
  double y0 = h*a;
  double z0 = 0;


  /* Build Balls */
  Ball *b;
  for (int y=0; y<2*N; y++) {
    for (int x=0; x<2*N; x++) {

      if (y%2 == 0) //even
	{ b = new Ball(x0 + x*L, y0 + y*h, z0); }    
      else //odd
	{ b = new Ball(x0+L/2 + x*L, y0 + y*h, z0); }

      if (isEdge(x,y,N)) {
	b->isEdge = true; 
	b->pid=2;

	b->edgeType = edgeType(x,y,N);
      }

      v_balls.push_back(*b);
    }
  }

  springInit();

  cout << "   mesh complete" << endl;
}


void springInit() {

  /* Build Springs */
  int j;
  int N = _nSys;
  int n = 2*N;
  for (int y=0; y<n; y++) {
    for (int x=0; x<n; x++) {

      j = n*y + x;

      if (x != n-1) {
	Spring s(j,j+1);
	v_springs.push_back(s);
      }

      if (y == n-1) continue;
      Spring s1(j, j+n);
      v_springs.push_back(s1);

      if (y%2 == 0) { /* y even */
	if (x == 0) continue;
	Spring s0(j, j+n - 1);
	v_springs.push_back(s0);

      } else { /* y odd */
	if (x == n-1) continue;
	Spring s0(j, j+n + 1);
	v_springs.push_back(s0);
      }
      
    }
  }
  cout << "proto-spring complete" << endl;
}


/* between mesh and spectrin, loop through existing balls 
   and write nodes 

   e1 and e0 stand for elemUP and elemDOWN respectively,
  where an element is the triangular surface element connecting 3 balls
*/
void elemInit() {

  Elem *e0, *e1;

  /* borrowed from Build Springs */
  int j;
  int N = _nSys;
  int n = 2*N;
  for (int y=0; y<n; y++) {
    for (int x=0; x<n; x++) {

      if (x%n == 0) continue;
      j = y*n + x;

      if (y%2 == 0) {
	e0 = new Elem(j,j+1,j+1-n);

	if (y==n) continue;

	e1 = new Elem(j,j+1,j+1+n);

      } else {
	e1 = new Elem(j,j+1,j+n);

	if (y==0) continue;
	e0 = new Elem(j,j+1,j-n);
      }

      if (e0) v_elems.push_back(*e0);
      if (e1) v_elems.push_back(*e1);
    }
    
  }

}


/* Build Edges, same basic algorithm as initSprings */
void edgeInit() {

  int j;
  int N = _nSys;
  int n = 2*N;
  for (int y=0; y<n; y++) {
    for (int x=0; x<n; x++) {

      j = n*y + x;

      if (x != n-1) {
	Edge e(j,j+1);
	v_edges.push_back(e);
      }

      if (y == n-1) continue;
      Edge e1(j, j+n);
      v_edges.push_back(e1);

      if (y%2 == 0) { /* y even */
	if (x == 0) continue;
	Edge e0(j, j+n - 1);
	v_edges.push_back(e0);

      } else { /* y odd */
	if (x == n-1) continue;
	Edge e0(j, j+n + 1);
	v_edges.push_back(e0);
      }
      
    }
  }
  /* still need to identify two adjacent elems!! */
  cout << "edge init complete" << endl;
}


/*  int n = _nSpectrin; //new springs :: 1 should be identity */
void spectrinInit(int n) {

  //makes a copy
  //vector<Spring> vs = v_springs;
  int N = v_springs.size();
  if (N < 1) {
    cout << "Error: spring vector is empty!" << endl;
    getc(stdin);
  }

  Spring *spr;
  Ball *b1, *b2;

  double r1[2],r2[2],dr[2]; //lives in 2D plane
  double x,y,z;

  double A = _lActin;
  double C = _Contour;
  double h = sqrt(C*C + A*A) / 2; //creates dagger strings

  vector<Spring> newSpring_v;
  //for each spring
  int j1,j2;
  for (int j=0; j<N; j++) {
    spr = &v_springs.at(j);

    j1 = spr->n1; 
    j2 = spr->n2; 

    b1 = &v_balls.at(j1);
    b2 = &v_balls.at(j2);

    for (int i=0; i<2; i++) {
      r1[i] = b1->r[i];
      r2[i] = b2->r[i];
      dr[i] = (r2[i] - r1[i]) / n;
    }

    //create n-1 new balls
    int k = v_balls.size(); //used later
    double a = A/2;
    for (int i=1; i<n; i++) {
      x = r1[0] + i*dr[0];
      y = r1[1] + i*dr[1];
      z = -h + h/a * abs( a - (i*A/n) );

      Ball b(x,y,z);
      b.pid = 0;
      v_balls.push_back(b);
    }

    //connect w N new springs (+2 old balls)
    if (n == 1) return;
    newSpring_v.push_back( Spring(j1,k) ); 
    for (int i=1; i<n-1; i++) 
    {
      newSpring_v.push_back( Spring(k, k+1) );
      k++;
    }
    newSpring_v.push_back( Spring(k, j2) );

  }

  //haben Sie memory leak??
  v_springs = newSpring_v;

  cout << "   spectrin spring complete " << n << endl;
       
}

void fileInit() {

  cout << " initializing file.." << endl;

  // is it ineffective to do all of these all the time?
  f1 = fopen("balls.dat", "w");
  f2 = fopen("springs.dat", "w");
  f3 = fopen("force.txt", "w");
  f4 = fopen("contour.txt", "w");
  f5 = fopen("part.txt", "w");
  f6 = fopen("zForce.txt", "w");
  f7 = fopen("WrappingEnergy.txt", "w");

  /* Writes nTime, nBalls, and nSprings */
  nBalls = v_balls.size();
  nSprings = v_springs.size();

  //  nTime = (int) tmax/dt;
  nTime = nTime/ts + 1;

  nTime = (int) nSteps/ts;
  cout << "   predicted time steps: " << nTime << endl;
  cout << "   physics steps per time step: " << ts << endl;
  cout << "   TOTAL differential steps: " << nSteps << endl;

  fprintf(f1, "%i\n" , nTime);
  fprintf(f1, "%i\n\n" , nBalls);
  fprintf(f2, "%i\n\n" , nSprings);

  initPID();
}

void filesClose() {
  fclose(f1);
  fclose(f2);
  fclose(f3);
  fclose(f4);
  fclose(f5);
  fclose(f6);
  fclose(f7);

  cout << "files written." << endl;
}

/* writes PID */
void initPID() {

  int N = nBalls;
  for(int j=0; j<N; j++) { 

    //write pid (1-anchor, 0-spectrin)
    fprintf(f1, " %i", v_balls[j].pid);
  }
  fprintf(f1, "\n\n");

}

void randInit() {

  cout << " initializing rand.." << endl;
  unsigned int saat = (unsigned int)time(0); 
  randi.seed(saat);
}


void hexInit() {

  int i = 1;
  Ball a0(0.001, 0 ,i);
  v_balls.push_back(a0);

  double sd = PI/3; //'Sixty Degrees'

  /* Can I automate this in a 2line for loop? */
  Ball a1(1, 0*sd, i);
  Ball a2(1, 1*sd, i);
  Ball a3(1, 2*sd, i);
  Ball a4(1, 3*sd, i);
  Ball a5(1, 4*sd, i);
  Ball a6(1, 5*sd, i);

  v_balls.push_back(a1);
  v_balls.push_back(a2);
  v_balls.push_back(a3);
  v_balls.push_back(a4);
  v_balls.push_back(a5);
  v_balls.push_back(a6);

  /* now for springs */
  Spring s01 = Spring(0,1);
  Spring s02 = Spring(0,2);
  Spring s03 = Spring(0,3);
  Spring s04 = Spring(0,4);
  Spring s05 = Spring(0,5);
  Spring s06 = Spring(0,6);

  Spring s12 = Spring(1,2);
  Spring s23 = Spring(2,3);
  Spring s34 = Spring(3,4);
  Spring s45 = Spring(4,5);
  Spring s56 = Spring(5,6);
  Spring s61 = Spring(6,1);

  v_springs.push_back(s01);
  v_springs.push_back(s02);
  v_springs.push_back(s03);
  v_springs.push_back(s04);
  v_springs.push_back(s05);
  v_springs.push_back(s06);
  v_springs.push_back(s12);
  v_springs.push_back(s23);
  v_springs.push_back(s34);
  v_springs.push_back(s45);
  v_springs.push_back(s56);
  v_springs.push_back(s61);

}

void toyInit() {
  Ball a0(-1,0,0.);
  Ball a1(0,0,0.);
  Ball a2(1,0,0.);
  v_balls.push_back(a0);
  v_balls.push_back(a1);
  v_balls.push_back(a2);

  Spring s1 = Spring(0,1);
  Spring s2 = Spring(1,2);
  v_springs.push_back(s1);
  v_springs.push_back(s2);

}

/* Equilibriates the system, 
   or add desired displacement */
void sysConfig() {
  
  int N = v_balls.size();
  int j;
  for (j=0; j<N; j++) {

    v_balls[j].pid = 1;
    v_balls[j].isEdge = false;
  }

  int n = v_springs.size();
  float f = 0.4;
  for (j=0; j<n; j++) {
       v_springs[j].eql = f * v_springs[j].L;
  }

}
