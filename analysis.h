/*
  This file contains the following Analysis functions
    measureEdge()
    measureEnergy()
    measureContour()
    show()

  and the retired msd functions
    writeKaranXY()
    calcMSD()
    calcMSDx()
 */

void measureEdge(FILE* f) {
  int N = nBalls;
  Ball *b;
  double F;
  for (size_t j=0; j<N; j++) {
    b = &v_balls.at(j);
    if (b->edgeType == 1) {
      //cout << b->getForce() << endl;
      F = b->getForce();
      fprintf(f, "%f ", F);
      //hmmm there's probably a c++ way to do this
      //      f << F;
    }
  }

  //  getc(stdin);
  //  cout << endl;
  fprintf(f, "\n");
}

void measureEnergy() {

  int N = nSprings;

  Spring *s;
  float E = 0;
  for (int j=0; j<N; j++) {

    s = &v_springs[j];
    E += s->getEnergy();
  }

  //should I normalize this?
  //  cout << E << endl;
}

void measureContour(FILE *f) {

  /* the springs are created N at a time
     so I can exploit this in reading the segments */

  int N = nSprings;
  int n = _nSpectrin;

  Spring *s;
  float C=0; //contour
  for(int j=0; j<N; j++) {

    s = &v_springs.at(j);
    C += s->getLength();

    if (j == 0) continue;
    if (j%n == 0) {
      //cout << C << endl;
      fprintf(f, "%f ", C);
      C = 0;
    }
  }
  //cout << endl;
  fprintf(f,"\n");
}

/* for misc purposes */


void writeKaranXY(FILE* fx, FILE* fy) {

  size_t N = v_balls.size();
  float x=0;
  float y=0;
  for (int j=0; j<N; j++) {
    x = v_balls[j].r[0];
    fprintf(fx,"%f ",x);

    y = v_balls[j].r[1];
    fprintf(fy,"%f ",y);
  }
  fprintf(fx,"\n");
  fprintf(fy,"\n");
}


/* new development: no spring, just brown*/
float calcMSD(float **r0) {

  size_t N = v_balls.size(); 

  float rt[N][3];
  float sum=0;  

  for (size_t j=0; j<N; j++) {
    /* Zeros Forces */
    for (int i=0; i<3; i++) {
      v_balls[j].F[i] = 0;
    }

    /* update physics (no springs) */
    updateBrownianPosition(v_balls[j]); 

    /* compute MSD */
    for (int i=0; i<3; i++) {
      rt[j][i] = v_balls.at(j).r[i];
      sum += (rt[j][i] - r0[j][i])
   	    *(rt[j][i] - r0[j][i]);
    }
  }

  sum = sum/N; //norm
  return sum;
}

float calcMSDx(float *x0) {

  size_t N = v_balls.size(); 

  float xt[N];
  float sum=0;  

  for (size_t j=0; j<N; j++) {
    /* Zeros Forces */
    for (int i=0; i<3; i++) {
      v_balls[j].F[i] = 0;
    }

    /* update physics (no springs) */
    updateBrownianPosition(v_balls[j]); 

    /* compute MSD */
    xt[j] = v_balls.at(j).r[0];
    sum += (xt[j] - x0[j])
          *(xt[j] - x0[j]);
  }
  //cout << sum << " / " << N;
  sum = sum / N; //norm
  //cout << " = " << sum << endl;
  return sum;
}



void show() {
  cout << "showing..." << endl;

  int n = v_balls.size();

  Ball* b;
  for (int i=0; i<n; i++) {
    b = &v_balls.at(i);

    float m; //magnitutde
    for(int j=0; j<3; j++) {
      m += b->F[i] * b->F[i];
    }
    m = sqrt(m);    
    cout << m << endl;
  }

  //  getc(stdin);
}
