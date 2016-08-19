/*
  This file contains the following Analysis functions
    writeBalls()
    writeSprings()
    measureEdge()
    measureSpringEnergy()
    measureContour()
    show()
 */


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


void measureSpringEnergy() {

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

// to be retired??
/* writes all the contours to file, at a given t */
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

    //    cout << j << endl;
    if (j == 0) continue;
    if ((j+1) % n == 0) {
      //cout << C << endl;
      fprintf(f, "%f ", C);

      C = 0;
    }
  }
  //cout << endl;
  fprintf(f,"\n");
}


/* Samples contour through ALL t
   for single entropic spring */
float sampleContour() {

  int N = v_springs.size();

  Spring *s;
  float C=0; //contour
  for(int j=0; j<N; j++) {
    s = &v_springs.at(j);
    C += s->getLength();
  }

  return C;
}

void sampleForceZ() {

  float z = Particle->F[2];
  statForce.push_back(z);
}


void sampleForce3D() {

  float F=0, Fi;
  for (int i=0; i<3; i++) {
    Fi = Particle->F[i];
    //    cout << Fi << endl;
    forceST.at(i).push_back(Fi);
    
    F += Fi * Fi;
  }
  F = sqrt(F);

  //    cout << " total force " << F << endl;
  //    cout << endl;
}


void writeForceZ(FILE *f) {

  pair<float,float> out;
  out = doStats(statForce);

  float fz,z,zdev;
  z = Particle->r[2];
  fz = out.first;
  zdev = out.second;
  //  cout << z << " " << fz << " " << zdev << endl;

  fprintf(f,"%f %f %f \n",-z,fz,zdev);
  statForce.clear();


  // pass to Energy integrator
  pair<float,float> in(z,fz);
  forceTime.push_back(in);

}


void writeForce3D() {
  
  int S = forceST.size();
  int T = forceST[0].size();
  //  cout << S << " " << T << endl;

  double F[3];
  double fx, fy, fz;
  for (int t=0; t<T; t++) {
    for (int j=0; j<3; j++) {
      F[j] = forceST[j].at(t);
    }
    fprintf(f5,"%f %f %f\n", F[0], F[1], F[2]);
  }
  
  // print time averaged force
  double Ft[3];
  for (int j=0; j<3; j++) {
    Ft[j] = 0;

    for (int t=0; t<T; t++) {
      Ft[j] += forceST[j].at(t);
    }
    Ft[j] = Ft[j] / T;
    //    cout << Ft[j] << endl;
  }

}

void integrateWrappingEnergy(){
  
  cout << "  integrating Wrapping Energy!" << endl;
  int N = forceTime.size();

  float E=0;
  float z,Fz;
  float dz = _dz;
  for (int j=0; j<N; j++) {
    
    z  = forceTime[j].first;
    Fz = forceTime[j].second;

    E += Fz * dz;

    //write to file E(z)
    fprintf(f7,"%.2f %f\n",-z,E);
  }

  
}

/* gets Deviation for data in HardCoded vector */
float getDev(vector<float> &v) {

  int N = v.size();

  //calc mean
  float sum, mean, dev;
  sum = 0;
  for (int j=0; j < N; j++) {
    sum += v.at(j);
  }
  mean = sum / N;

  //calc stdev
  sum = 0;
  for (int j=0; j < N; j++) {
    sum += (mean - v.at(j))
          *(mean - v.at(j));
  }
  dev = sqrt(sum/N);

  return dev;
}

pair<float,float> doStats(vector<float> &v) {
  
  pair<float,float> out;
  int N = v.size();
  /*  int buffer = 10; //
  float percent = buffer / (float)N;
  cout << "throwing out first " << percent 
       << "% of force time steps for equlibrium" << endl;
  */
  //calc mean
  float sum, mean, dev;
  sum = 0;
  for (int j=0; j < N; j++) {
    //if (j < buffer) continue;
    
    sum += v.at(j);
  }
  mean = sum / N;
  //  mean = sum / (N - buffer);

  //calc stdev
  sum = 0;
  for (int j=0; j < N; j++) {
    sum += (mean - v.at(j))
          *(mean - v.at(j));
  }
  dev = sqrt(sum/N);

  out.first = mean;
  out.second = dev;
  return out;
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
