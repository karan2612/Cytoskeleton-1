//#include "global.h"

#include <fstream>
#include <iostream>

using namespace std;

void process(char* in, std::ifstream &fin) {

  if (strcmp(in, "SysSize") == 0) {

    cout << "reading: " << in << endl;

    int x;
    fin >> x; //gets data
 
    cout << _nSys;
    _nSys = x;
    cout << " fate changed! " << _nSys << endl;
  }

  if (strcmp(in, "EQSamp") == 0) {

    cout << "reading: " << in << endl;

    int x;
    fin >> x; //gets data
 
    cout << _tEqSamp;
    _tEqSamp = x;
    cout << " fate changed! " << _tEqSamp << endl;
  }

  if (strcmp(in, "Radius") == 0) {

    cout << "reading: " << in << endl;

    float r;
    fin >> r; //gets data
    cout << _pRadius;
    _pRadius = r * _lActin;
    cout << " fate changed! " << _pRadius << endl;

  }

  if (strcmp(in, "Zara") == 0) {
    cout << "true love" << endl;
  }

}


void readInput() {

  cout << "Reading from file" << endl;   

  char data[100];
  ifstream infile; 
  /*
  if (!fname) {
    cout << "loading default" << endl;
    infile.open(fname);
  }
  */
  infile.open("fooIN.txt");


  while (infile >> data) {

    process(data, infile);

  }

  infile.close();
}
