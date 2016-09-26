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

  if (strcmp(in, "tSamp") == 0) {

    cout << "reading: " << in << endl;

    int x;
    fin >> x; //gets data
 
    cout << _tSamp;
    _tSamp = x;
    cout << " fate changed! " << _tSamp << endl;
  }

  if (strcmp(in, "nStep") == 0) {

    cout << "reading: " << in << endl;

    int x;
    fin >> x; //gets data
 
    cout << _nSteps;
    _nSteps = x;
    cout << " fate changed! " << _nSteps << endl;
  }

  if (strcmp(in, "Radius") == 0) {

    cout << "reading: " << in << endl;

    float r;
    fin >> r; //gets data
    cout << _pRadius;
    _pRadius = r * _lActin;
    cout << " fate changed! " << _pRadius << endl;
  }

  if (strcmp(in, "Print") == 0) {

    cout << "reading: " << in << endl;

    bool p;
    fin >> p;
    _printTime = p;
    cout << " fate changed! " << _printTime << endl;
  }

  if (strcmp(in, "Tag") == 0) {

    cout << "reading: " << in << endl;
   
    fin >> _tag;
    _fileTag = true;

    cout << " fate changed! " << _tag << endl;
  }

  if (strcmp(in, "Zara") == 0) {
    cout << "true love" << endl;
  }

}


void readInput() {

  cout << "Reading from file" << endl;   

  char data[100];
  ifstream infile; 

  infile.open("command.txt");
  while (infile >> data) {

    process(data, infile);
  }
  infile.close();

  /* recalc time */
  _tMax = _tEqGlobal + (_tSamp * _nSteps);
}
