/* This file includes
   process()
   checktime()
   readInput()
   we suggest the reader start from the bottom

   printLog()
 */

#include <fstream>
#include <iostream>
#include <ctime>

using namespace std;

void process(char* in, std::ifstream &fin) {

  if (strcmp(in, "SysSize") == 0) {

    cout << "reading: " << in << endl;
    int x; fin >> x; //gets data
 
    cout << "   " << _nSys << endl;
    _nSys = x;
    cout << "   " << _nSys << endl;
  }

  if (strcmp(in, "dt") == 0) {

    cout << "reading: " << in << endl;
    float x; fin >> x; //gets data
 
    cout << "   " << _dt << endl;
    _dt = x;
    cout << "   " << _dt << endl;
  }


  if (strcmp(in, "tSamp") == 0) {

    cout << "reading: " << in << endl;
    int x; fin >> x; //gets data
 
    cout << "   " << _tSamp << endl;
    _tSamp = x;
    cout << "   " << _tSamp << endl;
  }

  if (strcmp(in, "tEqLocal") == 0) {

    cout << "reading: " << in << endl;
    int x; fin >> x; //gets data
 
    cout << "   " << _tEqLocal << endl;
    _tEqLocal = x;
    cout << "   " << _tEqLocal << endl;
  }

  if (strcmp(in, "tEqGlobal") == 0) {

    cout << "reading: " << in << endl;
    int x; fin >> x; //gets data
 
    cout << "   " << _tEqGlobal << endl;
    _tEqGlobal = x;
    cout << "   " << _tEqGlobal << endl;
  }

  if (strcmp(in, "nStep") == 0) {

    cout << "reading: " << in << endl;
    int x; fin >> x; //gets data
 
    cout << "   " << _nSteps << endl;
    _nSteps = x;
    cout << "   " << _nSteps << endl;
  }

  if (strcmp(in, "Radius") == 0) {

    cout << "reading: " << in << endl;
    float r; fin >> r; //gets data

    cout << "   " << _pRadius << endl;
    _pRadius = r * _lActin;
    cout << "   " << _pRadius << endl;
  }

  if (strcmp(in, "nSpectrin") == 0) {

    cout << "reading: " << in << endl;
    int n; fin >> n; //gets data

    cout << "   " << _nSpectrin << endl;
    _nSpectrin = n;
    cout << "   " << _nSpectrin << endl;
  }

  if (strcmp(in, "SpectrinContour") == 0) {

    cout << "reading: " << in << endl;
    float c; fin >> c; //gets data

    cout << "   " << _Contour << endl;
    _Contour = c * _lActin;
    cout << "   " << _Contour << endl;
  }

  if (strcmp(in, "Tag") == 0) {

    cout << "reading: " << in << endl;
   
    _fileTag = true;
    fin >> _tag;
    cout << "   " << _tag << endl;
  }

  if (strcmp(in, "Print") == 0) {

    cout << "reading: " << in << endl;
    bool p; fin >> p;

    cout << "   " << _printTime << endl;
    _printTime = p;
    cout << "   " << _printTime << endl;
  }

  if (strcmp(in, "Sunrise") == 0) {

    cout << "reading: " << in << endl;
    bool p; fin >> p;

    cout << "   " << _sunrise << endl;
    _sunrise = p;
    cout << "   " << _sunrise << endl;
  }

  if (strcmp(in, "Zara") == 0) {
    cout << "true love" << endl;
  }

}


/* Timing Considerations and Recalculation */
void checkTime() {

  if (_tSamp <= _tEqLocal) {

    cout << "\nFatal Error: Equilibrium Time (local) exceeds sample Time" << endl;
    cout << "   local Equlibrium delay: " << _tEqLocal << endl;
    cout << "   allowed sample time: " << _tSamp << endl;
    cout << "** measurements will never be made!!" << endl;
    getc(stdin);
  }

  //recalc max time
  _tMax = _tEqGlobal + (_tSamp * _nSteps);
}

/* this is the main() file */

void readInput() {

  cout << "Reading from file" << endl;   

  char data[100];
  ifstream infile; 
  infile.open("command.txt");

  while (infile >> data) {

    process(data, infile);
  }
  infile.close();

  checkTime();
}


/* append */
void printLog() {

  cout << "Preparing .log file" << endl;
  char fname[64] = "";

  if (_fileTag) strcat(fname, _tag.c_str() );
  else strcat(fname, "temp");

  strcat(fname, ".log");
  cout << fname << endl;

  std::ofstream log (fname, std::ofstream::out);

  /* time */
  log << "Time Parameters" << endl;
  log << " nSteps: " << _nSteps << endl;
  log << " tSamp: " << _tSamp << endl;
  log << " tMax: " << _tMax << endl;
  log << " dt: " << _dt << endl;
  log << " tEqGlobal: " << _tEqGlobal << endl;
  log << " tEqLocal: " << _tEqLocal << endl;


  /* system */
  log << "System Parameters" << endl;
  log << " nSys: " << _nSys << endl;
  log << " lActin: " << _lActin << endl;
  log << " Contour: " << _Contour << endl;
  log << " nSpectrin: " << _nSpectrin << endl;
  log << " spectrinRadius: " << _sRadius << endl;
  log << " particleRadius: " << _pRadius << endl;
  log << " sunrise: " << _sunrise << endl;

  /* physics */
  log << "Physics Parameters" << endl;
  log << " D: " << _D << endl;
  log << " gamma: " << _gamma << endl;
  log << " epsilon: " << _epsilon << endl;
  log << " k: " << _k << endl;
  log << " m: " << _m << endl;
  log << " dz: " << _dz << endl;

  /* finish */
  time_t now = time(0);
  char* cstamp = ctime(&now);
  log << "Launch time (moscow): " << cstamp << endl;
  log.close();

}
