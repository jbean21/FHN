# include <cstdlib>
# include <iostream>
# include <iomanip>
# include <cmath>
#include <fstream>
#include <string>
#include <stdlib.h>
#include <ctime>
#include <time.h>
#include <omp.h>
#include <vector>
#include <algorithm>

using namespace std;

# include "ode.hpp"

int main();
void test01();
void f01(double t, double y[], double yp[]);
double coupling(int i, int neqn, double y[]);
//double thresh(int i, int neqn, double y[]);
double a2abidir(int i, int neqn, double y[]);
double nneigh(int i, int neqn, double y[]);
int** matrix_creation(int dimN, int neighbours);
//double I_t(double thresh);
vector<int> locations(int locy);
vector<double> PulseGenerator(double thresh);

//Global Variables
int N=9;
int dimN = (int)sqrt(N);
int neighbours;
int neqn = 2 * N;
double k = 0.113;
double I;
double thresh;
int locy;
vector<double> pulse_series;

//****************************************************************************80

int main()

//****************************************************************************80
//
//  Purpose:
//
//    MAIN is the main program for ODE_PRB.
//
//  Discussion:
//
//    ODE_PRB tests the ODE library.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    15 January 2012
//
//  Author:
//
//    John Burkardt
//
{
  
  
  
  cout << "\n";
  cout << "ODE_PRB\n";
  cout << "  C++ version\n";
  cout << "  Test the ODE library.\n";
  
  
  
  
  //Simulation Parameters
  //double k_min = 0;
  //double k_max = 1.0;
  int locmin = 2;
  int locmax = N-1;
  double thresh_min = 0;
  int neighmax;
  double t_start, t_end;
  int j;
  double I_min = -0.58;
  double I_max = -1.4;
  
  /*-------------*/
  /* Start timer */
  /*-------------*/
#ifndef _OPENMP
  /* Use clock_gettime */
  clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &ti);
  t_start = (double)ti.tv_sec + ti.tv_nsec / 1.0E9;
#endif
#ifdef _OPENMP
  /* Use OpenMP clock */
  t_start = omp_get_wtime();
#endif
  
  
  while (N <= 225) {
    if (N % 2 != 0) { neighmax = (N - 1) / 2; }
    else { neighmax = (N - 2) / 2; }
    for (neighbours = 1; neighbours <= neighmax; neighbours++) {
      for (locy = locmin; locy <= locmax; locy++) {
        //for (k = k_min; k <= k_max; k += 0.05) {
          thresh = thresh_min;
          //I = I_min;
#pragma omp parallel for default(shared) private(j) schedule(dynamic) num_threads(omp_get_max_threads())
          for (j = 0; j < 20; j++) {
            //pulse_series = PulseGenerator(thresh);
            
            test01();
            //cout << I << endl;
            thresh += 0.05;
            //I -= 0.05;
          }
        }
      //}
    }
    N = (int)(sqrt(N) + 1)*(sqrt(N) + 1);
    dimN = (int)sqrt(N);
    neqn = 2 * N;
  }
  
  /*------------*/
  /* Stop timer */
  /*------------*/
#ifndef _OPENMP
  /* Use clock_gettime */
  clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &tf);
  t_end = (double)tf.tv_sec + tf.tv_nsec / 1.0E9;
#endif
#ifdef _OPENMP
  /* Use OpenMP clock */
  t_end = omp_get_wtime();
#endif
  
  
  //
  //  Terminate.
  //
  cout << "\n";
  cout << "ODE_PRB\n";
  cout << "  Normal end of execution.\n";
  cout << "\n";
  
  
  return 0;
}
//****************************************************************************80

void test01()

//****************************************************************************80
//
//  Purpose:
//
//    TEST01 tests ODE.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    15 November 2012
//
//  Author:
//
//    John Burkardt
//
{
  double abserr;
  int i;
  int iflag;
  int iwork[5];
  double pi = 3.141592653589793;
  double relerr;
  double step_size = 1.0e-5;
  double step_num = 1 / step_size;
  double t;
  double tf = 150; //
  double tout; //
  double *work;
  double *y;
  
  
  //Create results file
  ofstream resultsFile;
  string filename = to_string(N) + "_" + to_string(neighbours) + "_" + to_string(locy) + "_" + to_string(thresh) + ".csv";
  resultsFile.open(filename);
  resultsFile << "t,";
  for (i = 0; i < neqn / 2; i++) {
    string x;
    string y;
    x = "x" + to_string(i + 1);
    y = "y" + to_string(i + 1);
    resultsFile << x << ",";
  }
  resultsFile << endl;
  
  
  
  abserr = 0.00001;
  relerr = 0.00001;
  
  iflag = 1;
  
  t = 0.0;
  y = new double[neqn];
  
  //INITIAL CONDITIONS
  for (i = 0; i < neqn; i++) { y[i] = 0.0; }
  
  
  /*cout << "  " << setw(8) << t
   << "  " << setw(14) << y[0]
   << "  " << setw(14) << y[1] << "\n";*/
  
  work = new double[100 + 21 * neqn];
  int fuckoff;
  for (i = 1; i <= step_num; i++)
  {
    if ((rand() % 100)*0.01 > thresh) { I = -1.0; }
    else { I = 0.0; }
    
    tout = (double)(i) * tf * step_size;
    
    ode(f01, neqn, y, t, tout, relerr, abserr, iflag, work, iwork);
    
    if (iflag != 2)
    {
      cout << "\n";
      cout << "TEST01 - Fatal error!\n";
      cout << "  ODE returned IFLAG = " << iflag << "\n";
      break;
    }
    /*cout << "  " << setw(8) << t
     << "  " << setw(14) << y[0]
     << "  " << setw(14) << y[1] << "\n";*/
    resultsFile << t << ",";
    for (fuckoff = 0; fuckoff < neqn - 1; fuckoff += 2) {
      resultsFile << y[fuckoff] << ",";
    }
    resultsFile << endl;
  }
  
  delete[] work;
  delete[] y;
  
  resultsFile.close();
  
  return;
}
//****************************************************************************80


//****************************************************************************80

void f01(double t, double y[], double yp[])

//****************************************************************************80
//
//  Purpose:
//
//    F01 supplies the right hand side of the ODE for problem 1.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    15 January 2012
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double T, the time.
//
//    Input, double Y[], the dependent variable.
//
//    Output, double YP[], the value of the derivative.
//
{
  
  int i;
  
  double a = 0.7;
  double b = 0.8;
  double c = 3.0;
  
  for (i = 0; i < neqn-1; i += 2) {
    yp[i] = c * (y[i] - 0.3333333333333*(y[i] * y[i] * y[i]) + y[i+1] + nneigh(i,neqn,y));
  }
  for (i = 1; i < neqn; i += 2) {
    yp[i] = -(y[i-1] - a + b * y[i]) / c;
  }
  return;
}
double coupling(int i, int neqn, double y[]) {
  double I = -0.58;
  double k = 0.113;
  if (i == 0) { return I +k*(y[i+2]-y[i]); }//+ k*(y[i+2] - 2*y[i] + y[neqn-2])
  else if (i == neqn - 2) { return k * ( - y[i] + y[i - 2]); }
  else { return k * (y[i-2] - 2*y[i] + y[i + 2]); }
}

//double thresh(int i, int neqn, double y[]) {
//  srand(time(NULL));
//  double xthresh = (rand() % 20) / 10;
//  double I = -0.58;
//  double II = -0.58;
//  if (i == 0) { return I; }
//  else {
//    if (y[i - 2] > xthresh) {
//      return II;
//    }
//    else { return 0; }
//  }
//}

double a2abidir(int i, int neqn, double y[]) {
  int j;
  double I = -0.58;
  double K = 0;
  double k;
  if (i != 0 ) {
    for (j = 0; j < neqn; j += 2) {
      if (j != i) {
        k = 0.01*(rand() % 100);
        K += k * (y[j] - y[i]);
        //cout << k << endl;
      }
      else { K += 0; }
    }
  }
  else { K = I; } //This means that the first neuron isn't affected by the activity of the rest.
  return K;
}

double nneigh(int i, int neqn, double y[]) {
  //Create matrix of connections
  int** Y = matrix_creation(dimN, neighbours);
  //create driving locations array
  vector<int> locs = locations(locy);
  int j;
  double K = 0;
  if (find(locs.begin(), locs.end(), i) != locs.end() == false) { //(find(locs.begin(), locs.end(), i) != locs.end())!=false
    //if statements tests whether i is in the driving locations vector, if false then execute following
    for (j = 0; j < dimN; j++) {
      if (Y[(i/2)%dimN][j] == 1) {
        K += k * (y[2 * j] - y[i]); //bidirectional coupling
      }
      else { K += 0; }
    }
  }
  else { K = I;} //edit to couple driven neuron,
  //delete connection matrix
  int deleteint;
  for (deleteint = 0; deleteint < dimN; deleteint++) {
    delete[] Y[deleteint];
  }
  delete[] Y;
  //return coupling
  return K;
}
int** matrix_creation(int dimN, int neighbours) {
  int** Y = new int*[dimN];
  for (int i = 0; i < dimN; ++i)
    Y[i] = new int[dimN];
  int rowM;
  int colM;
  int neighcounter;
  for (rowM = 0; rowM < dimN; rowM++) {
    for (colM = 0; colM < dimN; colM++) {
      Y[rowM][colM] = 0;
    }
  }
  for (rowM = 0; rowM < dimN; rowM++) {
    for (neighcounter = -neighbours; neighcounter < neighbours + 1; neighcounter++) {
      if (neighcounter != 0) {
        if (rowM + neighcounter < 0) { Y[rowM][dimN + rowM + neighcounter] = 1; }
        else if (rowM + neighcounter > dimN - 1) { Y[rowM][(rowM + neighcounter) % dimN] = 1; }
        else { Y[rowM][rowM + neighcounter] = 1; }
      }
      else { Y[rowM][rowM + neighcounter] = 0; }
    }
    
  }
  return Y;
}

//double I_t(double thresh) {
//  //Create a random rectangular pulse function
//  //To do this, generate a random number between 0 and 1
//  //If number>thresh, return pulse, else return 0.
//  double number = (rand() % 100) / 100;
//  if (number > thresh) { return -1.0; }
//  else { return 0; }
//}

vector<int> locations(int locy) {
  //Creates a vector containing the eqn numbers which will be driven
  vector<int> locs = {};
  
  int i;
  int r = 2*(rand()%N);
  
  for (i=0; i<locy; i++) {
    if (find(locs.begin(), locs.end(), r) != locs.end() == false) {
      locs.push_back(r);
    }
    else {
      //locs.push_back((r+dimN)%N);
    }
  }
  
  return locs;
}

vector<double> PulseGenerator(double thresh) {
  //Generate random seed
  srand(time(NULL));
  //Create vector to store random numbers in
  vector<double> pulse_series = {};
  int i;
  for (i = 0; i < 1e5; i++) { //1e5 = step num from f01, if it changes there it must also change here
    if ((rand() % 100) / 100 > thresh) { pulse_series.push_back(-1.0); } //if random number generated is greate than thresh pushback -1
    else { pulse_series.push_back(0); }
  }
  return pulse_series;
}

