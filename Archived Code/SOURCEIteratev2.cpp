# include <cstdlib>
# include <iostream>
# include <iomanip>
# include <cmath>
# include <fstream>
# include <string>
# include <stdlib.h>
# include <ctime>
# include <time.h>
# include <omp.h>
# include <vector>
# include <algorithm>
# include "ode.hpp"

using namespace std;

// Global Variables
int N=4;
int sqrtN = (int)sqrt(N);
int neighbours;
int neqn = 2 * N;
int cbrtN = cbrt(N);
int cbrtNsq = cbrtN*cbrtN;
double sqrt2 = sqrt(2);
double k, thresh;
int locy;
int timeindex;

double step_size = 1.0e-5;
double step_num = 1 / step_size;

// Make thresh threadprivate
#pragma omp threadprivate(thresh)

// Matrix of values
double **Y = new double*[N];

string vname;
string command;

// System micro parameters
double a = 0.75;
double b = 0.8;
double c = 3.0;
double c_inv = 1/c;

// Change to 1.5 for driving from middle of 2D lattice
// 3D: I used I=2
double I_ext = 0.58;

// Declare an array on the stack
// (hence do not need to free)
int main();
void test01();
double Iext_Thresh(int i);
void f01(double t, double y[], double yp[]);

//****************************************************************************80

int main() {
  // Simulation Parameters
  int Nmin = 4;
  int Nmax = 64;
  double N_inv;
  
  double k_min = 0;
  double k_max = 1.0;
  
  int locmin = 1;
  int locmax = N-1;
  
  double thresh_min = 0.0;
  
  int neighmax;
  
  double I_min = -0.58;
  double I_max = -1.4;
  
  int j;
  int jmax = 20;
  int jmin = 0;
  
  // Timing variables
  double t_start, t_end;
  
  // Create file of computation times
  ofstream timesFile;
  string timesfilename = "Times.csv";
  timesFile.open(timesfilename);
  
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
  
  N = Nmin;
  neqn = 2*N;
  
  // Vary number of neurones up to 10000
  while (N <= Nmax) {
    
    for (int i=0; i<N; i++){ Y[i] = new double[2]; }
    
    int neighcounter;
    
    if (N % 2 != 0) { neighmax = (N - 1) / 2; }
    else            { neighmax = (N - 2) / 2; }
    
    // Vary number of neighbours used in coupling
    for (neighbours = 1; neighbours <= neighmax; neighbours++) {
      
      //Create matrix of connections
      //M = matrix_creation(sqrtN, neighbours);
      
      // Vary density of drivers
      for (locy = locmin; locy <= locmax; locy++) {
        // Vary the coupling strength
        for (k = k_min; k <= k_max; k += 0.05) {
          thresh = thresh_min;
          //I = I_min;
#pragma omp parallel for default(shared) private(j, vname) schedule(dynamic) num_threads(omp_get_max_threads())
          // Vary threshold current
          for (j = jmin; j < jmax; j++) {
            thresh = j*0.05;
            //printf("iteration %d on thread %d : thresh = %f\n", j, omp_get_thread_num(), thresh);
            //pulse_series = PulseGenerator(thresh);
            
            test01();
          }
          
          //system("python PLOTIterate.py");
          //command = "rm -f *.csv";
          //system(command.c_str());
          
        }
      }
    }
    
    // This program is for 2D only
    sqrtN += 1;
    N = sqrtN*sqrtN;
    neqn = 2 * N;
    
    /* Use OpenMP clock */
    t_end = omp_get_wtime();
    
    // Output time taken for current N to CSV
    timesFile << "," << N << "," << (t_end-t_start) << ",";
    timesFile << endl;
    
    for (int j=0; j<N; j++) {
      
    }
    
  } // N
  
  timesFile.close();
  
  return 0;
	return 0;
}
//****************************************************************************80

void test01() {
	double abserr;
	int i;
	int iflag;
	int iwork[5];
	double pi = 3.141592653589793;
	double relerr;
	double tf = 100; //important
	double stepSize = 1e-5; //important
	double step_num = 1 / stepSize; //unsure about this
	double t;
	double tout;
	double *work;
	double *y;

	/*cout << "      T              ";
	for (int i = 1; i < N + 1; i++) {
		cout << "v";
		cout << i;
		cout << "              ";
		cout << "w";
		cout << i;
		cout << "              ";
	}
	cout << "\n\n";
  */
	abserr = 0.00001;
	relerr = 0.00001;

	iflag = 1;

	t = 0.0;
	y = new double[neqn];

	// Initial Conditions
	for (i = 0; i < neqn; i++) { y[i] = 0.0; }
	
	/* -------------------------------------------- */
  	//////////////////////////
	/* Initialise CVS Files */
	//////////////////////////

	//Create CSV file for results to be written to
	ofstream resultsFile;
	resultsFile.open("FHN.csv"); //change filename to something more descriptive


	/* -------------------------------------------- */
  	//////////////////////////////
	/* Write initial conditions */
	//////////////////////////////
	//cout << "  " << setw(8) << t;
	//resultsFile << t;
	//for (int i = 0; i < neqn; i++) {
	//	cout << "  " << setw(14) << y[i];
	//	resultsFile << ",";
	//	resultsFile << y[i];
	//}
	//cout << "\n";
	//resultsFile << endl;


	// Output 2D matrix
	//cout << "\n";
	
  resultsFile << t;
	for(i=0; i<N;i++) {
		for(int j=0; j<2; j++) {
			//cout << setw(5) << "  " << Y[i][j];
			resultsFile << ",";
			resultsFile << Y[i][j];
		}
	//cout << "\n";
	}
	resultsFile << endl;
  

  
	/* -------------------------------------------- */

	work = new double[100 + 21 * neqn]; //want to work out what this does and why it uses 100 and 21
  	
	/* -------------------------------------------- */
	//////////////////////
	/* All other values */
	//////////////////////
	for (timeindex = 1; timeindex <= step_num; timeindex++)
	{	
		tout = (double)(timeindex)* tf * stepSize;

		ode(f01, neqn, y, t, tout, relerr, abserr, iflag, work, iwork);

		if (iflag != 2)
		{
			cout << "\n";
			cout << "TEST01 - Fatal error!\n";
			cout << "  ODE returned IFLAG = " << iflag << "\n";
			break;
		}

		// Print the output to console and to csv
		//cout << "  " << setw(8) << t;
		//resultsFile << t;
		//for (int i = 0; i < neqn; i++) {
		//	cout << "  " << setw(14) << y[i];
		//	resultsFile << ",";
		//	resultsFile << y[i];
		//}


		// Output 2D matrix
		//cout << "\n";
#pragma omp critical
    {
    for(i=0; i<N; i++) {
      for(int j=0; j<2; j++) {
        //cout << setw(5) << "  " << Y[i][j];
        resultsFile << ",";
        resultsFile << Y[i][j];
      }
    }
    //cout << "\n";
    resultsFile << endl;
    }
	}
	/* -------------------------------------------- */
	
	delete[] work;
	delete[] y;

	return;
}

//****************************************************************************80



double bi_nearest_neighbours_3D_non_periodic(double k, double y[], int i) {
	i = i/2; // Convert i from "equation number" to "neuron number"

  int n, j;

  // Define the locations of the corners
  int c0, c1, c2, c3, c4, c5, c6, c7;
  c0 = 0;
  c1 = cbrtN-1;
  c2 = cbrtNsq-cbrtN;
  c3 = cbrtNsq-1;
  c4 = N-1-cbrtNsq+1;
  c5 = N-1-cbrtNsq+cbrtN;
  c6 = N-1-cbrtN;
  c7 = N-1;

	// Coupling is currently not wrapped (non-periodic)

	// Cube corners (3 coupling terms)
	if	    (i == c0) { return k*(Y[1  ][0] - Y[i][0]) + k*(Y[  cbrtN][0] - Y[i][0]) + k*(Y[  cbrtNsq][0] - Y[i][0]); } //0
	else if (i == c1) { return         k*(Y[i-1][0] - Y[i][0]) + k*(Y[i+cbrtN][0] - Y[i][0]) + k*(Y[i+cbrtNsq][0] - Y[i][0]); } //2
  else if (i == c2) { return         k*(Y[i+1][0] - Y[i][0]) + k*(Y[i-cbrtN][0] - Y[i][0]) + k*(Y[i+cbrtNsq][0] - Y[i][0]); } //6
  else if (i == c3) { return         k*(Y[i-1][0] - Y[i][0]) + k*(Y[i-cbrtN][0] - Y[i][0]) + k*(Y[i+cbrtNsq][0] - Y[i][0]); } //8
  else if (i == c4) { return         k*(Y[i+1][0] - Y[i][0]) + k*(Y[i+cbrtN][0] - Y[i][0]) + k*(Y[i-cbrtNsq][0] - Y[i][0]); } //18
  else if (i == c5) { return         k*(Y[i-1][0] - Y[i][0]) + k*(Y[i+cbrtN][0] - Y[i][0]) + k*(Y[i-cbrtNsq][0] - Y[i][0]); } //20
  else if (i == c6) { return         k*(Y[i+1][0] - Y[i][0]) + k*(Y[i-cbrtN][0] - Y[i][0]) + k*(Y[i-cbrtNsq][0] - Y[i][0]); } //24
  else if (i == c7) { return         k*(Y[i-1][0] - Y[i][0]) + k*(Y[i-cbrtN][0] - Y[i][0]) + k*(Y[i-cbrtNsq][0] - Y[i][0]); } //26
  
	// Cube edges (4 coupling terms)
	else if (i > 0            && i < c1) { return k*(Y[i+1][0] - Y[i][0]) + k*(Y[i-1    ][0] - Y[i][0]) + k*(Y[i+cbrtN][0] - Y[i][0]) + k*(Y[i+cbrtNsq][0] - Y[i][0]); }   // Top layer. back edge
  else if (i > c2           && i < c3) { return k*(Y[i+1][0] - Y[i][0]) + k*(Y[i-1    ][0] - Y[i][0]) + k*(Y[i-cbrtN][0] - Y[i][0]) + k*(Y[i+cbrtNsq][0] - Y[i][0]); }   // Top layer, front edge
  else if (i%cbrtN == 0     && i < c2) { return k*(Y[i+1][0] - Y[i][0]) + k*(Y[i-cbrtN][0] - Y[i][0]) + k*(Y[i+cbrtN][0] - Y[i][0]) + k*(Y[i+cbrtNsq][0] - Y[i][0]); }   // Top layer, left edge
  else if ((i+1)%cbrtN == 0 && i < c3) { return k*(Y[i-1][0] - Y[i][0]) + k*(Y[i-cbrtN][0] - Y[i][0]) + k*(Y[i+cbrtN][0] - Y[i][0]) + k*(Y[i+cbrtNsq][0] - Y[i][0]); }   // Top layer, right edge

  else if (i > c4           && i < c5          ) { return k*(Y[i+1][0] - Y[i][0]) + k*(Y[i-1    ][0] - Y[i][0]) + k*(Y[i+cbrtN][0] - Y[i][0]) + k*(Y[i-cbrtNsq][0] - Y[i][0]); }   // Bottom layer. back edge (18, 20)
  else if (i > c6           && i < c7          ) { return k*(Y[i+1][0] - Y[i][0]) + k*(Y[i-1    ][0] - Y[i][0]) + k*(Y[i-cbrtN][0] - Y[i][0]) + k*(Y[i-cbrtNsq][0] - Y[i][0]); }   // Bottom layer, front edge (24, 26)
  else if (i%cbrtN == 0     && i > c4 && i < c6) { return k*(Y[i+1][0] - Y[i][0]) + k*(Y[i-cbrtN][0] - Y[i][0]) + k*(Y[i+cbrtN][0] - Y[i][0]) + k*(Y[i-cbrtNsq][0] - Y[i][0]); }   // Bottom layer, left edge
  else if ((i+1)%cbrtN == 0 && i > c5 && i < c7) { return k*(Y[i-1][0] - Y[i][0]) + k*(Y[i-cbrtN][0] - Y[i][0]) + k*(Y[i+cbrtN][0] - Y[i][0]) + k*(Y[i-cbrtNsq][0] - Y[i][0]); }   // Bottom layer, right edge
  
  /* ---------------------------- */

  // Cube face points (5 connections) -- cube edges caught by upper if statements
  else if (i>cbrtN && i<cbrtNsq-cbrtN)       { return k*(Y[i-1][0] - Y[i][0]) + k*(Y[i+1][0] - Y[i][0]) + k*(Y[i-cbrtN][0] - Y[i][0]) + k*(Y[i+cbrtN][0] - Y[i][0]) + k*(Y[i+cbrtNsq][0] - Y[i][0]); }   // Top layer, face points
  else if (i>N-cbrtNsq+cbrtN && i<N-cbrtN-1) { return k*(Y[i-1][0] - Y[i][0]) + k*(Y[i+1][0] - Y[i][0]) + k*(Y[i-cbrtN][0] - Y[i][0]) + k*(Y[i+cbrtN][0] - Y[i][0]) + k*(Y[i-cbrtNsq][0] - Y[i][0]); }   // Bottom layer, face points


  // Middle layers
  else {
      int cn0, cn1, cn2, cn3;

      // Loop over interior layers
      for (j=1; j<cbrtN; j++) {
          // Find which layer the neurone is in
          if (i < (j+1)*cbrtNsq && i > j*cbrtNsq-1) {
              // Neurone's layer
              n = j;

              // Define the corners of the layer
              // (which are edges of the cube)
              cn0 = n*cbrtNsq;
              cn1 = n*cbrtNsq+cbrtN-1;
              cn2 = (n+1)*cbrtNsq-cbrtN;
              cn3 = (n+1)*cbrtNsq-1;

              break;
              // check that this break only exits this loop!!!
          }
      }
      
      // Use this print statement to see which neurones each
      // if statement evaluates
      // printf("i=%d in layer %d : cn0=%d, cn1=%d, cn2=%d, cn3=%d\n", i, n, cn0, cn1, cn2, cn3);

      // Middle layer corners / Cube edges (4 connections)
      if      (i == cn0) { return k*(Y[i+1][0] - Y[i][0]) + k*(Y[i+cbrtN][0] - Y[i][0]) + k*(Y[i-cbrtNsq][0] - Y[i][0]) + k*(Y[i+cbrtNsq][0] - Y[i][0]); }   // Middle layers. back edge
      else if (i == cn1) { return k*(Y[i-1][0] - Y[i][0]) + k*(Y[i+cbrtN][0] - Y[i][0]) + k*(Y[i-cbrtNsq][0] - Y[i][0]) + k*(Y[i+cbrtNsq][0] - Y[i][0]); }   // Middle layers, front edge
      else if (i == cn2) { return k*(Y[i+1][0] - Y[i][0]) + k*(Y[i-cbrtN][0] - Y[i][0]) + k*(Y[i-cbrtNsq][0] - Y[i][0]) + k*(Y[i+cbrtNsq][0] - Y[i][0]); }   // Middle layers, left edge
      else if (i == cn3) { return k*(Y[i-1][0] - Y[i][0]) + k*(Y[i-cbrtN][0] - Y[i][0]) + k*(Y[i-cbrtNsq][0] - Y[i][0]) + k*(Y[i+cbrtNsq][0] - Y[i][0]); }   // Middle layers, right edge

      // Middle layer sides / Cube face point (5 connections)
      else if (i > cn0 && i < cn1                    ) { return k*(Y[i+1][0] - Y[i][0]) + k*(Y[i-1    ][0] - Y[i][0]) + k*(Y[i+cbrtN][0] - Y[i][0]) + k*(Y[i+cbrtNsq][0] - Y[i][0]) + k*(Y[i-cbrtNsq][0] - Y[i][0]); }   // Middle layers. back edge
      else if (i > cn2 && i < cn3                    ) { return k*(Y[i+1][0] - Y[i][0]) + k*(Y[i-1    ][0] - Y[i][0]) + k*(Y[i-cbrtN][0] - Y[i][0]) + k*(Y[i+cbrtNsq][0] - Y[i][0]) + k*(Y[i-cbrtNsq][0] - Y[i][0]); }   // Middle layers, front edge
      else if (i%cbrtN == 0     && i > cn0 && i < cn2) { return k*(Y[i+1][0] - Y[i][0]) + k*(Y[i-cbrtN][0] - Y[i][0]) + k*(Y[i+cbrtN][0] - Y[i][0]) + k*(Y[i+cbrtNsq][0] - Y[i][0]) + k*(Y[i-cbrtNsq][0] - Y[i][0]); }   // Middle layers, left edge
      else if ((i+1)%cbrtN == 0 && i > cn1 && i < cn3) { return k*(Y[i-1][0] - Y[i][0]) + k*(Y[i-cbrtN][0] - Y[i][0]) + k*(Y[i+cbrtN][0] - Y[i][0]) + k*(Y[i+cbrtNsq][0] - Y[i][0]) + k*(Y[i-cbrtNsq][0] - Y[i][0]); }   // Middle layers, right edge

      // Interior Neurons (6 coupling terms)
      else {
        //printf("i=%d in layer %d : cn0=%d, cn1=%d, cn2=%d, cn3=%d\n", i, n, cn0, cn1, cn2, cn3);
        return I_ext + k*(Y[i+1][0] - Y[i][0]) + k*(Y[i-1][0] - Y[i][0]) + k*(Y[i-cbrtN][0] - Y[i][0]) + k*(Y[i+cbrtN][0] - Y[i][0]) + k*(Y[i-cbrtNsq][0] - Y[i][0]) + k*(Y[i+cbrtNsq][0] - Y[i][0]); }
  }
}

double bi_nearest_neighbours_3D_periodic(double k, double y[], int i) {
	i = i/2; // Convert i from "equation number" to "neuron number"

  int n, j;
  double I_1;

  // Define the locations of the corners
  int c0, c1, c2, c3, c4, c5, c6, c7;
  c0 = 0;
  c1 = cbrtN-1;
  c2 = cbrtNsq-cbrtN;
  c3 = cbrtNsq-1;
  c4 = N-1-cbrtNsq+1;
  c5 = N-1-cbrtNsq+cbrtN;
  c6 = N-1-cbrtN;
  c7 = N-1;

	// Coupling is currently not wrapped (non-periodic)

	// Cube corners (6 coupling terms)
  if	    (i == 0 ) { return         k*(Y[c1 ][0] - Y[i][0]) + k*(Y[i+1][0] - Y[i][0]) + k*(Y[c2     ][0] - Y[i][0]) + k*(Y[  cbrtN][0] - Y[i][0]) + k*(Y[c4       ][0] - Y[i][0]) + k*(Y[  cbrtNsq][0] - Y[i][0]); } //0
  else if (i == c1) { return         k*(Y[i-1][0] - Y[i][0]) + k*(Y[c0 ][0] - Y[i][0]) + k*(Y[c3     ][0] - Y[i][0]) + k*(Y[i+cbrtN][0] - Y[i][0]) + k*(Y[c5       ][0] - Y[i][0]) + k*(Y[i+cbrtNsq][0] - Y[i][0]); } //2
  else if (i == c2) { return         k*(Y[c3 ][0] - Y[i][0]) + k*(Y[i+1][0] - Y[i][0]) + k*(Y[i-cbrtN][0] - Y[i][0]) + k*(Y[c0     ][0] - Y[i][0]) + k*(Y[c6       ][0] - Y[i][0]) + k*(Y[i+cbrtNsq][0] - Y[i][0]); } //6
  else if (i == c3) { return         k*(Y[i-1][0] - Y[i][0]) + k*(Y[c2 ][0] - Y[i][0]) + k*(Y[i-cbrtN][0] - Y[i][0]) + k*(Y[c1     ][0] - Y[i][0]) + k*(Y[c7       ][0] - Y[i][0]) + k*(Y[i+cbrtNsq][0] - Y[i][0]); } //8
  else if (i == c4) { return         k*(Y[c5 ][0] - Y[i][0]) + k*(Y[i+1][0] - Y[i][0]) + k*(Y[c6     ][0] - Y[i][0]) + k*(Y[i+cbrtN][0] - Y[i][0]) + k*(Y[i-cbrtNsq][0] - Y[i][0]) + k*(Y[c0       ][0] - Y[i][0]); } //18
  else if (i == c5) { return         k*(Y[i-1][0] - Y[i][0]) + k*(Y[c4 ][0] - Y[i][0]) + k*(Y[c7     ][0] - Y[i][0]) + k*(Y[i+cbrtN][0] - Y[i][0]) + k*(Y[i-cbrtNsq][0] - Y[i][0]) + k*(Y[c1       ][0] - Y[i][0]); } //20
  else if (i == c6) { return         k*(Y[c7 ][0] - Y[i][0]) + k*(Y[i+1][0] - Y[i][0]) + k*(Y[i-cbrtN][0] - Y[i][0]) + k*(Y[c4     ][0] - Y[i][0]) + k*(Y[i-cbrtNsq][0] - Y[i][0]) + k*(Y[c2       ][0] - Y[i][0]); } //24
  else if (i == c7) { return         k*(Y[i-1][0] - Y[i][0]) + k*(Y[c6 ][0] - Y[i][0]) + k*(Y[i-cbrtN][0] - Y[i][0]) + k*(Y[c5     ][0] - Y[i][0]) + k*(Y[i-cbrtNsq][0] - Y[i][0]) + k*(Y[c3       ][0] - Y[i][0]); } //26
  
	// Cube edges (6 coupling terms)
	else if (i > 0            && i < c1)           { return k*(Y[i-1      ][0] - Y[i][0]) + k*(Y[i+1      ][0] - Y[i][0]) + k*(Y[i-cbrtN+cbrtNsq][0] - Y[i][0]) + k*(Y[i+cbrtN        ][0] - Y[i][0]) + k*(Y[i+cbrtNsq][0] - Y[i][0]) + k*(Y[i+(cbrtN-1)*cbrtNsq][0] - Y[i][0]); }   // Top layer. back edge
  else if (i > c2           && i < c3)           { return k*(Y[i-1      ][0] - Y[i][0]) + k*(Y[i+1      ][0] - Y[i][0]) + k*(Y[i-cbrtN        ][0] - Y[i][0]) + k*(Y[i+cbrtN-cbrtNsq][0] - Y[i][0]) + k*(Y[i+cbrtNsq][0] - Y[i][0]) + k*(Y[i+(cbrtN-1)*cbrtNsq][0] - Y[i][0]); }   // Top layer, front edge
  else if (i%cbrtN == 0     && i < c2)           { return k*(Y[i-1+cbrtN][0] - Y[i][0]) + k*(Y[i+1      ][0] - Y[i][0]) + k*(Y[i-cbrtN        ][0] - Y[i][0]) + k*(Y[i+cbrtN        ][0] - Y[i][0]) + k*(Y[i+cbrtNsq][0] - Y[i][0]) + k*(Y[i+(cbrtN-1)*cbrtNsq][0] - Y[i][0]); }   // Top layer, left edge
  else if ((i+1)%cbrtN == 0 && i < c3)           { return k*(Y[i-1      ][0] - Y[i][0]) + k*(Y[i-cbrtN+1][0] - Y[i][0]) + k*(Y[i-cbrtN        ][0] - Y[i][0]) + k*(Y[i+cbrtN        ][0] - Y[i][0]) + k*(Y[i+cbrtNsq][0] - Y[i][0]) + k*(Y[i+(cbrtN-1)*cbrtNsq][0] - Y[i][0]); }   // Top layer, right edge

  else if (i > c4           && i < c5          ) { return k*(Y[i-1      ][0] - Y[i][0]) + k*(Y[i+1      ][0] - Y[i][0]) + k*(Y[i-cbrtN-cbrtNsq][0] - Y[i][0]) + k*(Y[i+cbrtN        ][0] - Y[i][0]) + k*(Y[i-cbrtNsq][0] - Y[i][0]) + k*(Y[i-(cbrtN-1)*cbrtNsq][0] - Y[i][0]); }   // Bottom layer. back edge (18, 20)
  else if (i > c6           && i < c7          ) { return k*(Y[i-1      ][0] - Y[i][0]) + k*(Y[i+1      ][0] - Y[i][0]) + k*(Y[i-cbrtN        ][0] - Y[i][0]) + k*(Y[i+cbrtN-cbrtNsq][0] - Y[i][0]) + k*(Y[i-cbrtNsq][0] - Y[i][0]) + k*(Y[i-(cbrtN-1)*cbrtNsq][0] - Y[i][0]); }   // Bottom layer, front edge (24, 26)
  else if (i%cbrtN == 0     && i > c4 && i < c6) { return k*(Y[i-1+cbrtN][0] - Y[i][0]) + k*(Y[i+1      ][0] - Y[i][0]) + k*(Y[i-cbrtN        ][0] - Y[i][0]) + k*(Y[i+cbrtN        ][0] - Y[i][0]) + k*(Y[i-cbrtNsq][0] - Y[i][0]) + k*(Y[i-(cbrtN-1)*cbrtNsq][0] - Y[i][0]); }   // Bottom layer, left edge
  else if ((i+1)%cbrtN == 0 && i > c5 && i < c7) { return k*(Y[i-1      ][0] - Y[i][0]) + k*(Y[i+1-cbrtN][0] - Y[i][0]) + k*(Y[i-cbrtN        ][0] - Y[i][0]) + k*(Y[i+cbrtN        ][0] - Y[i][0]) + k*(Y[i-cbrtNsq][0] - Y[i][0]) + k*(Y[i-(cbrtN-1)*cbrtNsq][0] - Y[i][0]); }   // Bottom layer, right edge
  
  /* ---------------------------- */

  // Cube face points (5 connections) -- cube edges caught by upper if statements
  else if (i>cbrtN && i<cbrtNsq-cbrtN)       { return k*(Y[i-1][0] - Y[i][0]) + k*(Y[i+1][0] - Y[i][0]) + k*(Y[i-cbrtN][0] - Y[i][0]) + k*(Y[i+cbrtN][0] - Y[i][0]) + k*(Y[i+(cbrtN-1)*cbrtNsq][0] - Y[i][0]) + k*(Y[i+cbrtNsq          ][0] - Y[i][0]); }   // Top layer, face points
  else if (i>N-cbrtNsq+cbrtN && i<N-cbrtN-1) { return k*(Y[i-1][0] - Y[i][0]) + k*(Y[i+1][0] - Y[i][0]) + k*(Y[i-cbrtN][0] - Y[i][0]) + k*(Y[i+cbrtN][0] - Y[i][0]) + k*(Y[i-cbrtNsq          ][0] - Y[i][0]) + k*(Y[i-(cbrtN-1)*cbrtNsq][0] - Y[i][0]); }   // Bottom layer, face points


  // Middle layers
  else {
      int cn0, cn1, cn2, cn3;

      // Loop over interior layers
      for (j=1; j<cbrtN; j++) {
          // Find which layer the neurone is in
          if (i < (j+1)*cbrtNsq && i > j*cbrtNsq-1) {
              // Neurone's layer
              n = j;

              // Define the corners of the layer
              // (which are edges of the cube)
              cn0 = n*cbrtNsq;
              cn1 = n*cbrtNsq+cbrtN-1;
              cn2 = (n+1)*cbrtNsq-cbrtN;
              cn3 = (n+1)*cbrtNsq-1;

              break;
              // check that this break only exits this loop!!!
          }
      }
      
      // Use this print statement to see which neurones each
      // if statement evaluates
      // printf("i=%d in layer %d : cn0=%d, cn1=%d, cn2=%d, cn3=%d\n", i, n, cn0, cn1, cn2, cn3);

      // Middle layer corners / Cube edges (6 connections)
      if      (i == cn0) { return k*(Y[cn1][0] - Y[i][0]) + k*(Y[i+1][0] - Y[i][0]) + k*(Y[cn2    ][0] - Y[i][0]) + k*(Y[i+cbrtN][0] - Y[i][0]) + k*(Y[i-cbrtNsq][0] - Y[i][0]) + k*(Y[i+cbrtNsq][0] - Y[i][0]); }   // Middle layers. back edge
      else if (i == cn1) { return k*(Y[i-1][0] - Y[i][0]) + k*(Y[cn0][0] - Y[i][0]) + k*(Y[cn3    ][0] - Y[i][0]) + k*(Y[i+cbrtN][0] - Y[i][0]) + k*(Y[i-cbrtNsq][0] - Y[i][0]) + k*(Y[i+cbrtNsq][0] - Y[i][0]); }   // Middle layers, front edge
      else if (i == cn2) { return k*(Y[cn3][0] - Y[i][0]) + k*(Y[i+1][0] - Y[i][0]) + k*(Y[i-cbrtN][0] - Y[i][0]) + k*(Y[cn0    ][0] - Y[i][0]) + k*(Y[i-cbrtNsq][0] - Y[i][0]) + k*(Y[i+cbrtNsq][0] - Y[i][0]); }   // Middle layers, left edge
      else if (i == cn3) { return k*(Y[i-1][0] - Y[i][0]) + k*(Y[cn2][0] - Y[i][0]) + k*(Y[i-cbrtN][0] - Y[i][0]) + k*(Y[cn1    ][0] - Y[i][0]) + k*(Y[i-cbrtNsq][0] - Y[i][0]) + k*(Y[i+cbrtNsq][0] - Y[i][0]); }   // Middle layers, right edge

      // Middle layer sides / Cube face point (6 connections)
      else if (i > cn0 && i < cn1                    ) { return k*(Y[i-1      ][0] - Y[i][0]) + k*(Y[i+1      ][0] - Y[i][0]) + k*(Y[i-cbrtN+cbrtNsq][0] - Y[i][0]) + k*(Y[i+cbrtN        ][0] - Y[i][0]) + k*(Y[i-cbrtNsq][0] - Y[i][0]) + k*(Y[i+cbrtNsq][0] - Y[i][0]); }   // Middle layers. back edge
      else if (i > cn2 && i < cn3                    ) { return k*(Y[i-1      ][0] - Y[i][0]) + k*(Y[i+1      ][0] - Y[i][0]) + k*(Y[i-cbrtN        ][0] - Y[i][0]) + k*(Y[i+cbrtN-cbrtNsq][0] - Y[i][0]) + k*(Y[i-cbrtNsq][0] - Y[i][0]) + k*(Y[i+cbrtNsq][0] - Y[i][0]); }   // Middle layers, front edge
      else if (i%cbrtN == 0     && i > cn0 && i < cn2) { return k*(Y[i-1+cbrtN][0] - Y[i][0]) + k*(Y[i+1      ][0] - Y[i][0]) + k*(Y[i-cbrtN        ][0] - Y[i][0]) + k*(Y[i+cbrtN        ][0] - Y[i][0]) + k*(Y[i-cbrtNsq][0] - Y[i][0]) + k*(Y[i+cbrtNsq][0] - Y[i][0]); }   // Middle layers, left edge
      else if ((i+1)%cbrtN == 0 && i > cn1 && i < cn3) { return k*(Y[i-1      ][0] - Y[i][0]) + k*(Y[i+1-cbrtN][0] - Y[i][0]) + k*(Y[i-cbrtN        ][0] - Y[i][0]) + k*(Y[i+cbrtN        ][0] - Y[i][0]) + k*(Y[i-cbrtNsq][0] - Y[i][0]) + k*(Y[i+cbrtNsq][0] - Y[i][0]); }   // Middle layers, right edge

      // Interior Neurons (6 coupling terms)
      else {
        //printf("i=%d in layer %d : cn0=%d, cn1=%d, cn2=%d, cn3=%d\n", i, n, cn0, cn1, cn2, cn3);
        return I_ext + k*(Y[i+1][0] - Y[i][0]) + k*(Y[i-1][0] - Y[i][0]) + k*(Y[i-cbrtN][0] - Y[i][0]) + k*(Y[i+cbrtN][0] - Y[i][0]) + k*(Y[i-cbrtNsq][0] - Y[i][0]) + k*(Y[i+cbrtNsq][0] - Y[i][0]); }
  }
}




double bi_all_to_all_inward(double y[], int i) {
	int j;
	double K = 0.0;
	double k;
	if (i!=0) {
		for (j=0; j<neqn && i!=j; j+=2) {
			// Random coupling strength between 0 and 1
      k = 0.01*(rand() % 100);

			// Inwards (ie everyone else - current neurone)
			K += k * (y[j] - y[i]);
		}
	}
	// First neurone unaffected by others
	else {K = I_ext;}
  //printf("%f\n", K);
  
  // Decrease K by sqrt(2) for every diagonal
  // it is away from the current neurone
  //if (N>4) {K /= ((sqrtN-2)*4*sqrt2);}
  
	return K;
}

double mean_field_1D_circular(double k, double y[], int i) {
  double sum = 0.0;
  int j;
  
  // Change for nearest/next nearest neighbour
  // Nearest neighbour = 1, next = 2 etc
  int neigh = 8;
  
  if (neigh>N-1) {neigh = N-1;}
  
  // PERIODIC
  for (j=-neigh; j<neigh+1; j++) {sum += y[(j+N)%N];}
  
  sum /= (neigh+1);
  
  if (i/2!=0) {return sum - y[i];}
  else        {return I_ext + sum - y[i];}
}


double mean_field_1D(double k, double y[], int i) {
  double sum = 0.0;
  int j;
  
  // Change for nearest/next nearest neighbour
  // Nearest neighbour = 1, next = 2 etc
  int neigh = 1;
  int n;
  
  // NON PERIODIC
  if (i-2*neigh<0) {
    j = 0;
    n = 2*neigh;
  }
  else if (i+2*neigh>neqn-2) {
    j = neqn-2*neigh;
    n = neqn-2;
  }
  else {
    j = i-2*neigh;
    n = i+2*neigh;
  }
  
  for (j=j; j<n; j+=2) {sum += y[j];}
  
  sum /= (1+2*neigh);
  
  // PERIODIC - 2D, not connected on diagonals
  /*int p = i/2;
  for (j=0; j<neigh; j++) {
    sum += y[(p+1*(j+1))%N];
    sum += y[(p-1*(j+1))%N];
    sum += y[(p-(j+1)*sqrtN)%N];
    sum += y[(p+(j+1)*sqrtN)%N];
  }*/
  
  //printf("sum : %lf\n", sum);
  if (i/2!=0) {return sum - y[i];}
  else        {return I_ext + sum - y[i];}
}

double mean_field_2D_non_periodic(double k, double y[], int i) {
  i=i/2; // Convert i from "equation number" to "neuron number"
  double sum = 0.0;
  
  // Corners (2 coupling terms)
  // Coupling is currently not wrapped
  if      (i==0)          {sum = Y[1        ][0] + Y[sqrtN    ][0];}
  else if (i==sqrtN-1)    {sum = Y[sqrtN-2  ][0] + Y[2*sqrtN-1][0];}
  else if (i==N-sqrtN)    {sum = Y[N-2*sqrtN][0] + Y[N-sqrtN+1][0];}
  else if (i==N-1)        {sum = Y[N-2      ][0] + Y[N-sqrtN-1][0];}
  
  // Edges (2 coupling terms)
  else if (i>0 && i < sqrtN-1)   {sum = Y[i+1][0] + Y[i-1    ][0] + Y[i+sqrtN][0];}   // Top Edge
  else if (i>N-sqrtN && i<N-1)   {sum = Y[i+1][0] + Y[i-1    ][0] + Y[i-sqrtN][0];}   // Bottom Edge
  else if (i%sqrtN==0)           {sum = Y[i+1][0] + Y[i-sqrtN][0] + Y[i+sqrtN][0];}   // Left Edge
  else if ((i+1)%sqrtN==0)       {sum = Y[i-1][0] + Y[i-sqrtN][0] + Y[i+sqrtN][0];}   // Right Edge
  
  // Interior Neurons (4 coupling terms)
  else {sum = Y[i+1][0] + Y[i-1][0] + Y[i-sqrtN][0] + Y[i+sqrtN][0];}
  
  sum += y[i];
  
  sum = 5.0;
  
  // Change the number in the condition below to
  // alter where the driving current is applied
  // Remember i is neurone number, not index
  if (i!=1) {return sum - y[i];}
  else      {return I_ext + sum - y[i];}
}

double mean_field_2D_periodic(double k, double y[], int i) {
  i=i/2; // Convert i from "equation number" to "neuron number"
  double sum = 0.0;
  
  // Corners
  if      (i==0)          {sum = Y[1        ][0] + Y[sqrtN    ][0] + Y[sqrtN-1][0] + Y[N-sqrtN][0];}
  else if (i==sqrtN-1)    {sum = Y[sqrtN-2  ][0] + Y[2*sqrtN-1][0] + Y[0      ][0] + Y[N-1    ][0];}
  else if (i==N-sqrtN)    {sum = Y[N-2*sqrtN][0] + Y[N-sqrtN+1][0] + Y[N-1    ][0] + Y[0      ][0];}
  else if (i==N-1)        {sum = Y[N-2      ][0] + Y[N-sqrtN-1][0] + Y[N-sqrtN][0] + Y[sqrtN-1][0];}
  
  // Edges (4 coupling terms)
  else if (i>0 && i < sqrtN-1)   {sum = Y[i+1][0]   +   Y[i-1    ][0]   +   Y[i+sqrtN][0]   +   Y[i+N-sqrtN][0];}   // Top Edge
  else if (i>N-sqrtN && i<N-1)   {sum = Y[i+1][0]   +   Y[i-1    ][0]   +   Y[i-sqrtN][0]   +   Y[i-N+sqrtN][0];}   // Bottom Edge
  else if (i%sqrtN==0)           {sum = Y[i+1][0]   +   Y[i-sqrtN][0]   +   Y[i+sqrtN][0]   +   Y[i+sqrtN-1][0];}   // Left Edge
  else if ((i+1)%sqrtN==0)       {sum = Y[i-1][0]   +   Y[i-sqrtN][0]   +   Y[i+sqrtN][0]   +   Y[i-sqrtN+1][0];}   // Right Edge
  
  // Interior Neurons (4 coupling terms)
  else {sum = Y[i+1][0]   +   Y[i-1][0]   +   Y[i-sqrtN][0]   +   Y[i+sqrtN][0];}
  
  sum += y[i];
  
  sum /= 5.0;
  
  // Change the number in the condition below to
  // alter where the driving current is applied
  // Remember i is neurone number, not index
  if (i!=4) {return sum - y[i];}
  else      {return I_ext + sum - y[i];}
}

double uni_reverse(double k, double y[], int i) {
	if (i != 0) { return k * (y[i - 2] - y[i]); }
	else { return I_ext; }
}

double threshold_1D(double k, double y[], int i) {
  // Generates a random threshold between 0 and 2
  srand(time(NULL));
  double thresh = (rand() % 20) / 10;
  
  if (i==0) {return I_ext;}
  else if (i > 0 && y[i-2] > thresh) {return I_ext;}
  else {return 0;}
}

double bi_nearest_neighbour_1D(double k, double y[], int i) {
	//1D chain of neurones with free standing boundary
	if (i == 0) { return I_ext + k * (y[2] - y[0]); }
	else if (i == neqn-2) { return k * (y[neqn-2 - 2] - y[neqn-2]); }
	else { return k * (y[i+2] - y[i] + y[i-2] - y[i]); }
}

double bi_1D(double k, double y[], int i) {
  
  //1D chain of neurones with free standing boundary
  int j;
  
  // No. neurones either side of current neurone
  int neigh = 3;
  
  // If N even, final neurone is coupled to once only (14/1/19)
  if      (N%2==0 && neigh > ((N/2)-1)) {neigh = N/2;}
  else if (neigh > ((N-1)/2)) {neigh = (N-1)/2;}
  
  double sum = 0.0;
  
  // Set if (0==0) here and change if statement below
  // to get same as circular case
  if (i!=0) {
    for (j=-2*neigh; j<2*(neigh+1); j+=2) {
      if (i != (i+j+neqn)%neqn) {
        sum += y[(i+j+neqn)%neqn];
        //printf("%d : %d\n", i, (i+j+neqn)%neqn);
      }
    }
    
    sum -= 2*neigh*y[i];
    sum *= k;
  }
  else {sum = I_ext;}
  
  // Uncomment this and set 0==0 in above if
  // statement for same as bi_circular function
  //if (i==0) {sum += I_ext;}
  
  return sum;
}

double bi_nearest_neighbour_circular(double k, double y[], int i) {
	//1D chain of neurones in a ring
	if      (i == 0)      { return I_ext + k * (y[i+2] - y[i]) + k * (y[neqn-2] - y[i]); }
	else if (i == neqn-2) { return         k * (y[i-2] - y[i]) + k * (y[0     ] - y[i]); }
	else                  { return         k * (y[i+2] - y[i]) + k * (y[i-2   ] - y[i]); }
}

double uni_nearest_neighbour_circular(double k, double y[], int i) {
  //1D chain of neurones in a ring
  if (i == 0) { return I_ext + k * (y[neqn-2] - y[i]); }
  else        { return k * (y[i-2] - y[i]); }
}


double bi_nearest_neighbours_2D_periodic_inward(double k, double y[], int i) {
  i=i/2; // Convert i from "equation number" to "neuron number"
  
  // Corners (4 coupling terms)
  if      (i==0)          {return         k*(Y[1        ][0] - Y[0      ][0]) + k*(Y[sqrtN    ][0] - Y[0      ][0]) + k*(Y[sqrtN-1][0] - Y[0      ][0]) + k*(Y[N-sqrtN][0] - Y[0      ][0]);}
  else if (i==sqrtN-1)    {return         k*(Y[sqrtN-2  ][0] - Y[sqrtN-1][0]) + k*(Y[2*sqrtN-1][0] - Y[sqrtN-1][0]) + k*(Y[0      ][0] - Y[sqrtN-1][0]) + k*(Y[N-1    ][0] - Y[sqrtN-1][0]);}
  else if (i==N-sqrtN)    {return         k*(Y[N-2*sqrtN][0] - Y[N-sqrtN][0]) + k*(Y[N-sqrtN+1][0] - Y[N-sqrtN][0]) + k*(Y[N-1    ][0] - Y[N-sqrtN][0]) + k*(Y[0      ][0] - Y[N-sqrtN][0]);}
  else if (i==N-1)        {return         k*(Y[N-2      ][0] - Y[N-1    ][0]) + k*(Y[N-sqrtN-1][0] - Y[N-1    ][0]) + k*(Y[N-sqrtN][0] - Y[N-1    ][0]) + k*(Y[sqrtN-1][0] - Y[N-1    ][0]);}
  
  // Edges (4 coupling terms)
  else if (i>0 && i < sqrtN-1)   {return k*(Y[i+1][0] - Y[i][0])   +   k*(Y[i-1    ][0] - Y[i][0])   +   k*(Y[i+sqrtN][0] - Y[i][0])   +   k*(Y[i+N-sqrtN][0] - Y[0][0]);}   // Top Edge
  else if (i>N-sqrtN && i<N-1)   {return k*(Y[i+1][0] - Y[i][0])   +   k*(Y[i-1    ][0] - Y[i][0])   +   k*(Y[i-sqrtN][0] - Y[i][0])   +   k*(Y[i-N+sqrtN][0] - Y[0][0]);}   // Bottom Edge
  else if (i%sqrtN==0)           {return k*(Y[i+1][0] - Y[i][0])   +   k*(Y[i-sqrtN][0] - Y[i][0])   +   k*(Y[i+sqrtN][0] - Y[i][0])   +   k*(Y[i+sqrtN-1][0] - Y[0][0]);}   // Left Edge
  else if ((i+1)%sqrtN==0)       {return k*(Y[i-1][0] - Y[i][0])   +   k*(Y[i-sqrtN][0] - Y[i][0])   +   k*(Y[i+sqrtN][0] - Y[i][0])   +   k*(Y[i-sqrtN+1][0] - Y[0][0]);}   // Right Edge
  
  // Interior Neurons (4 coupling terms)
  else {return I_ext + k*(Y[i+1][0] - Y[i][0])   +   k*(Y[i-1][0] - Y[i][0])   +   k*(Y[i-sqrtN][0] - Y[i][0])   +   k*(Y[i+sqrtN][0] - Y[i][0]);}
}

double bi_nearest_neighbours_2D_non_periodic(double k, double y[], int i) {
	i=i/2; // Convert i from "equation number" to "neuron number"

	// Corners (2 coupling terms)
  // Coupling is currently not wrapped
	if 		  (i==0)          {return       k*(Y[1        ][0] - Y[0      ][0])   +   k*(Y[sqrtN    ][0] - Y[0      ][0]);}
	else if (i==sqrtN-1)		{return		    k*(Y[sqrtN-2  ][0] - Y[sqrtN-1][0])   +   k*(Y[2*sqrtN-1][0] - Y[sqrtN-1][0]);}
	else if (i==N-sqrtN)		{return		    k*(Y[N-2*sqrtN][0] - Y[N-sqrtN][0])   +   k*(Y[N-sqrtN+1][0] - Y[N-sqrtN][0]);}
	else if (i==N-1)        {return		    k*(Y[N-2      ][0] - Y[N-1    ][0])   +   k*(Y[N-sqrtN-1][0] - Y[N-1    ][0]);}
	
	// Edges (2 coupling terms)
	else if (i>0 && i < sqrtN-1) {return k*(Y[i+1][0] - Y[i][0])   +   k*(Y[i-1    ][0] - Y[i][0])   +   k*(Y[i+sqrtN][0] - Y[i][0]);}   // Top Edge
	else if (i>N-sqrtN && i<N-1) {return k*(Y[i+1][0] - Y[i][0])   +   k*(Y[i-1    ][0] - Y[i][0])   +   k*(Y[i-sqrtN][0] - Y[i][0]);}   // Bottom Edge
	else if (i%sqrtN==0) 		     {return k*(Y[i+1][0] - Y[i][0])   +   k*(Y[i-sqrtN][0] - Y[i][0])   +   k*(Y[i+sqrtN][0] - Y[i][0]);}   // Left Edge
	else if ((i+1)%sqrtN==0) 	   {return k*(Y[i-1][0] - Y[i][0])   +   k*(Y[i-sqrtN][0] - Y[i][0])   +   k*(Y[i+sqrtN][0] - Y[i][0]);}   // Right Edge
	
  // Curent is somewhere within the edges of the grid
  // Remember i is neurone number
  else if ((N%2==0 && ((i==(N+sqrtN)/2) || (i==(N+sqrtN)/2-1) || (i==(N+sqrtN)/2-sqrtN) || (i==(N+sqrtN)/2+-sqrtN-1))) || (N%2!=0 && i==(N-1)/2)) {
    return I_ext + k*(Y[i+1][0] - Y[i][0])   +   k*(Y[i-1][0] - Y[i][0])   +   k*(Y[i-sqrtN][0] - Y[i][0])   +   k*(Y[i+sqrtN][0] - Y[i][0]);}
  
	// Interior Neurons (4 coupling terms)
	else {return k*(Y[i+1][0] - Y[i][0])   +   k*(Y[i-1][0] - Y[i][0])   +   k*(Y[i-sqrtN][0] - Y[i][0])   +   k*(Y[i+sqrtN][0] - Y[i][0]);}
}

double Iext_Thresh(int i){
	// Use this coupling function to vary I_ext across the neurones to find a "threshold value" to see where the neuron starts firing
	return 0.015;
}

double white_noise() {
  return 0.0;
}

void f01(double t, double y[], double yp[]) {
	const double a = 0.75;
	const double b = 0.8;
	const double c = 3.0;
	const double c_inv = 1 / c;
	// I_ext now declared at top of file
	const double k = 0.113;
	int i;

	for (i=0; i<neqn; i+=2) {

		//coupling[i/2] = bi_all_to_all_inward(y, i);
		//coupling[i/2] = mean_field_1D(k, y, i);
    //coupling[i/2] = mean_field_1D_circular(k, y, i);
		//coupling[i/2] = mean_field_2D_non_periodic(k, y, i);
		//coupling[i/2] = mean_field_2D_periodic(k, y, i);
		//coupling[i/2] = uni_reverse(k, y, i);
		//coupling[i/2] = threshold_1D(k, y, i);
		//coupling[i/2] = bi_nearest_neighbour_1D(k, y, i);
		//coupling[i/2] = bi_nearest_neighbour_circular(k, y, i);
    //coupling[i/2] = bi_1D(k, y, i);
		//coupling[i/2] = uni_nearest_neighbour_circular(k, y, i);
		//coupling[i/2] = bi_nearest_neighbours_2D_non_periodic(k, y, i);
		//coupling[i/2] = bi_nearest_neighbours_2D_periodic_inward(k, y, i);
    //coupling[i/2] = bi_nearest_neighbours_3D_non_periodic(k, y, i);
    //coupling[i/2] = bi_nearest_neighbours_3D_periodic(k, y, i);
			
		yp[i] = c * (y[i] - 0.33333333*y[i] * y[i] * y[i] - y[i + 1] + uni_reverse(k, y, i));
		yp[i + 1] = c_inv * (y[i] + a - b * y[i + 1]);
		
		// Write new values to the matrix
		// i->i/2 converts i from "equation number" to "neuron number"
		Y[i/2][0] = y[i  ];
		Y[i/2][1] = y[i+1];
	}
  return;
}
