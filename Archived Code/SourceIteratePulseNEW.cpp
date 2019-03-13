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

// Functions
int main();
void test01();
void f01(double t, double y[], double yp[]);
//double coupling(int i, int neqn, double y[]);
//double thresh(int i, int neqn, double y[]);
//double a2abidir(int i, int neqn, double y[]);
//double nneigh(int i, int neqn, double y[]);
//int** matrix_creation(int sqrtN, int neighbours);
//double I_t(double thresh);
//ector<int> locations(int locy);
//vector<double> PulseGenerator(double thresh);

// Global Variables
int N=4;
int sqrtN = (int)sqrt(N);
int neighbours;
int neqn = 2 * N;
int cbrtN = cbrt(N);
int cbrtNsq = cbrtN*cbrtN;
double sqrt2 = sqrt(2);
double k, I_ext, thresh;
int locy;
int systemtype;
double Y[4][2];
//vector<double> pulse_series;

// Variable controlling output frequency of data
// timejump = 1 => output every single data point
// Useful if you have low disk space
int timejump = 4;

// Make thresh threadprivate
#pragma omp threadprivate(thresh)

// Matrix of connections
string vname;
string command;


// System micro parameters
double a = 0.75;
double b = 0.8;
double c = 3.0;
double c_inv = 1/c;


//****************************************************************************

int main() {
  // Simulation Parameters
  int Nmin = 4;
  int Nmax = 4;
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
    
    N_inv = 1/N;
    
    if (N % 2 != 0) { neighmax = (N - 1) / 2; }
    else            { neighmax = (N - 2) / 2; }
    
    for (neighbours=1; neighbours<=neighmax; neighbours++) {      // Vary number of neighbours used in coupling
      //M = matrix_creation(sqrtN, neighbours);                     // Create matrix of connections
      for (locy=locmin; locy<=locmax; locy++) {                   // Vary density of drivers
        for (k=k_min; k<=k_max; k+=0.05) {                        // Vary the coupling strength
          
          thresh = thresh_min;
          for (j = jmin; j < jmax; j++) {                         // Vary threshold current
            thresh = j*0.05;                                      // Get thresh on each thread
            //pulse_series = PulseGenerator(thresh);                // Generate random pulse
            test01();                                             // Run simulation
          }
          
          //system("python PLOTIterate.py");                        // Call plotting program
          //command = "rm -f *.csv";                                // Delete all CSVs
          //system(command.c_str());
        }
      }
    }
    
    // Increment sqrtN
    sqrtN += 1;
    
    // N is now next square number
    N = (sqrtN+1)*(sqrtN+1);
    neqn = 2*N;
  
  /*----------*/
  /* Get time */
  /*----------*/
#ifndef _OPENMP
  /* Use clock_gettime */
  clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &tf);
  t_end = (double)tf.tv_sec + tf.tv_nsec / 1.0E9;
#endif
#ifdef _OPENMP
  /* Use OpenMP clock */
  t_end = omp_get_wtime();
#endif
  
  // Output time taken for current N to CSV
  timesFile << systemtype << "," << N << "," << (t_end-t_start) << ",";
  timesFile << endl;
  
  } // N

  return 0;
}
//****************************************************************************

void test01() {
  double abserr;
  int i, j;
  int iflag;
  int iwork[5];
  double pi = 3.141592653589793;
  double relerr;
  double step_size = 1.0e-5;
  double step_num = 1 / step_size;
  double t;
  double tf = 150;
  double tout;
  double *work;
  double *y;
  
  abserr = 0.00001;
  relerr = 0.00001;
  
  iflag = 1;
  
  t = 0.0;
  y = new double[neqn];
  
  //INITIAL CONDITIONS
  for (i = 0; i < neqn; i++) { y[i] = 0.0; }
  
  work = new double[100 + 21 * neqn];
  
  // Create results file
  ofstream vFile;
  vname = to_string(N) + "_" + to_string(neighbours) + "_" + to_string(locy) + "_" + to_string(k) + "_" + to_string(thresh) + ".csv";
  vFile.open(vname);
  
  // Time loop
  for (i=1; i<=step_num; i++) {
    
    if (0.01*(rand() % 100) > thresh) { I_ext = -1.0; };
    
    tout = (double)(i) * tf * step_size;
    
    ode(f01, neqn, y, t, tout, relerr, abserr, iflag, work, iwork);
    
    /* Assume our program is running correctly
    if (iflag != 2) {
      cout << "\n";
      cout << "TEST01 - Fatal error!\n";
      cout << "  ODE returned IFLAG = " << iflag << "\n";
      break;
    }
    */
    
    // Skip data points
    // Maintains accuracy of simulation whilst
    // keeping storage cost low
    //if (i%timejump==0) {
    
    // Output each neurone's v to CSV
    for (j=0; j<neqn; j+=2) { vFile << y[j] << ","; }
    vFile << endl;
    
    //}
  }
  
  delete[] work;
  delete[] y;
  
  vFile.close();
  
  return;
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

/* --------------------- 1D systems --------------------- */

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


// Bidirectional Coupling
double coupling(int i, int neqn, double y[]) {
  double I = -0.58;
  double k = 0.113;
  if (i == 0)             { return I + k*(y[i+2]-y[i]); }//+ k*(y[i+2] - 2*y[i] + y[neqn-2])
  else if (i == neqn - 2) { return     k*(y[i-2]-y[i]); }
  else                    { return     k*(y[i-2]-2*y[i]+y[i + 2]); }
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


/* --------------------- 1D Circular systems --------------------- */

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




/* --------------------- 2D systems --------------------- */

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
  if      (i==0)          {return        I_ext + k*(Y[1        ][0] - Y[0      ][0])   +   k*(Y[sqrtN    ][0] - Y[0      ][0]);}
  else if (i==sqrtN-1)    {return        k*(Y[sqrtN-2  ][0] - Y[sqrtN-1][0])   +   k*(Y[2*sqrtN-1][0] - Y[sqrtN-1][0]);}
  else if (i==N-sqrtN)    {return        k*(Y[N-2*sqrtN][0] - Y[N-sqrtN][0])   +   k*(Y[N-sqrtN+1][0] - Y[N-sqrtN][0]);}
  else if (i==N-1)        {return        k*(Y[N-2      ][0] - Y[N-1    ][0])   +   k*(Y[N-sqrtN-1][0] - Y[N-1    ][0]);}
  
  // Edges (2 coupling terms)
  else if (i>0 && i < sqrtN-1)  {return k*(Y[i+1][0] - Y[i][0])   +   k*(Y[i-1    ][0] - Y[i][0])   +   k*(Y[i+sqrtN][0] - Y[i][0]);}   // Top Edge
  else if (i>N-sqrtN && i<N-1)  {return k*(Y[i+1][0] - Y[i][0])   +   k*(Y[i-1    ][0] - Y[i][0])   +   k*(Y[i-sqrtN][0] - Y[i][0]);}   // Bottom Edge
  else if (i%sqrtN==0)          {return k*(Y[i+1][0] - Y[i][0])   +   k*(Y[i-sqrtN][0] - Y[i][0])   +   k*(Y[i+sqrtN][0] - Y[i][0]);}   // Left Edge
  else if ((i+1)%sqrtN==0)      {return k*(Y[i-1][0] - Y[i][0])   +   k*(Y[i-sqrtN][0] - Y[i][0])   +   k*(Y[i+sqrtN][0] - Y[i][0]);}   // Right Edge
  
  // Curent is somewhere within the edges of the grid
  // Remember i is neurone number
  else if ((N%2==0 && ((i==(N+sqrtN)/2) || (i==(N+sqrtN)/2-1) || (i==(N+sqrtN)/2-sqrtN) || (i==(N+sqrtN)/2+-sqrtN-1))) || (N%2!=0 && i==(N-1)/2)) {
    return k*(Y[i+1][0] - Y[i][0])   +   k*(Y[i-1][0] - Y[i][0])   +   k*(Y[i-sqrtN][0] - Y[i][0])   +   k*(Y[i+sqrtN][0] - Y[i][0]);}
  
  // Interior Neurons (4 coupling terms)
  else {return k*(Y[i+1][0] - Y[i][0])   +   k*(Y[i-1][0] - Y[i][0])   +   k*(Y[i-sqrtN][0] - Y[i][0])   +   k*(Y[i+sqrtN][0] - Y[i][0]);}
}




/* --------------------- 3D systems --------------------- */
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
  if      (i == c0) { return k*(Y[1  ][0] - Y[i][0]) + k*(Y[  cbrtN][0] - Y[i][0]) + k*(Y[  cbrtNsq][0] - Y[i][0]); } //0
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
  if      (i == 0 ) { return         k*(Y[c1 ][0] - Y[i][0]) + k*(Y[i+1][0] - Y[i][0]) + k*(Y[c2     ][0] - Y[i][0]) + k*(Y[  cbrtN][0] - Y[i][0]) + k*(Y[c4       ][0] - Y[i][0]) + k*(Y[  cbrtNsq][0] - Y[i][0]); } //0
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



/* --------------------- Mean Field systems --------------------- */

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




/* --------------------- More-than-nearest-neighbours systems --------------------- */


/* --------------------- Function to solve equations --------------------- */

void f01(double t, double y[], double yp[]) {
  
  int i;
  
#pragma omp parallel for private(i, y) schedule(dynamic) num_threads(omp_get_max_threads())
  for (i=0; i<neqn; i+=2) {
    yp[i] = c * (y[i] - 0.3333333333333*(y[i] * y[i] * y[i]) + y[i+1] + uni_reverse(k,y,i));//uni_reverse(k,y,i));
    yp[i+1] = -c_inv*(y[i] - a + b * y[i+1]);
  }
  return;
}
