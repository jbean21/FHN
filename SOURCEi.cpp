# include <cstdlib>
# include <iostream>
# include <iomanip>
# include <cmath>
# include <fstream>
# include <stdlib.h>
# include <string.h>
# include <ctime>
# include <omp.h>

using namespace std;

#include "ode.hpp"

int N;
int sqrtN;
int cbrtN;
int cbrtNsq;
double sqrt2 = sqrt(2);
int neqn;

// Neurone micro-parameters
const double a = 0.7;
const double b = 0.8;
const double c = 3.0;
const double c_inv = 1 / c;

double k;

// Change to 1.5 for driving from middle of 2D lattice
// 3D: I used I=2
double I_ext = 1;
#pragma omp threadprivate(I_ext)

int ndriv;
int timejump = 4;
// Declare an array on the stack
// (hence do not need to free)
int *pos;
double *Y;
double *coupling;
#pragma omp threadprivate(Y, coupling)

int neigh = 1;
double thresh = 0.0;
#pragma omp threadprivate(thresh)

int main();
void test01();
double Iext_Thresh(int i);
void f01(double t, double y[], double yp[]);

// Create file of failed runs
ofstream failedFile;

//****************************************************************************80

int main() {
  /* Change this depending on system */
  // 1D : 1
  // 2D : 2
  // 3D : 3
  /*----------------*/
   int dimension = 2;
  /*----------------*/
  
  int size;
  int nroot, dim;
  int e, i, j, r;
  int maxdriv, neighmax;
  srand(time(NULL));
  
  // Timing variables
  double t_start, t_end;
  
  // Create file of computation times
  //ofstream timesFile;
  //string timesfilename = "/media/jay/0FD90FF80FD90FF83/PROJECTDATA/2DPERIODIC/Times.csv";
  //timesFile.open(timesfilename);
  
  
  string failedfilename = "/media/jay/0FD90FF80FD90FF83/PROJECTDATA/2DPERIODIC_2/Failedruns.csv";
  failedFile.open(failedfilename);
  
  printf("\n---You are running a %dD system---\n", dimension);
  printf("Calculating ...\n");
  
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
  
  // Increments the nth root of the system
  // Raise nroot to the power of dimensions above to get N
  for (nroot=6; nroot<14; nroot++) {
    
    // Find N and neqn
    N = 1;
    for (dim=0; dim<dimension; dim++) { N *= nroot; }
    neqn = 2*N;
    
    // No need for an if statement or sqrt()/cbrt() : will slow computation
    // Therefore, set both sqrtN and cbrtN to be nroot as only 2D systems
    // will use sqrtN while 3D uses cbrtN
    sqrtN   = nroot;
    cbrtN   = nroot;
    cbrtNsq = nroot*nroot;
    
    //printf("Iteration %d, N : %d\n", nroot, N);
    
    // Define a sensible maximum number of drivers
    if (N%2==0) {
      maxdriv  = (int)N*0.5;
      neighmax = (int)N*0.5;
    }
    else        {
      maxdriv  = (int)(N-1)*0.5;
      neighmax = (int)(N-1)*0.5;
    }
    
    //for (neigh=1; neigh<neighmax; neigh++) {
      for (ndriv=1; ndriv<maxdriv; ndriv++) {
        
        // Memory allocation
        pos      = (int *)   malloc(ndriv *sizeof(int   )); // Driver positions
        
        // Initialise all driver positions to be outside the system
        memset(pos, -1, ndriv*sizeof(int));
        
        bool inArray = false;
        for (e=0; e<ndriv; e++) {
          while (true) {
            r = rand()%N;                                  // Generate a random driver location
            //printf("r : %d\n", r);
            for (int element=0; element<ndriv; element++) {
              if (r == pos[element]) {
                //printf("YES, r : %d, neurone : %d\n", r, pos[element]);
                inArray = true;
                break;
              }
              inArray = false;
            }
            if (inArray == false) {
              for (int element=0; element<ndriv; element++) {
                if (pos[element] == -1) {
                  pos[element] = r;
                  break;
                }
              }
              break;
            }
            //printf("Help\n");// Append to the array
          }
        }
        for (i=0; i<10; i++) {
          k = i*0.1;
#pragma omp parallel for private(j) schedule(dynamic) num_threads(omp_get_max_threads())
          for (j=0; j<10; j++) {
            thresh = j*0.1;
            
            // Can I move this memory allocation and below freeing to the N loop?
            Y        = (double*) calloc(N, sizeof(double)); // Storage for potentials
            coupling = (double*) calloc(N, sizeof(double)); // Storage for coupling
            
            // Pulse height:
            // 1D    : +1
            // 2D    : 
            // 3D    :
            // MF 2D :
            if ((rand()%10)*0.1>thresh) { I_ext = 1.0; }
            else                        { I_ext = 0.0; }
            //printf("thresh iteration : %d, thresh : %f\n", j, thresh);
            //for (int e=0; e<ndriv; e++) {printf("pos%d : %d\n", e, pos[e]);}
            test01();
            
            free(Y);
            free(coupling);
            
          } /* thresh */
        } /* k */
        
        free(pos);
        
      } /* ndriv */
    //} /* neigh */
    
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
    //timesFile << N << "," << (t_end-t_start);
    //timesFile << endl;
    
  }
  
  //timesFile.close();
  failedFile.close();
  
  printf("Ding! I'm done!\n");
  printf("Largest system - N : %d\n", N);
  printf("Time Elapsed : %.2f seconds\n", t_end-t_start);
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
	double tf = 150; //important
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
	for (i=0; i<neqn; i++) { y[i] = 0.0; }
	
	/* -------------------------------------------- */
  //////////////////////////
	/* Initialise CVS Files */
	//////////////////////////

  char vname[255], couplname[255], drivname[255];
  
  //printf("N : %d, neigh : %d, ndriv : %d, thresh : %f, thread : %d\n", N, neigh, ndriv, thresh, omp_get_thread_num());
	
  
  // Create CSV file for results to be written to
  ofstream vFile;
  sprintf(vname, "/media/jay/0FD90FF80FD90FF83/PROJECTDATA/2DPERIODIC_2/Data/v_%d_%d_%d_%.3f_%.3f.csv", N, neigh, ndriv, k, thresh);
  vFile.open(vname);
  //vFile << N << "," << neigh << "," << ndriv << "," << thresh << endl;

	// Create CSV file for coupling to be written to
	//ofstream couplingFile;
  //sprintf(couplname, "/media/jay/0FD90FF80FD90FF83/PROJECTDATA/2DNONPERIODIC/Data/coupl_%d_%d_%d_%.3f_%.3f.csv", N, neigh, ndriv, k, thresh);
  //couplingFile.open(couplname);
  //vFile << N << "," << neigh << "," << ndriv << "," << thresh << endl;
  
  // Create CSV file for driver locations
  ofstream driversFile;
  sprintf(drivname, "/media/jay/0FD90FF80FD90FF83/PROJECTDATA/2DPERIODIC_2/Data/driv_%d_%d_%d_%.3f_%.3f.csv", N, neigh, ndriv, k, thresh);
  driversFile.open(drivname);
  //vFile << N << "," << neigh << "," << ndriv << "," << thresh << endl;
  
  // Output coupling to CSV
  for (i=0; i<ndriv; i++) {
    driversFile << pos[i];
    driversFile << endl;
  }
  driversFile.close();
  
	/* -------------------------------------------- */

	work = new double[100 + 21 * neqn]; //want to work out what this does and why it uses 100 and 21
  	
	/* -------------------------------------------- */
	////////////////
	/* Time Loop */
	///////////////
	for (int timeindex = 1; timeindex <= step_num; timeindex++)
	{	
		tout = (double)(timeindex)* tf * stepSize;

		ode(f01, neqn, y, t, tout, relerr, abserr, iflag, work, iwork);

		if (iflag != 2)
		{
			cout << "\n";
			cout << "TEST01 - Fatal error!\n";
			cout << "  ODE returned IFLAG = " << iflag << "\n";
      
      failedFile << N << "," << neigh << "," << ndriv << "," << k << "," << thresh << "," << "FLAG:" << iflag << endl;
      
			break;
		}
    
    
		// Output 2D matrix to the next columns
    if (timeindex%timejump==0) {
      // Outputs NORMALISED TIME
      vFile << t*c << ",";
      for(i=0; i<N; i++) {
        vFile << Y[i] << ",";
      }
      vFile << endl;
      
      // Output coupling to CSV
      //for (i=0; i<N; i++) {
      //  couplingFile << k*coupling[i];
      //  couplingFile << ",";
      //}
      //couplingFile << endl;
    }
    
	}
	/* -------------------------------------------- */
  
  vFile.close();
  //couplingFile.close();
	
	delete[] work;
	delete[] y;

	return;
}

//****************************************************************************80

double bi_all_to_all_inward(double y[], int i) {
	int j;
	double K = 0.0;
	double k;
  
  for (j=0; j<neqn && i!=j; j+=2) {
    // Random coupling strength between 0 and 1
    k = 0.01*(rand() % 100);

    // Inwards (ie everyone else - current neurone)
    K += k * (y[j] - y[i]);
  }
  
  // Decrease K by sqrt(2) for every diagonal
  // it is away from the current neurone
  //if (N>4) {K /= ((sqrtN-2)*4*sqrt2);}
  
	return K;
}


double mean_field_2D_non_periodic(double k, double y[], int i) {
  i=i/2; // Convert i from "equation number" to "neuron number"
  double sum = 0.0;
  
  // Corners (2 coupling terms)
  // Coupling is currently not wrapped
  if      (i==0)          {sum = Y[1        ] + Y[sqrtN    ];}
  else if (i==sqrtN-1)    {sum = Y[sqrtN-2  ] + Y[2*sqrtN-1];}
  else if (i==N-sqrtN)    {sum = Y[N-2*sqrtN] + Y[N-sqrtN+1];}
  else if (i==N-1)        {sum = Y[N-2      ] + Y[N-sqrtN-1];}
  
  // Edges (2 coupling terms)
  else if (i>0 && i < sqrtN-1)   {sum = Y[i+1] + Y[i-1    ] + Y[i+sqrtN];}   // Top Edge
  else if (i>N-sqrtN && i<N-1)   {sum = Y[i+1] + Y[i-1    ] + Y[i-sqrtN];}   // Bottom Edge
  else if (i%sqrtN==0)           {sum = Y[i+1] + Y[i-sqrtN] + Y[i+sqrtN];}   // Left Edge
  else if ((i+1)%sqrtN==0)       {sum = Y[i-1] + Y[i-sqrtN] + Y[i+sqrtN];}   // Right Edge
  
  // Interior Neurons (4 coupling terms)
  else {sum = Y[i+1] + Y[i-1] + Y[i-sqrtN] + Y[i+sqrtN];}
  
  sum += y[i];
  
  // Average
  sum *= 0.2;
  
  return sum - y[i];
}

double mean_field_2D_periodic(double k, double y[], int i) {
  i=i/2; // Convert i from "equation number" to "neuron number"
  double sum = 0.0;
  
  // Corners
  if      (i==0)          {sum = Y[1        ] + Y[sqrtN    ] + Y[sqrtN-1] + Y[N-sqrtN];}
  else if (i==sqrtN-1)    {sum = Y[sqrtN-2  ] + Y[2*sqrtN-1] + Y[0      ] + Y[N-1    ];}
  else if (i==N-sqrtN)    {sum = Y[N-2*sqrtN] + Y[N-sqrtN+1] + Y[N-1    ] + Y[0      ];}
  else if (i==N-1)        {sum = Y[N-2      ] + Y[N-sqrtN-1] + Y[N-sqrtN] + Y[sqrtN-1];}
  
  // Edges (4 coupling terms)
  else if (i>0 && i < sqrtN-1)   {sum = Y[i+1]   +   Y[i-1    ]   +   Y[i+sqrtN]   +   Y[i+N-sqrtN];}   // Top Edge
  else if (i>N-sqrtN && i<N-1)   {sum = Y[i+1]   +   Y[i-1    ]   +   Y[i-sqrtN]   +   Y[i-N+sqrtN];}   // Bottom Edge
  else if (i%sqrtN==0)           {sum = Y[i+1]   +   Y[i-sqrtN]   +   Y[i+sqrtN]   +   Y[i+sqrtN-1];}   // Left Edge
  else if ((i+1)%sqrtN==0)       {sum = Y[i-1]   +   Y[i-sqrtN]   +   Y[i+sqrtN]   +   Y[i-sqrtN+1];}   // Right Edge
  
  // Interior Neurons (4 coupling terms)
  else {sum = Y[i+1]   +   Y[i-1]   +   Y[i-sqrtN]   +   Y[i+sqrtN];}
  
  sum += y[i];
  
  // Average
  sum *= 0.2;
  
  return sum - y[i];
}

double bi_1D_vary_num_neighbours(double k, double y[], int i) {
  
  //1D chain of neurones with free standing boundary
  int j;
  
  // neigh is no. neurones either side of current neurone
  
  double sum = 0.0;
  
  // Set if (0==0) here and change if statement below
  // to get same as circular case
  for (j=-2*neigh; j<=2*neigh; j+=2) {
    if (i != (i+j+neqn)%neqn) {
      sum += y[(i+j+neqn)%neqn];
      //printf("%d : %d\n", i, (i+j+neqn)%neqn);
    }
  }
  
  sum -= 2*neigh*y[i];
  
  return sum;
}


double bi_nearest_neighbours_2D_periodic(double k, double y[], int i) {
  i=i/2; // Convert i from "equation number" to "neuron number"
  
  // Corners (4 coupling terms)
  if      (i==0)          {return         (Y[1        ] - Y[0      ]) + (Y[sqrtN    ] - Y[0      ]) + (Y[sqrtN-1] - Y[0      ]) + (Y[N-sqrtN] - Y[0      ]);}
  else if (i==sqrtN-1)    {return         (Y[sqrtN-2  ] - Y[sqrtN-1]) + (Y[2*sqrtN-1] - Y[sqrtN-1]) + (Y[0      ] - Y[sqrtN-1]) + (Y[N-1    ] - Y[sqrtN-1]);}
  else if (i==N-sqrtN)    {return         (Y[N-2*sqrtN] - Y[N-sqrtN]) + (Y[N-sqrtN+1] - Y[N-sqrtN]) + (Y[N-1    ] - Y[N-sqrtN]) + (Y[0      ] - Y[N-sqrtN]);}
  else if (i==N-1)        {return         (Y[N-2      ] - Y[N-1    ]) + (Y[N-sqrtN-1] - Y[N-1    ]) + (Y[N-sqrtN] - Y[N-1    ]) + (Y[sqrtN-1] - Y[N-1    ]);}
  
  // Edges (4 coupling terms)
  else if (i>0 && i < sqrtN-1)   {return (Y[i+1] - Y[i])   +   (Y[i-1    ] - Y[i])   +   (Y[i+sqrtN] - Y[i])   +   (Y[i+N-sqrtN] - Y[0]);}   // Top Edge
  else if (i>N-sqrtN && i<N-1)   {return (Y[i+1] - Y[i])   +   (Y[i-1    ] - Y[i])   +   (Y[i-sqrtN] - Y[i])   +   (Y[i-N+sqrtN] - Y[0]);}   // Bottom Edge
  else if (i%sqrtN==0)           {return (Y[i+1] - Y[i])   +   (Y[i-sqrtN] - Y[i])   +   (Y[i+sqrtN] - Y[i])   +   (Y[i+sqrtN-1] - Y[0]);}   // Left Edge
  else if ((i+1)%sqrtN==0)       {return (Y[i-1] - Y[i])   +   (Y[i-sqrtN] - Y[i])   +   (Y[i+sqrtN] - Y[i])   +   (Y[i-sqrtN+1] - Y[0]);}   // Right Edge
  
  // Interior Neurons (4 coupling terms)
  else {return  (Y[i+1] - Y[i])   +   (Y[i-1] - Y[i])   +   (Y[i-sqrtN] - Y[i])   +   (Y[i+sqrtN] - Y[i]);}
}

double bi_nearest_neighbours_2D_non_periodic(double k, double y[], int i) {
	i=i/2; // Convert i from "equation number" to "neuron number"

	// Corners (2 coupling terms)
  // Coupling is currently not wrapped
	if 		  (i==0)          {return       (Y[1        ] - Y[0      ])   +   (Y[sqrtN    ] - Y[0      ]);}
	else if (i==sqrtN-1)		{return		    (Y[sqrtN-2  ] - Y[sqrtN-1])   +   (Y[2*sqrtN-1] - Y[sqrtN-1]);}
	else if (i==N-sqrtN)		{return		    (Y[N-2*sqrtN] - Y[N-sqrtN])   +   (Y[N-sqrtN+1] - Y[N-sqrtN]);}
	else if (i==N-1)        {return		    (Y[N-2      ] - Y[N-1    ])   +   (Y[N-sqrtN-1] - Y[N-1    ]);}
	
	// Edges (2 coupling terms)
	else if (i>0 && i < sqrtN-1) {return (Y[i+1] - Y[i])   +   (Y[i-1    ] - Y[i])   +   (Y[i+sqrtN] - Y[i]);}   // Top Edge
	else if (i>N-sqrtN && i<N-1) {return (Y[i+1] - Y[i])   +   (Y[i-1    ] - Y[i])   +   (Y[i-sqrtN] - Y[i]);}   // Bottom Edge
	else if (i%sqrtN==0) 		     {return (Y[i+1] - Y[i])   +   (Y[i-sqrtN] - Y[i])   +   (Y[i+sqrtN] - Y[i]);}   // Left Edge
	else if ((i+1)%sqrtN==0) 	   {return (Y[i-1] - Y[i])   +   (Y[i-sqrtN] - Y[i])   +   (Y[i+sqrtN] - Y[i]);}   // Right Edge
	
  // Curent is somewhere within the edges of the grid
  // Remember i is neurone number
  else if ((N%2==0 && ((i==(N+sqrtN)/2) || (i==(N+sqrtN)/2-1) || (i==(N+sqrtN)/2-sqrtN) || (i==(N+sqrtN)/2+-sqrtN-1))) || (N%2!=0 && i==(N-1)/2)) {
    return (Y[i+1] - Y[i])   +   (Y[i-1] - Y[i])   +   (Y[i-sqrtN] - Y[i])   +   (Y[i+sqrtN] - Y[i]);}
  
	// Interior Neurons (4 coupling terms)
	else {return (Y[i+1] - Y[i])   +   (Y[i-1] - Y[i])   +   (Y[i-sqrtN] - Y[i])   +   (Y[i+sqrtN] - Y[i]);}
}

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
  if      (i == c0) { return         (Y[1  ] - Y[i]) + (Y[  cbrtN] - Y[i]) + (Y[  cbrtNsq] - Y[i]); } //0
  else if (i == c1) { return         (Y[i-1] - Y[i]) + (Y[i+cbrtN] - Y[i]) + (Y[i+cbrtNsq] - Y[i]); } //2
  else if (i == c2) { return         (Y[i+1] - Y[i]) + (Y[i-cbrtN] - Y[i]) + (Y[i+cbrtNsq] - Y[i]); } //6
  else if (i == c3) { return         (Y[i-1] - Y[i]) + (Y[i-cbrtN] - Y[i]) + (Y[i+cbrtNsq] - Y[i]); } //8
  else if (i == c4) { return         (Y[i+1] - Y[i]) + (Y[i+cbrtN] - Y[i]) + (Y[i-cbrtNsq] - Y[i]); } //18
  else if (i == c5) { return         (Y[i-1] - Y[i]) + (Y[i+cbrtN] - Y[i]) + (Y[i-cbrtNsq] - Y[i]); } //20
  else if (i == c6) { return         (Y[i+1] - Y[i]) + (Y[i-cbrtN] - Y[i]) + (Y[i-cbrtNsq] - Y[i]); } //24
  else if (i == c7) { return         (Y[i-1] - Y[i]) + (Y[i-cbrtN] - Y[i]) + (Y[i-cbrtNsq] - Y[i]); } //26
  
  // Cube edges (4 coupling terms)
  else if (i > 0            && i < c1) { return (Y[i+1] - Y[i]) + (Y[i-1    ] - Y[i]) + (Y[i+cbrtN] - Y[i]) + (Y[i+cbrtNsq] - Y[i]); }   // Top layer. back edge
  else if (i > c2           && i < c3) { return (Y[i+1] - Y[i]) + (Y[i-1    ] - Y[i]) + (Y[i-cbrtN] - Y[i]) + (Y[i+cbrtNsq] - Y[i]); }   // Top layer, front edge
  else if (i%cbrtN == 0     && i < c2) { return (Y[i+1] - Y[i]) + (Y[i-cbrtN] - Y[i]) + (Y[i+cbrtN] - Y[i]) + (Y[i+cbrtNsq] - Y[i]); }   // Top layer, left edge
  else if ((i+1)%cbrtN == 0 && i < c3) { return (Y[i-1] - Y[i]) + (Y[i-cbrtN] - Y[i]) + (Y[i+cbrtN] - Y[i]) + (Y[i+cbrtNsq] - Y[i]); }   // Top layer, right edge
  
  else if (i > c4           && i < c5          ) { return (Y[i+1] - Y[i]) + (Y[i-1    ] - Y[i]) + (Y[i+cbrtN] - Y[i]) + (Y[i-cbrtNsq] - Y[i]); }   // Bottom layer. back edge (18, 20)
  else if (i > c6           && i < c7          ) { return (Y[i+1] - Y[i]) + (Y[i-1    ] - Y[i]) + (Y[i-cbrtN] - Y[i]) + (Y[i-cbrtNsq] - Y[i]); }   // Bottom layer, front edge (24, 26)
  else if (i%cbrtN == 0     && i > c4 && i < c6) { return (Y[i+1] - Y[i]) + (Y[i-cbrtN] - Y[i]) + (Y[i+cbrtN] - Y[i]) + (Y[i-cbrtNsq] - Y[i]); }   // Bottom layer, left edge
  else if ((i+1)%cbrtN == 0 && i > c5 && i < c7) { return (Y[i-1] - Y[i]) + (Y[i-cbrtN] - Y[i]) + (Y[i+cbrtN] - Y[i]) + (Y[i-cbrtNsq] - Y[i]); }   // Bottom layer, right edge
  
  /* ---------------------------- */
  
  // Cube face points (5 connections) -- cube edges caught by upper if statements
  else if (i>cbrtN && i<cbrtNsq-cbrtN)       { return (Y[i-1] - Y[i]) + (Y[i+1] - Y[i]) + (Y[i-cbrtN] - Y[i]) + (Y[i+cbrtN] - Y[i]) + (Y[i+cbrtNsq] - Y[i]); }   // Top layer, face points
  else if (i>N-cbrtNsq+cbrtN && i<N-cbrtN-1) { return (Y[i-1] - Y[i]) + (Y[i+1] - Y[i]) + (Y[i-cbrtN] - Y[i]) + (Y[i+cbrtN] - Y[i]) + (Y[i-cbrtNsq] - Y[i]); }   // Bottom layer, face points
  
  
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
    if      (i == cn0) { return (Y[i+1] - Y[i]) + (Y[i+cbrtN] - Y[i]) + (Y[i-cbrtNsq] - Y[i]) + (Y[i+cbrtNsq] - Y[i]); }   // Middle layers. back edge
    else if (i == cn1) { return (Y[i-1] - Y[i]) + (Y[i+cbrtN] - Y[i]) + (Y[i-cbrtNsq] - Y[i]) + (Y[i+cbrtNsq] - Y[i]); }   // Middle layers, front edge
    else if (i == cn2) { return (Y[i+1] - Y[i]) + (Y[i-cbrtN] - Y[i]) + (Y[i-cbrtNsq] - Y[i]) + (Y[i+cbrtNsq] - Y[i]); }   // Middle layers, left edge
    else if (i == cn3) { return (Y[i-1] - Y[i]) + (Y[i-cbrtN] - Y[i]) + (Y[i-cbrtNsq] - Y[i]) + (Y[i+cbrtNsq] - Y[i]); }   // Middle layers, right edge
    
    // Middle layer sides / Cube face point (5 connections)
    else if (i > cn0 && i < cn1                    ) { return (Y[i+1] - Y[i]) + (Y[i-1    ] - Y[i]) + (Y[i+cbrtN] - Y[i]) + (Y[i+cbrtNsq] - Y[i]) + (Y[i-cbrtNsq] - Y[i]); }   // Middle layers. back edge
    else if (i > cn2 && i < cn3                    ) { return (Y[i+1] - Y[i]) + (Y[i-1    ] - Y[i]) + (Y[i-cbrtN] - Y[i]) + (Y[i+cbrtNsq] - Y[i]) + (Y[i-cbrtNsq] - Y[i]); }   // Middle layers, front edge
    else if (i%cbrtN == 0     && i > cn0 && i < cn2) { return (Y[i+1] - Y[i]) + (Y[i-cbrtN] - Y[i]) + (Y[i+cbrtN] - Y[i]) + (Y[i+cbrtNsq] - Y[i]) + (Y[i-cbrtNsq] - Y[i]); }   // Middle layers, left edge
    else if ((i+1)%cbrtN == 0 && i > cn1 && i < cn3) { return (Y[i-1] - Y[i]) + (Y[i-cbrtN] - Y[i]) + (Y[i+cbrtN] - Y[i]) + (Y[i+cbrtNsq] - Y[i]) + (Y[i-cbrtNsq] - Y[i]); }   // Middle layers, right edge
    
    // Interior Neurons (6 coupling terms)
    else {
      //printf("i=%d in layer %d : cn0=%d, cn1=%d, cn2=%d, cn3=%d\n", i, n, cn0, cn1, cn2, cn3);
      return  (Y[i+1] - Y[i]) + (Y[i-1] - Y[i]) + (Y[i-cbrtN] - Y[i]) + (Y[i+cbrtN] - Y[i]) + (Y[i-cbrtNsq] - Y[i]) + (Y[i+cbrtNsq] - Y[i]); }
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
  if      (i == 0 ) { return         (Y[c1 ] - Y[i]) + (Y[i+1] - Y[i]) + (Y[c2     ] - Y[i]) + (Y[  cbrtN] - Y[i]) + (Y[c4       ] - Y[i]) + (Y[  cbrtNsq] - Y[i]); } //0
  else if (i == c1) { return         (Y[i-1] - Y[i]) + (Y[c0 ] - Y[i]) + (Y[c3     ] - Y[i]) + (Y[i+cbrtN] - Y[i]) + (Y[c5       ] - Y[i]) + (Y[i+cbrtNsq] - Y[i]); } //2
  else if (i == c2) { return         (Y[c3 ] - Y[i]) + (Y[i+1] - Y[i]) + (Y[i-cbrtN] - Y[i]) + (Y[c0     ] - Y[i]) + (Y[c6       ] - Y[i]) + (Y[i+cbrtNsq] - Y[i]); } //6
  else if (i == c3) { return         (Y[i-1] - Y[i]) + (Y[c2 ] - Y[i]) + (Y[i-cbrtN] - Y[i]) + (Y[c1     ] - Y[i]) + (Y[c7       ] - Y[i]) + (Y[i+cbrtNsq] - Y[i]); } //8
  else if (i == c4) { return         (Y[c5 ] - Y[i]) + (Y[i+1] - Y[i]) + (Y[c6     ] - Y[i]) + (Y[i+cbrtN] - Y[i]) + (Y[i-cbrtNsq] - Y[i]) + (Y[c0       ] - Y[i]); } //18
  else if (i == c5) { return         (Y[i-1] - Y[i]) + (Y[c4 ] - Y[i]) + (Y[c7     ] - Y[i]) + (Y[i+cbrtN] - Y[i]) + (Y[i-cbrtNsq] - Y[i]) + (Y[c1       ] - Y[i]); } //20
  else if (i == c6) { return         (Y[c7 ] - Y[i]) + (Y[i+1] - Y[i]) + (Y[i-cbrtN] - Y[i]) + (Y[c4     ] - Y[i]) + (Y[i-cbrtNsq] - Y[i]) + (Y[c2       ] - Y[i]); } //24
  else if (i == c7) { return         (Y[i-1] - Y[i]) + (Y[c6 ] - Y[i]) + (Y[i-cbrtN] - Y[i]) + (Y[c5     ] - Y[i]) + (Y[i-cbrtNsq] - Y[i]) + (Y[c3       ] - Y[i]); } //26
  
  // Cube edges (6 coupling terms)
  else if (i > 0            && i < c1)           { return (Y[i-1      ] - Y[i]) + (Y[i+1      ] - Y[i]) + (Y[i-cbrtN+cbrtNsq] - Y[i]) + (Y[i+cbrtN        ] - Y[i]) + (Y[i+cbrtNsq] - Y[i]) + (Y[i+(cbrtN-1)*cbrtNsq] - Y[i]); }   // Top layer. back edge
  else if (i > c2           && i < c3)           { return (Y[i-1      ] - Y[i]) + (Y[i+1      ] - Y[i]) + (Y[i-cbrtN        ] - Y[i]) + (Y[i+cbrtN-cbrtNsq] - Y[i]) + (Y[i+cbrtNsq] - Y[i]) + (Y[i+(cbrtN-1)*cbrtNsq] - Y[i]); }   // Top layer, front edge
  else if (i%cbrtN == 0     && i < c2)           { return (Y[i-1+cbrtN] - Y[i]) + (Y[i+1      ] - Y[i]) + (Y[i-cbrtN        ] - Y[i]) + (Y[i+cbrtN        ] - Y[i]) + (Y[i+cbrtNsq] - Y[i]) + (Y[i+(cbrtN-1)*cbrtNsq] - Y[i]); }   // Top layer, left edge
  else if ((i+1)%cbrtN == 0 && i < c3)           { return (Y[i-1      ] - Y[i]) + (Y[i-cbrtN+1] - Y[i]) + (Y[i-cbrtN        ] - Y[i]) + (Y[i+cbrtN        ] - Y[i]) + (Y[i+cbrtNsq] - Y[i]) + (Y[i+(cbrtN-1)*cbrtNsq] - Y[i]); }   // Top layer, right edge
  
  else if (i > c4           && i < c5          ) { return (Y[i-1      ] - Y[i]) + (Y[i+1      ] - Y[i]) + (Y[i-cbrtN-cbrtNsq] - Y[i]) + (Y[i+cbrtN        ] - Y[i]) + (Y[i-cbrtNsq] - Y[i]) + (Y[i-(cbrtN-1)*cbrtNsq] - Y[i]); }   // Bottom layer. back edge (18, 20)
  else if (i > c6           && i < c7          ) { return (Y[i-1      ] - Y[i]) + (Y[i+1      ] - Y[i]) + (Y[i-cbrtN        ] - Y[i]) + (Y[i+cbrtN-cbrtNsq] - Y[i]) + (Y[i-cbrtNsq] - Y[i]) + (Y[i-(cbrtN-1)*cbrtNsq] - Y[i]); }   // Bottom layer, front edge (24, 26)
  else if (i%cbrtN == 0     && i > c4 && i < c6) { return (Y[i-1+cbrtN] - Y[i]) + (Y[i+1      ] - Y[i]) + (Y[i-cbrtN        ] - Y[i]) + (Y[i+cbrtN        ] - Y[i]) + (Y[i-cbrtNsq] - Y[i]) + (Y[i-(cbrtN-1)*cbrtNsq] - Y[i]); }   // Bottom layer, left edge
  else if ((i+1)%cbrtN == 0 && i > c5 && i < c7) { return (Y[i-1      ] - Y[i]) + (Y[i+1-cbrtN] - Y[i]) + (Y[i-cbrtN        ] - Y[i]) + (Y[i+cbrtN        ] - Y[i]) + (Y[i-cbrtNsq] - Y[i]) + (Y[i-(cbrtN-1)*cbrtNsq] - Y[i]); }   // Bottom layer, right edge
  
  /* ---------------------------- */
  
  // Cube face points (5 connections) -- cube edges caught by upper if statements
  else if (i>cbrtN && i<cbrtNsq-cbrtN)       { return (Y[i-1] - Y[i]) + (Y[i+1] - Y[i]) + (Y[i-cbrtN] - Y[i]) + (Y[i+cbrtN] - Y[i]) + (Y[i+(cbrtN-1)*cbrtNsq] - Y[i]) + (Y[i+cbrtNsq          ] - Y[i]); }   // Top layer, face points
  else if (i>N-cbrtNsq+cbrtN && i<N-cbrtN-1) { return (Y[i-1] - Y[i]) + (Y[i+1] - Y[i]) + (Y[i-cbrtN] - Y[i]) + (Y[i+cbrtN] - Y[i]) + (Y[i-cbrtNsq          ] - Y[i]) + (Y[i-(cbrtN-1)*cbrtNsq] - Y[i]); }   // Bottom layer, face points
  
  
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
    if      (i == cn0) { return (Y[cn1] - Y[i]) + (Y[i+1] - Y[i]) + (Y[cn2    ] - Y[i]) + (Y[i+cbrtN] - Y[i]) + (Y[i-cbrtNsq] - Y[i]) + (Y[i+cbrtNsq] - Y[i]); }   // Middle layers. back edge
    else if (i == cn1) { return (Y[i-1] - Y[i]) + (Y[cn0] - Y[i]) + (Y[cn3    ] - Y[i]) + (Y[i+cbrtN] - Y[i]) + (Y[i-cbrtNsq] - Y[i]) + (Y[i+cbrtNsq] - Y[i]); }   // Middle layers, front edge
    else if (i == cn2) { return (Y[cn3] - Y[i]) + (Y[i+1] - Y[i]) + (Y[i-cbrtN] - Y[i]) + (Y[cn0    ] - Y[i]) + (Y[i-cbrtNsq] - Y[i]) + (Y[i+cbrtNsq] - Y[i]); }   // Middle layers, left edge
    else if (i == cn3) { return (Y[i-1] - Y[i]) + (Y[cn2] - Y[i]) + (Y[i-cbrtN] - Y[i]) + (Y[cn1    ] - Y[i]) + (Y[i-cbrtNsq] - Y[i]) + (Y[i+cbrtNsq] - Y[i]); }   // Middle layers, right edge
    
    // Middle layer sides / Cube face point (6 connections)
    else if (i > cn0 && i < cn1                    ) { return (Y[i-1      ] - Y[i]) + (Y[i+1      ] - Y[i]) + (Y[i-cbrtN+cbrtNsq] - Y[i]) + (Y[i+cbrtN        ] - Y[i]) + (Y[i-cbrtNsq] - Y[i]) + (Y[i+cbrtNsq] - Y[i]); }   // Middle layers. back edge
    else if (i > cn2 && i < cn3                    ) { return (Y[i-1      ] - Y[i]) + (Y[i+1      ] - Y[i]) + (Y[i-cbrtN        ] - Y[i]) + (Y[i+cbrtN-cbrtNsq] - Y[i]) + (Y[i-cbrtNsq] - Y[i]) + (Y[i+cbrtNsq] - Y[i]); }   // Middle layers, front edge
    else if (i%cbrtN == 0     && i > cn0 && i < cn2) { return (Y[i-1+cbrtN] - Y[i]) + (Y[i+1      ] - Y[i]) + (Y[i-cbrtN        ] - Y[i]) + (Y[i+cbrtN        ] - Y[i]) + (Y[i-cbrtNsq] - Y[i]) + (Y[i+cbrtNsq] - Y[i]); }   // Middle layers, left edge
    else if ((i+1)%cbrtN == 0 && i > cn1 && i < cn3) { return (Y[i-1      ] - Y[i]) + (Y[i+1-cbrtN] - Y[i]) + (Y[i-cbrtN        ] - Y[i]) + (Y[i+cbrtN        ] - Y[i]) + (Y[i-cbrtNsq] - Y[i]) + (Y[i+cbrtNsq] - Y[i]); }   // Middle layers, right edge
    
    // Interior Neurons (6 coupling terms)
    else {
      //printf("i=%d in layer %d : cn0=%d, cn1=%d, cn2=%d, cn3=%d\n", i, n, cn0, cn1, cn2, cn3);
      return  (Y[i+1] - Y[i]) + (Y[i-1] - Y[i]) + (Y[i-cbrtN] - Y[i]) + (Y[i+cbrtN] - Y[i]) + (Y[i-cbrtNsq] - Y[i]) + (Y[i+cbrtNsq] - Y[i]); }
  }
}


void f01(double t, double y[], double yp[]) {
  // Neurone number, Equation number
  int i, i2;
  
  // Current variable
  double I_var;
  int element;
  
  // Loop over neurone number
	for (i=0; i<N; i++) {
    
    // Find v equation number
    i2 = 2*i;
    
    /* Uncomment to choose coupling type */
		//coupling[i] = bi_all_to_all_inward(y, i2);
    //coupling[i] = bi_1D_vary_num_neighbours(k, y, i2);
    coupling[i] = bi_nearest_neighbours_2D_periodic(k, y, i2);        /* DONE */
		//coupling[i] = bi_nearest_neighbours_2D_non_periodic(k, y, i2);
    //coupling[i] = bi_nearest_neighbours_3D_non_periodic(k, y, i2);
    //coupling[i] = bi_nearest_neighbours_3D_periodic(k, y, i2);
		//coupling[i] = mean_field_2D_non_periodic(k, y, i2);
		//coupling[i] = mean_field_2D_periodic(k, y, i2);
    
    I_var = 0.0;
    
    // Determine which neruones should have the driving current applied to them
    for (element=0; element<ndriv; element++) {
      if (i == pos[element]) {
        I_var = I_ext;
      }
    }
			
		yp[i2] = c * (y[i2] - 0.33333333*y[i2] * y[i2] * y[i2] - y[i2 + 1] + I_var + k*coupling[i]);
		yp[i2 + 1] = c_inv * (y[i2] + a - b * y[i2 + 1]);
		
		// Write new values to the matrix
		// i2->i converts i2 from "equation number" to "neuron number"
		Y[i] = y[i2];
	}
  return;
}
