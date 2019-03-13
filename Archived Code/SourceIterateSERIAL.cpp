# include <cstdlib>
# include <iostream>
# include <iomanip>
# include <cmath>
#include <fstream>
#include <string>
#include <stdlib.h>
#include <ctime>
#include <time.h>

using namespace std;

# include "ode.hpp"

int main();
void test01();
void f01(double t, double y[], double yp[]);
double coupling(int i, int neqn, double y[]);
double thresh(int i, int neqn, double y[]);
double a2abidir(int i, int neqn, double y[]);
double nneigh(int i, int neqn, double y[]);
int** matrix_creation(int N, int neighbours);


//Global Variables
int N = 9;
int neighbours=1;
int neqn = 2*N;
double k;
double I;

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
	

	/*for (int rowM = 0; rowM < 9; rowM++)
	{
		for (int colM = 0; colM < 9; colM++)
		{
			cout << Y[rowM][colM] << ' ';
		}
		cout << endl;
	}*/
	
	cout << "\n";
	cout << "ODE_PRB\n";
	cout << "  C++ version\n";
	cout << "  Test the ODE library.\n";
	//Simulation Parameters
	double k_min = 0;
	double k_max = 1.0;
	double I_min = -0.34;
	double I_max = -1.4;
	int neighmax;



	//test01();
	while (N <= 225) {
		if (N % 2 != 0) { neighmax = (N - 1) / 2; }
		else { neighmax = (N - 2) / 2; }
		for (neighbours = 1; neighbours <= neighmax; neighbours++) {
			for (k = k_min; k <= k_max; k += 0.05) {
				for (I = I_min; I > I_max; I = I - 0.05) {
					test01();
				}
			}
		}
		N = (int) (sqrt(N) + 1)*(sqrt(N) + 1);
	}
	
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
	//int neqn = 18; //
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
	string filename = "FHN" + to_string(N) + to_string(neighbours)+ to_string(k) + to_string(I) + ".csv";
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
	
	

	cout << "\n";
	cout << "TEST01\n";
	cout << "  ODE solves a system of ordinary differential\n";
	cout << "  equations.\n";
	cout << "\n";
	cout << "      T           Y(1)         Y(2)\n";
	cout << "\n";

	abserr = 0.00001;
	relerr = 0.00001;

	iflag = 1;

	t = 0.0;
	y = new double[neqn];
	//INITIAL CONDITIONS
	for (i = 0; i < neqn; i++) { y[i] = 0.0; }
	

	cout << "  " << setw(8) << t
		<< "  " << setw(14) << y[0]
		<< "  " << setw(14) << y[1] << "\n";

	work = new double[100 + 21 * neqn];
	int fuckoff;
	for (i = 1; i <= step_num; i++)
	{
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
	srand(time(NULL));
	int i;
	
	double a = 0.75;
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

double thresh(int i, int neqn, double y[]) {
	srand(time(NULL));
	double xthresh = (rand() % 20) / 10;
	double I = -0.58;
	double II = -0.58;
	if (i == 0) { return I; }
	else {
		if (y[i - 2] > xthresh) {
			return II;
		}
		else { return 0; }
	}
}

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
	int** Y = matrix_creation(N, neighbours);
	int j;
	double K = 0;
	if (i != 0) {
		for (j = 0; j < N; j++) {
			if (Y[i/2][j] == 1) {
				K += k * (y[2 * j] - y[i]);
			}
			else { K += 0; }
		}
	}
	else { K = I; }
	//delete connection matrix
	int deleteint;
	for (deleteint = 0; deleteint < N; deleteint++) {
		delete[] Y[deleteint];
	}
	delete[] Y;
	//return coupling
	return K;
}
int** matrix_creation(int N, int neighbours) {
	int** Y = new int*[N];
	for (int i = 0; i < N; ++i)
		Y[i] = new int[N];
	int rowM;
	int colM;
	int neighcounter;
	for (rowM = 0; rowM < N; rowM++) {
		for (colM = 0; colM < N; colM++) {
			Y[rowM][colM] = 0;
		}
	}
	for (rowM = 0; rowM < N; rowM++) {
		for (neighcounter = -neighbours; neighcounter < neighbours + 1; neighcounter++) {
			if (neighcounter != 0) {
				if (rowM + neighcounter < 0) { Y[rowM][N + rowM + neighcounter] = 1; }
				else if (rowM + neighcounter > N - 1) { Y[rowM][(rowM + neighcounter) % N] = 1; }
				else { Y[rowM][rowM + neighcounter] = 1; }
			}
			else { Y[rowM][rowM + neighcounter] = 0; }
		}

	}
	return Y;
}
