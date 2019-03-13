# include <cstdlib>
# include <iostream>
# include <iomanip>
# include <cmath>
#include <ctime>
# include <fstream>

using namespace std;

# include "ode.hpp"

int Nneurons = 4;
int neqn = 2 * Nneurons;
double random_driver(double t);
int main();
void test01();
void f01(double t, double y[], double yp[]);

//****************************************************************************80

int main() {
	
	cout << "\n";
	cout << "ODE_PRB\n";
	cout << "  C++ version\n";
	cout << "  Test the ODE library.\n";

	test01();
	//
	//  Terminate.
	//
	  //GOTTA GO FAST
	cout << "\n";
	cout << "ODE_PRB\n";
	cout << "  Normal end of execution.\n";
	cout << "\n";
	

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
	double stepSize = 1e-4; //important
	double step_num = 1 / stepSize; //unsure about this
	double t;
	double tout;
	double *work;
	double *y;

	cout << "\n";
	cout << "TEST01\n";
	cout << "  ODE solves a system of ordinary differential\n";
	cout << "  equations.\n";
	cout << "\n";
	cout << "      T              ";
	for (int i = 1; i < Nneurons + 1; i++) {
		cout << "v";
		cout << i;
		cout << "              ";
		cout << "w";
		cout << i;
		cout << "              ";
	}
	cout << "\n\n";

	abserr = 0.00001;
	relerr = 0.00001;

	iflag = 1;

	t = 0.0;
	y = new double[neqn];

	// Initial Conditions
	for (i = 0; i < neqn; i++) { y[i] = 0.0; }

	//Create CSV file for results to be written to
	ofstream resultsFile;
	resultsFile.open("FHN.csv"); //change filename to something more descriptive

	// Print the initial conditions
	cout << "  " << setw(8) << t;
	resultsFile << t;
	for (int i = 0; i < neqn; i++) {
		cout << "  " << setw(14) << y[i];
		resultsFile << ",";
		resultsFile << y[i];
	}
	cout << "\n";
	resultsFile << endl;

	work = new double[100 + 21 * neqn]; //want to work out what this does and why it uses 100 and 21



	for (i = 1; i <= step_num; i++)
	{
		tout = (double)(i)* tf * stepSize;

		ode(f01, neqn, y, t, tout, relerr, abserr, iflag, work, iwork);

		if (iflag != 2)
		{
			cout << "\n";
			cout << "TEST01 - Fatal error!\n";
			cout << "  ODE returned IFLAG = " << iflag << "\n";
			break;
		}

		// Print the output for to console and to csv
		cout << "  " << setw(8) << t;
		resultsFile << t;
		for (int i = 0; i < neqn; i++) {
			cout << "  " << setw(14) << y[i];
			resultsFile << ",";
			resultsFile << y[i];
		}
		cout << "\n";
		resultsFile << endl;
	}

	system("pause");
	delete[] work;
	delete[] y;

	return;
}
//****************************************************************************80

//****************************************************************************80


double uni_reverse(double I_ext, double k, double y[], int i, int N, int loop_number) {
	if (i != 0) { return k * (y[i - 2] - y[i]); }
	else { return I_ext; }
}

double bi_nearest_neighbour(double I_ext, double k, double y[], int i, int Nneurons, int loop_number) {
	//1D chain of neurones with free standing boundary
	if (i == 0) { return I_ext + k * (y[2] - y[0]); }
	else if (i == Nneurons) { return k * (y[Nneurons - 2] - y[Nneurons]); }
	else { return k * (y[i + 2] - y[i] + y[i] - y[i - 2]); }
}

double bi_nearest_neighbour_circular(double I_ext, double k, double y[], int i, int Nneurons, int loop_number) {

	// I_ext is only input to the first neurone of the first cycle
	// Removing this statement gives periodic driving
	if (loop_number != 0) { I_ext = 0; }

	//1D chain of neurones in a ring
	if (i == 0) { return I_ext + k * (y[2] - y[0]) + k * (y[0] - y[Nneurons]); }
	else if (i == Nneurons) { return k * (y[Nneurons - 2] - y[Nneurons]) + k * (y[Nneurons] - y[0]); }
	else { return k * (y[i + 2] - y[i] + y[i] - y[i - 2]); }
}

double random_driving_bi_nearest_neighbour(double k, double y[], int i, int Nneurons, double t) {

	//1D chain of neurones in a ring
	if (i == 0) { return random_driver(t) + k * (y[2] - y[0]) + k * (y[0] - y[Nneurons]); }
	else if (i == Nneurons) { return random_driver(t) + k * (y[Nneurons - 2] - y[Nneurons]) + k * (y[Nneurons] - y[0]); }
	else { return random_driver(t) + k * (y[i + 2] - y[i] + y[i] - y[i - 2]); }

}

double random_driver(double t) {
	srand((int)time(0));
	double I = -(rand() % 70) / 100;
	return I;
}

double bi_direction_2D_lattice(double I_ext, double k, double y[], int i, int Nneurons, int loop_number) {
	//Rectangular 2D lattice with periodic boundary conditions
	
	// I_ext is only input to the first neurone of the first cycle
	// Removing this statement gives periodic driving
	if (loop_number != 0) { I_ext = -0.58; }
	//Number of columns, must be compatible with Nneurons
	int Nc = sqrt(Nneurons);

	//Edge cases first

	//Left Vertical edge (non-corner)
	if (i % (2*Nc) == 0 && i != 0 && i != (neqn-2*Nc) ) {
		return k * (y[i + 2] - y[i + (2*Nc) - 2] + y[i+2*Nc] - y[i-2*Nc]);
	}

	//Right vertical edge (non-corner)
	else if (i % (2*Nc) == 2*Nc - 2 && i != 2*Nc-2 && i != neqn-2) {
		return k * (y[i-2*Nc+2]-y[i-2] + y[i + 2 * Nc] - y[i - 2 * Nc]);
	}

	//Top horizontal
	else if (0<i && i<2*Nc-2) {
		return k * (y[i + 2] - y[i - 2] + y[i+2*Nc] - y[neqn-2*Nc+i]);
	}

	//Bottom horizontal
	else if (neqn-2*Nc<i && i<neqn-2) {
		return k * (y[i + 2] - y[i - 2] - y[i-2*Nc] + y[i-(neqn-2*Nc)]);
	}

	//Top left corner
	else if (i==0) {
		return I_ext + k * (y[i+2]-y[2*Nc-2]+y[2*Nc]-y[neqn-2*Nc]);
	}

	//Top right corner
	else if (i == 2 * Nc - 2) {
		return k * (-y[i-2] + y[0] + y[i+2*Nc] - y[neqn-2]);
	}

	//Bottom left corner
	else if (i == neqn - 2 * Nc) {
		return k * (y[i + 2] - y[neqn - 2] + y[0] - y[i - 2 * Nc]);
	}

	//Bottom right corner
	else if (i == neqn - 2) {
		return k * (-y[i - 2] + y[neqn - 2 * Nc] + y[2 * Nc - 2] - y[i - 2 * Nc]);
	}
	
	//middle bits
	else{
		return k * (y[i + 2] - y[i - 2] + y[i + 2 * Nc] - y[i - 2 * Nc]);
	}
}

void f01(double t, double y[], double yp[]) {
	const double a = 0.75;
	const double b = 0.8;
	const double c = 3.0;
	const double c_inv = 1 / c;
	const double I_ext = -0.58;
	const double k = 0.113;

	int loop_number = 0;
	int Nloops = 1;

	for (loop_number = 0; loop_number < Nloops; loop_number++) {
		// Each iteration is one cycle around the ring of neurones
		for (int i = 0; i < neqn; i += 2) {
			yp[i] = c * (y[i] - 0.33333333*y[i] * y[i] * y[i] + y[i + 1] + bi_direction_2D_lattice(I_ext,k,y,i,Nneurons,loop_number));
			yp[i + 1] = -c_inv * (y[i] - a + b * y[i + 1]);
		}
	}

	return;
}