#include <cmath>
#include <vector>

#include "odeint.hpp"


typedef std::vector< double > state_type;

//define constants

const double a = 0.7;
const double b = 0.8;
const double c = 12.5;
const double gamma = 1;
const int N = 20;

void neuron_network(const state_type &x, state_type &dxdt, const double /* t */) {
	int i;
	for (i = 0; i < N; i++) {
		if (i % 2 == 0 && i < N - 1) {
			dxdt[i] = c * (x[i] - ((1 / 3)*(x[i] - pow(x[i], 3))) + x[i + 1]);
		}
		else{
			dxdt[i] = (-1 / c)*(x[i - 1] - a + b * x[i]);
		}
	}
}



struct push_back_state_and_time
{
	std::vector< state_type >& m_states;
	std::vector< double >& m_times;

	push_back_state_and_time(std::vector< state_type > &states, std::vector< double > &times)
		: m_states(states), m_times(times) { }

	void operator()(const state_type &x, double t)
	{
		m_states.push_back(x);
		m_times.push_back(t);
	}
};

int main() {

	using namespace std;
	using namespace boost::numeric::odeint;

	vector<state_type> x_vec;
	vector<double> times;

	state_type x(20); //make sure number in the brackets = N
	int j;
	for (j = 0; j < N; j++) {
		x[j] = 0;
	}

	size_t steps = integrate(neuron_network,
		x, 0.0, 10.0, 0.1,
		push_back_state_and_time(x_vec, times)); //variable step RK4 5th order


	return 0;
}