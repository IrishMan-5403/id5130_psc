
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <omp.h>

using namespace std;

double calcSin(double value) {
	return sin(5 * value);
}

int main(int argc, char* argv[]) {

	int num_threads = (argc > 1 ? atoi(argv[1]) : 8);

	double lower_limit = 0, upper_limit = 3;
	int n = 1000;
	double h = (upper_limit - lower_limit) / (n - 1);

	double l_diag[n], m_diag[n], u_diag[n], results[n], solutions[n];
	for(int index = 0; index < n; index++) {
		l_diag[index] = 1;
		m_diag[index] = 4;
		u_diag[index] = 1;
		results[index] = 3 / h * (calcSin(h * (index + 1)) - calcSin(h * (index - 1)));
	}
	l_diag[0] = 0; l_diag[n-1] = 2;
	m_diag[0] = 1; m_diag[n-1] = 1;
	u_diag[0] = 2; u_diag[n-1] = 0;
	results[0]  = 1 / h * (-5.0 / 2 * calcSin(lower_limit) + 2 * calcSin(lower_limit + h) + 1.0 / 2 * calcSin(lower_limit + 2 * h));
	results[n-1] = 1 / h * (5.0 / 2 * calcSin(upper_limit) - 2 * calcSin(upper_limit - h) - 1.0 / 2 * calcSin(upper_limit - 2 * h));

	double alpha, beta;
	int limit = ceil(log2(n));

	double computation_start = omp_get_wtime();
	int idx;
	#pragma omp parallel  num_threads(num_threads) private(idx, alpha, beta)
	for(int i = 0; i < limit; i++) {
		double temp_lower[n+1], temp_main[n+1], temp_upper[n+1], temp_results[n+1];
		int stride = pow(2, i), idx;

		#pragma omp for
		for(idx = 0; idx < n; idx++) {
			temp_main[idx] = m_diag[idx];
			temp_results[idx] = results[idx];
			if(idx >= stride) {
				alpha = -l_diag[idx] / m_diag[idx - stride];
				temp_lower[idx] = alpha * l_diag[idx - stride];
				temp_main[idx] += alpha * u_diag[idx - stride];
				temp_results[idx] += alpha * results[idx - stride];
			}
			if(idx < n - stride) {
				beta = -u_diag[idx] / m_diag[idx + stride];
				temp_upper[idx] = beta * u_diag[idx + stride];
				temp_main[idx] += beta * l_diag[idx + stride];
				temp_results[idx] += beta * results[idx + stride];
			}
		}

		#pragma omp for
		for( idx = 0; idx < n; idx++) {
			if(idx >= stride) {
				l_diag[idx] = temp_lower[idx];
			}
			if(idx < n - stride) {
				u_diag[idx] = temp_upper[idx];
			}
			m_diag[idx] = temp_main[idx];
			results[idx] = temp_results[idx];
		}
	}

	double computation_end = omp_get_wtime();
	#pragma omp for
	for(int idx = 0; idx < n; idx++) {
		solutions[idx] = results[idx] / m_diag[idx];
	}



	cout << "Computation Time: " << computation_end - computation_start << endl;
}




// Computation Time
//p=2  t=0.000147089
//p=4  t=0.00013314
//p=8  t=0.000200241