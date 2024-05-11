#include <cmath>
#include <iostream>

using namespace std;

double calcSin(double x)
{
	return sin(5 * x);
}

int main()
{

	double lower = 0, upper = 3;
	int n = 25;
	double h = (upper - lower) / (n - 1);

	double a[n], b[n], c[n], y[n], l[n], u[n], x[n];
	for (int i = 0; i < n; i++)
	{
		a[i] = 1;
		b[i] = 4;
		c[i] = 1;
		y[i] = 3 / h * (calcSin(h * (i + 1)) - calcSin(h * (i - 1)));
	}
	a[0] = 0;
	a[n - 1] = 2;
	b[0] = 1;
	b[n - 1] = 1;
	c[0] = 2;
	c[n - 1] = 0;
	y[0] = 1 / h * (-5.0 / 2 * calcSin(lower) + 2 * calcSin(lower + h) + 1.0 / 2 * calcSin(lower + 2 * h));
	y[n - 1] = 1 / h * (5.0 / 2 * calcSin(upper) - 2 * calcSin(upper - h) - 1.0 / 2 * calcSin(upper - 2 * h));

	u[0] = b[0];
	for (int i = 1; i < n; i++)
	{
		l[i] = a[i] / u[i - 1];
		u[i] = b[i] - l[i] * c[i - 1];
		y[i] = y[i] - l[i] * y[i - 1];
	}

	x[n - 1] = y[n - 1] / u[n - 1];
	for (int i = n - 2; i >= 0; i--)
	{
		x[i] = (y[i] - c[i] * x[i + 1]) / u[i];
	}

	cout << "[";
	for (int i = 0; i < n - 1; i++)
	{
		cout << x[i] << ",";
	}
	cout << x[n - 1] << "]" << endl;
}

//[5.12378,4.01686,1.58442,-1.4988,-4.0015,-4.99505,-4.09911,-1.65366,1.41705,3.952,4.99281,4.14596,1.73164,-1.33738,-3.90076,-4.98937,-4.19164,-1.80912,1.25725,3.84874,4.98346,4.24024,1.87087,-1.11988,-4.00803]