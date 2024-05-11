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
#pragma acc enter data create(a[ : n], b[ : n], c[ : n], y[ : n], l[ : n], u[ : n], x[ : n])

// Initialize arrays a, b, c, and calculate y in parallel
#pragma acc parallel loop present(a, b, c, y)
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

    // Forward elimination and backward substitution
    u[0] = b[0];
#pragma acc parallel loop present(l, u, y)
    for (int i = 1; i < n; i++)
    {
        l[i] = a[i] / u[i - 1];
        u[i] = b[i] - l[i] * c[i - 1];
        y[i] = y[i] - l[i] * y[i - 1];
    }

    x[n - 1] = y[n - 1] / u[n - 1];
#pragma acc parallel loop present(x, y, c, u)
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

#pragma acc exit data delete (a[ : n], b[ : n], c[ : n], y[ : n], l[ : n], u[ : n], x[ : n])

    return 0;
}
