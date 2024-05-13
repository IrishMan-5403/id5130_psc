#include <time.h>
#include <math.h>
#include <stdio.h>

double function(double x)
{
    return sin(5 * x);
}

int main()
{

    double lower_limit = 0, upper_limit = 3;
    int num_intervals = 1000;
    double step_size = (upper_limit - lower_limit) / (num_intervals - 1);

    double a[num_intervals], b[num_intervals], c[num_intervals], y[num_intervals], l[num_intervals], u[num_intervals], x[num_intervals];

    clock_t timer;
    timer = clock();

// Initialization
#pragma acc parallel loop num_gangs(1000) create(a[0 : num_intervals], b[0 : num_intervals], c[0 : num_intervals], y[0 : num_intervals])
    for (int i = 0; i < num_intervals; i++)
    {
        a[i] = 1;
        b[i] = 4;
        c[i] = 1;
        y[i] = 3 / step_size * (function(step_size * (i + 1)) - function(step_size * (i - 1)));
    }
    a[0] = 0;
    a[num_intervals - 1] = 2;
    b[0] = 1;
    b[num_intervals - 1] = 1;
    c[0] = 2;
    c[num_intervals - 1] = 0;
    y[0] = 1 / step_size * (-5.0 / 2 * function(lower_limit) + 2 * function(lower_limit + step_size) + 1.0 / 2 * function(lower_limit + 2 * step_size));
    y[num_intervals - 1] = 1 / step_size * (5.0 / 2 * function(upper_limit) - 2 * function(upper_limit - step_size) - 1.0 / 2 * function(upper_limit - 2 * step_size));

    // Forward Elimination
    u[0] = b[0];
#pragma acc parallel loop num_gangs(1000) present(a[0 : num_intervals], b[0 : num_intervals], c[0 : num_intervals], y[0 : num_intervals]) copyout(l[0 : num_intervals], u[0 : num_intervals], y[0 : num_intervals])
    for (int i = 1; i < num_intervals; i++)
    {
        l[i] = a[i] / u[i - 1];
        u[i] = b[i] - l[i] * c[i - 1];
        y[i] = y[i] - l[i] * y[i - 1];
    }

    // Backward Substitution
    x[num_intervals - 1] = y[num_intervals - 1] / u[num_intervals - 1];
#pragma acc kernels
    for (int i = num_intervals - 2; i >= 0; i--)
    {
        x[i] = (y[i] - c[i] * x[i + 1]) / u[i];
    }

    timer = clock() - timer;
    printf("Time taken: %fs\n", ((float)timer) / CLOCKS_PER_SEC);

    // Print the solution
    printf("[");
    for (int i = 0; i < num_intervals; i++)
    {
        printf("%f,", x[i]);
    }
    printf("]\n");

    return 0;
}
