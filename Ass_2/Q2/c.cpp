#include <bits/stdc++.h>

using namespace std;
const double pi = 3.14159265358979323846;

void init_grid(double x_grid[], double y_grid[], double delta, int N, double x_min, double x_max, double y_min, double y_max)
{
    x_grid[0] = x_min;
    y_grid[0] = y_min;

    for (int i = 1; i < N + 1; i++)
    {
        x_grid[i] = x_grid[i - 1] + delta;
        y_grid[i] = y_grid[i - 1] + delta;
    }
    x_grid[10] = 0;
    y_grid[10] = 0;

    return;
}

double q_func(double x, double y)
{
    return ((pow(x, 2) + pow(y, 2)));
}

int main(int argc, char *argv[])
{

    double delta = 0.01;

    double xmin = -1.000;
    double xmax = 1.000;
    double ymin = -1.000;
    double ymax = 1.000;

    int N = int((xmax - xmin) / delta);

    double tolerance = 0.0001;
    double g_error = 1.00;

    double x_grid[N + 1];
    double y_grid[N + 1];
    double phi1[N + 1][N + 1];
    double phi2[N + 1][N + 1];
    double q[N + 1][N + 1];
    double error[N + 1][N + 1];

    init_grid(x_grid, y_grid, delta, N, xmin, xmax, ymin, ymax);

    for (int i = 0; i < N + 1; i++)
    {
        for (int j = 0; j < N + 1; j++)
        {
            q[i][j] = q_func(x_grid[i], y_grid[j]);
            phi1[i][j] = 0.00;
            phi2[i][j] = 0.00;
        }
    }
    for (int i = 0; i < N + 1; i++)
    {
        phi1[i][0] = sin(2.0 * pi * y_grid[i]);
        phi2[i][0] = sin(2.0 * pi * y_grid[i]);
    }

    int iteration_count = 0;
    int i, j;

    while (g_error >= tolerance)
    {

        for (i = 1; i < N; i++)
        {
            for (j = 1; j < N; j++)
            {
                if ((i + j) % 2 == 1)
                {
                    phi2[i][j] = (phi1[i + 1][j] + phi1[i - 1][j] + phi1[i][j + 1] + phi1[i][j - 1] + q[i][j] * pow(delta, 2)) / 4.00;
                    error[i][j] = fabs(phi1[i][j] - phi2[i][j]);
                    phi1[i][j] = phi2[i][j];
                }
            }
        }

        for (i = 1; i < N; i++)
        {
            for (j = 1; j < N; j++)
            {
                if ((i + j) % 2 == 0)
                {
                    phi2[i][j] = (phi1[i + 1][j] + phi1[i - 1][j] + phi1[i][j + 1] + phi1[i][j - 1] + q[i][j] * pow(delta, 2)) / 4.00;
                    error[i][j] = fabs(phi1[i][j] - phi2[i][j]);
                    phi1[i][j] = phi2[i][j];
                }
            }
        }

        iteration_count++;

        g_error = error[1][1];

        for (i = 1; i < N; i++)
        {
            for (j = 1; j < N; j++)
            {
                if (g_error < error[i][j])
                {
                    g_error = error[i][j];
                }
            }
        }
        if (iteration_count % 10 == 0)
        {
            cout << endl
                 << "Iteration Count : " << iteration_count << " Error : " << g_error << endl;
        }
    }
    cout << "Iteration Count : " << iteration_count << " Error : " << g_error << endl;

    cout << "phi (y=0) - [";
    for (int i = 0; i < N + 1; i++)
    {
        cout << phi1[100][i] << ",";
    }
    cout << "]" << endl;

    return 0;
}