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

    double delta = 0.1;

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

    int itr_count = 0;

    while (g_error >= tolerance)
    {

        if (itr_count % 2 == 0)
        {
            for (int i = 1; i < N; i++)
            {
                for (int j = 1; j < N; j++)
                {
                    phi1[i][j] = (phi2[i + 1][j] + phi2[i - 1][j] + phi2[i][j + 1] + phi2[i][j - 1] + q[i][j] * pow(delta, 2)) / 4.00;
                    error[i][j] = fabs((phi1[i][j] - phi2[i][j]));
                }
                phi1[i][N] = (4.0 * phi1[i][N - 1] - phi1[i][N - 2]) / 3.0;
                error[i][N] = fabs((phi1[i][N] - phi2[i][N]));
            }
        }
        else
        {
            for (int i = 1; i < N; i++)
            {
                for (int j = 1; j < N; j++)
                {
                    phi2[i][j] = (phi1[i + 1][j] + phi1[i - 1][j] + phi1[i][j + 1] + phi1[i][j - 1] + q[i][j] * pow(delta, 2)) / 4.00;
                    error[i][j] = fabs((phi1[i][j] - phi2[i][j]));
                }
                phi2[i][N] = (4.0 * phi2[i][N - 1] - phi2[i][N - 2]) / 3.0;
                error[i][N] = fabs((phi1[i][N] - phi2[i][N]));
            }
        }

        g_error = error[1][1];
        for (int i = 1; i < N; i++)
        {
            for (int j = 1; j < N + 1; j++)
            {
                if (g_error < error[i][j])
                {
                    g_error = error[i][j];
                }
            }
        }

        itr_count++;
        if (itr_count % 10 == 0)
        {
            cout << endl
                 << "Iteration Count : " << itr_count << " Error : " << g_error << endl;
        }
    }

    cout << endl
         << "Iteration Count : " << itr_count << " Error : " << g_error << endl;

    cout << "phi (y=0) - [";
    for (int i = 0; i < N + 1; i++)
    {
        if (itr_count % 2 == 0)
            cout << phi1[10][i] << ",";

        else
        {
            cout << phi2[10][i] << ",";
        }
    }
    cout << "]" << endl;
    cout << "phi (x=0) - [";
    for (int i = 0; i < N + 1; i++)
    {
        if (itr_count % 2 == 0)
            cout << phi1[i][10] << ",";

        else
        {
            cout << phi2[i][10] << ",";
        }
    }
    cout << "]" << endl;

    return 0;
}

// Iteration Count : 387 Error : 9.94133e-05
// phi (y=0) - [0,0.0329707,0.0577094,0.0759284,0.0892349,0.0990919,0.106791,0.113436,0.119934,0.126992,0.135121,0.144632,0.155642,0.16807,0.181636,0.195856,0.210035,0.223261,0.23439,0.242042,0.244593,]
// phi (x=0) - [0,0.0389741,0.0688925,0.0908486,0.106163,0.11635,0.122966,0.127398,0.130654,0.133233,0.135121,0.135907,0.134981,0.131725,0.12564,0.11635,0.10349,0.0865222,0.064566,0.0363002,0,]
// [-1,-0.9,-0.8,-0.7,-0.6,-0.5,-0.4,-0.3,-0.2,-0.1,0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1]