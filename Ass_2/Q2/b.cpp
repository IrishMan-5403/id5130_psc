#include <bits/stdc++.h>
#include <mpi.h>

// mpic++ -o b b.cpp
// mpirun -np 4 ./b

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
    int rank, size;
    double *all_error = NULL;
    int I, J;
    MPI_Init(NULL, NULL);

    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    double delta = 0.1;

    double xmin = -1.000;
    double xmax = 1.000;
    double ymin = -1.000;
    double ymax = 1.000;

    int N = int((xmax - xmin) / delta);
    int local_N = N / size;

    double tolerance = 0.0001;
    double g_error = 1.00;

    double x_grid[N + 1];
    double y_grid[N + 1];
    double phi1[N + 1][N + 1];
    double phi2[N + 1][N + 1];
    double q[N + 1][N + 1];
    double error[N + 1][N + 1];

    if (rank == 0)
    {
        init_grid(x_grid, y_grid, delta, N, xmin, xmax, ymin, ymax);
        all_error = (double *)malloc(size * sizeof(double));
    }

    MPI_Bcast(x_grid, N + 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);

    int start_index = rank * local_N;
    int end_index = (rank + 1) * local_N;
    if (rank == size - 1)
    {
        end_index = N;
    }

    for (int i = start_index; i < end_index + 1; i++)
    {
        for (int j = 0; j < N + 1; j++)
        {
            q[i][j] = q_func(x_grid[j], y_grid[i]);
            phi1[i][j] = 0.00;
            phi2[i][j] = 0.00;
        }
    }

    for (int i = start_index; i < end_index + 1; i++)
    {
        phi1[i][0] = sin(2.0 * pi * i);
        phi2[i][0] = sin(2.0 * pi * i);
    }

    int itr_count = 0;

    while (g_error >= tolerance)
    {

        if (itr_count % 2 == 0)
        {
            for (int i = start_index; i < end_index + 1; i++)
            {
                if (i == 0)
                    continue;
                for (int j = 1; j < N; j++)
                {
                    phi1[i][j] = (phi2[i + 1][j] + phi2[i - 1][j] + phi2[i][j + 1] + phi2[i][j - 1] + q[i][j] * pow(delta, 2)) / 4.00;
                    error[i][j] = fabs((phi1[i][j] - phi2[i][j]));
                }
            }
            for (int i = start_index; i < end_index + 1; i++)
            {
                phi1[i][N] = (4.0 * phi1[i][N - 1] - phi1[i][N - 2]) / 3.0;
                error[i][N] = fabs(phi1[i][N] - phi2[i][N]);
            }
        }
        else
        {
            for (int i = start_index; i < end_index + 1; i++)
            {
                for (int j = 1; j < N; j++)
                {
                    phi2[i][j] = (phi1[i + 1][j] + phi1[i - 1][j] + phi1[i][j + 1] + phi1[i][j - 1] + q[i][j] * pow(delta, 2)) / 4.00;
                    error[i][j] = fabs((phi1[i][j] - phi2[i][j]));
                }
            }
            for (int i = 1; i < N; i++)
            {
                phi2[i][N] = (4.0 * phi2[i][N - 1] - phi2[i][N - 2]) / 3.0;
                error[i][N] = fabs(phi1[i][N] - phi2[i][N]);
            }
        }

        g_error = error[start_index + 1][1];
        for (int i = start_index + 1; i < end_index; i++)
        {
            for (int j = 1; j < N; j++)
            {
                if (g_error < error[i][j])
                {
                    g_error = error[i][j];
                    I = i;
                    J = j;
                }
            }
        }
        // cout << "(" << I << ","
        //      << J
        //      << ")" << endl;

        itr_count++;
        MPI_Gather(&g_error, 1, MPI_DOUBLE, all_error, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);

        if (rank == 0)
        {
            for (int i = 0; i < size; i++)
            {
                if (all_error[i] > g_error)
                {
                    g_error = all_error[i];
                }
            }

            // if (itr_count == 1)
            cout << "Iteration Count : " << itr_count << " Error : " << g_error << " " << I << " " << J << endl;
            // break;
        }
        MPI_Bcast(&g_error, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);

        MPI_Barrier(MPI_COMM_WORLD);
    }

    //   double y_exact[N+1];
    //   double y[N+1];
    //   for (int i=0;i<N+1;i++){
    //      y[i]=phi[i][5];
    //      y_exact[i]=phi_exact[i][5];

    //   }
    if (rank == 0)
        cout << "Iteration Count : " << itr_count << " Error : " << g_error << endl;

    // cout << "phi (y=0) - ";
    // for (int i = 0; i < N + 1; i++)
    // {
    //     if (itr_count % 2 == 0)
    //         cout << phi1[i][10] << ",";

    //     else
    //     {
    //         cout << phi2[i][10] << ",";
    //     }
    // }
    // cout << endl;
    // cout << "phi (x=0) - ";
    // for (int i = 0; i < N + 1; i++)
    // {
    //     if (itr_count % 2 == 0)
    //         cout << phi1[10][i] << ",";

    //     else
    //     {
    //         cout << phi2[10][i] << ",";
    //     }
    // }
    // cout << endl;

    MPI_Finalize();

    return 0;
}
