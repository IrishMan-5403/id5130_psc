#include <bits/stdc++.h>
#include <mpi.h>

using namespace std;

const double pi = 3.14159265358979323846;

void init_grid(double x_grid[], double delta, int N, double x_min, double x_max)
{
    x_grid[0] = x_min;

    for (int i = 1; i < N + 1; i++)
    {
        x_grid[i] = x_grid[i - 1] + delta;
    }

    return;
}

int main(int argc, char **argv)
{
    double dx = 0.002;
    double dt = 0.0001;

    double xmin = 0;
    double xmax = 2.0; // L = 2.0
    double T = 1.0;

    int N = int((xmax - xmin) / dx);
    int Nt = int((T - 0.0) / dt);

    int rank, size;
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    int ln = ceil(double(N - 1) / size);
    int ln_max = ln;
    int la = rank * ln;
    int lb = la + ln;
    lb = lb > N - 1 ? N - 1 : lb;
    ln = lb - la;

    double *u, *nu;
    u = (double *)malloc(ln_max * size * sizeof(double));
    nu = (double *)malloc(ln_max * sizeof(double));

    double u_grid_upw[N + 1 + 1];
    double u_grid_upw2[ln_max + 1];
    double u_grid_quick[N + 1 + 1];
    double u_grid_quick2[ln_max + 1];

    if (rank == 0)
    {
        double x_grid[N + 1];
        init_grid(x_grid, dx, N, xmin, xmax);

        for (int i = 0; i < N / 4; ++i)
        {
            u_grid_quick[i] = sin(4 * pi * (x_grid[i]));
            u_grid_upw[i] = sin(4 * pi * (x_grid[i]));
        }

        for (int i = N / 4; i < N + 1; i++)
        {
            u_grid_quick[i] = 0.0;
            u_grid_upw[i] = 0.0;
        }
    }
    MPI_Bcast(u_grid_quick, N + 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast(u_grid_upw, N + 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);

    for (int j = 0; j < Nt + 1; j++)
    {

        for (int i = 0; i < ln; ++i)
        {
            u_grid_upw2[i] = u_grid_upw[la + i] - dt * ((u_grid_upw[la + i] - u_grid_upw[la + i - 1]) / dx);
        }

        if (rank == 0)
        {
            u_grid_quick2[1] = u_grid_quick[1] - dt * (u_grid_quick[1] - u_grid_quick[0]) / dx;
        }
        else
        {
            u_grid_quick2[0] = u_grid_quick[la] - dt * ((3.0 * u_grid_quick[la] - 7.0 * u_grid_quick[la - 1] + u_grid_quick[la - 2] + 3.0 * u_grid_quick[la + 1]) / (8.0 * dx));
            u_grid_quick2[1] = u_grid_quick[la + 1] - dt * ((3.0 * u_grid_quick[la + 1] - 7.0 * u_grid_quick[la] + u_grid_quick[la - 1] + 3.0 * u_grid_quick[la + 2]) / (8.0 * dx));
        }
        for (int i = 2; i < N; ++i)
        {
            u_grid_quick2[i] = u_grid_quick[i] - dt * ((3.0 * u_grid_quick[i] - 7.0 * u_grid_quick[i - 1] + u_grid_quick[i - 2] + 3.0 * u_grid_quick[i + 1]) / (8.0 * dx));
        }

        MPI_Allgather(u_grid_quick2, ln_max, MPI_DOUBLE, u_grid_quick, ln_max, MPI_DOUBLE, MPI_COMM_WORLD);
        MPI_Allgather(u_grid_upw2, ln_max, MPI_DOUBLE, u_grid_upw, ln_max, MPI_DOUBLE, MPI_COMM_WORLD);

        if (j == (Nt / 2) && rank == 0)
        {
            cout << "t=0.5" << endl;
            cout << "upwind - [";
            for (int i = 0; i < N + 1; i++)
            {
                cout << u_grid_upw[i] << ",";
            }
            cout << "]" << endl;
            cout << "quick - [";
            for (int i = 0; i < N + 1; i++)
            {
                cout << u_grid_quick[i] << ",";
            }
            cout << "]" << endl;
        }
        MPI_Barrier(MPI_COMM_WORLD);
    }

    if (rank == 0)
    {

        cout << "t=1.0" << endl;
        cout << "upwind - [";
        for (int i = 0; i < N + 1; i++)
        {
            cout << u_grid_upw[i] << ",";
        }
        cout << "]" << endl;

        cout << "quick - [";
        for (int i = 0; i < N + 1; i++)
        {
            cout << u_grid_quick[i] << ",";
        }
        cout << "]" << endl;
    }

    MPI_Finalize();

    return 0;
}
