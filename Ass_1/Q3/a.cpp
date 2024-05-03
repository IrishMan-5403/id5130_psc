#include<bits/stdc++.h>
#include <omp.h>

using namespace std;

void init_grid(double x_grid [], double y_grid [], double delta, int N, double x_min, double x_max, double y_min, double y_max){
    x_grid[0] = x_min;
    y_grid[0] = y_min;

    for(int i = 1; i < N+1; i++){
        x_grid[i] = x_grid[i-1] + delta;
        y_grid[i] = y_grid[i-1] + delta;
    }

    return;
}

double exact_soln(double x, double y){
    return ((pow(x,2)-1)*(pow(y,2)-1));
}

double q_func(double x, double y){
    return (2*(2-pow(x,2)-pow(y,2)));
}

int main(int argc, char* argv[])
{
  
  int thread_count = (argc > 1 ? atoi(argv[1]) : 8);
  

  double delta = 0.01;


  double xmin = -1.000;
  double xmax = 1.000;
  double ymin = -1.000;
  double ymax = 1.000;

  int N = int((xmax - xmin)/delta);

  double tolerance = 1e-2;
  double g_error = 1.00;

  double x_grid [N+1];
  double y_grid [N+1];
  double phi    [N+1][N+1];
  double q      [N+1][N+1];
  double phi_exact [N+1][N+1];
  double error [N+1][N+1];

 

  init_grid(x_grid,y_grid,delta,N,xmin,xmax,ymin,ymax);

  for(int i = 0; i < N+1; i++){
    for(int j = 0; j < N+1; j++){
        phi_exact[i][j] = exact_soln(x_grid[i],y_grid[j]);
        q[i][j] = q_func(x_grid[i], y_grid[j]);
        phi[i][j] = 0.00;
    }
  }

  for(int i = 0; i < N+1; i++){
    phi[i][0] = 0.00;
    phi[i][N] = 0.00;
  }

  for(int j = 0; j < N+1; j++){
    phi[0][j] = 0.00;
    phi[N][j] = 0.00;
  }

  double start_time = omp_get_wtime();

  int itr_count = 0;

  while(g_error >= tolerance){

    for(int i = 1; i < N; i++){
        for(int j = 1; j < N; j++){
            phi[i][j] = (phi[i+1][j]+phi[i-1][j]+phi[i][j+1]+phi[i][j-1] + q[i][j]*pow(delta,2))/4.00;
            error[i][j] = fabs((phi_exact[i][j] - phi[i][j])/phi_exact[i][j]);
        }
    }

    g_error=error[1][1];
    for(int i = 1; i < N; i++){
        for(int j = 1; j < N; j++){
            if(g_error < error[i][j]){
                g_error = error[i][j];
            }
        }
    }

    itr_count++;
    // if(itr_count % 10== 0){
    //     cout << endl << "Iteration Count : " << itr_count << " Error : " << global_error << endl;
    // }
  }

  double end_time = omp_get_wtime();
  

//   double y_exact[N+1];
//   double y[N+1];
//   for (int i=0;i<N+1;i++){
//      y[i]=phi[i][5];
//      y_exact[i]=phi_exact[i][5];

//   }

  cout << "Time Taken : " << end_time - start_time  << endl;
//   cout<<"y - ";
//     for (int i=0;i<N+1;i++){
//      cout<<y[i]<<",";
//     }
//     cout<<endl;
//     cout<<"y_exact - ";
//       for (int i=0;i<N+1;i++){
    
//        cout<<y_exact[i]<<",";

//   }
//     cout<<endl;

  return 0;
}

