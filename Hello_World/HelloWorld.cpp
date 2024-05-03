#include <iostream>
#include <omp.h>
#include <stdlib.h>
using namespace std;
void hello();
int main(int argc, char *argv[])
{
    int thread_count = (argc > 1 ? atoi(argv[1]) : 8);
#pragma omp parallel num_threads(thread_count)
    {
        hello();
    }
}
void hello()
{
    int my_rank = omp_get_thread_num();
    int thread_count = omp_get_num_threads();
#pragma omp critical
    cout << "Hello from thread " << my_rank << " of " << thread_count << endl;
};
