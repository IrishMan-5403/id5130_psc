#include <iostream>
#include <omp.h>
#include <stdlib.h>
using namespace std;


int global_var = 0; // Global variable shared among threads

int main(int argc , char*argv[]) {
    int threads = ( argc > 1 ? atoi (argv[1]): 8);

    int private_var = 0; // Variable private to each thread
    int firstprivate_var = 0; // Variable private, but initialized once by the master thread
    int threadprivate_var; // Variable private to each thread, initialized separately for each thread

    #pragma omp parallel  num_threads (threads) private(private_var) firstprivate(firstprivate_var) 
    {
        int thread_id = omp_get_thread_num();

        private_var += 10 + thread_id;
        firstprivate_var += 20;

        #pragma omp critical
        {
            global_var += 1;
            threadprivate_var = 30 + thread_id; // Initialize threadprivate variable
            cout<<global_var<<"  "<<threadprivate_var<<endl;
        }

        // Display the values of each variable
        #pragma omp critical
        {
            cout << "Thread " << thread_id << ": private_var = " << private_var
                      << ", firstprivate_var = " << firstprivate_var
                      << ", threadprivate_var = " << threadprivate_var
                      << ", global_var = " << global_var << std::endl;
        }
    }

    // Display the final values of variables after the parallel region
    cout << "After parallel region: global_var = " << global_var << std::endl;

    return 0;
}


//    private and firstprivate create a private copy for each thread, but firstprivate 
//    initializes the private copies with the initial value from the master thread.
   
//    threadprivate allows each thread to have its own private copy with a potentially 
//    different initial value for the same variable.