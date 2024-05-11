#include <iostream>
#include <omp.h>
#include <stdlib.h>
using namespace std;


int main(int argc, char *argv[])
{
    int thread_count = (argc > 1 ? atoi(argv[1]) : 8);


#pragma omp parallel num_threads(thread_count)
    {
        #pragma for 
    }

};