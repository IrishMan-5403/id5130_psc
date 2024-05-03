#include<iostream>
#include <bits/stdc++.h>
#include <omp.h>

using namespace std;

int main(int argc , char*argv[])
{
    int threads = ( argc > 1 ? atoi (argv[1]): 8);
    vector<int> arr = { 7,3,6,9,10,11,23,434,545,656,787,878,323};
    int size=arr.size();
    int res =33333333;
    #pragma omp parallel for num_threads (threads) reduction (/:res)
    
        for (int i=0;i<size;i++)
        {
            res/=arr[i];
        }
    
    cout<<res<<endl;


}