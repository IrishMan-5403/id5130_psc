#include<iostream>
#include <bits/stdc++.h>
#include <omp.h>

using namespace std;

int main(int argc , char*argv[])
{
    int threads = ( argc > 1 ? atoi (argv[1]): 8);
    vector<int> arr = { 41984,444,78,525,11984,15164,5251,65,49,84,8,49849,454,8484,44};
    int size=arr.size();
    int tmp;
    int phase;
    int i;
    #pragma omp parallel num_threads (threads) default(none) shared (arr,size) private (i,tmp,phase)
        
        for (phase=0;phase<size;phase++)
        {
            if (phase % 2 == 0){
                #pragma omp for
                for ( i = 1; i < size; i+=2)
                {
                    if (arr[i-1]>arr[i])
                    {
                        tmp = arr[i-1];
                        arr[i-1]=arr[i];
                        arr[i]=tmp;
                    }
                }                
            }
            else{
                #pragma omp for 
                for (i=1;i<size;i+=2)
                {
                    if (arr[i]>arr[i+1])
                    {
                        tmp = arr[i];
                        arr[i]=arr[i+1];
                        arr[i+1]=tmp;
                    }
                }
            }
        }
    
    for(i=0;i<size;i++)
    {
        cout<<arr[i]<<endl;
    }


}