#include <cuda.h>
#include <cuda_runtime_api.h>
#include <stdio.h>
#include <iostream>
#include <string.h>
#include <algorithm>
#include <stdlib.h>

//#define N 4
#define BLOCK_SIZE 4
#define GRID_SIZE 2


using namespace std;

__device__ volatile int g_mutex;

void cuda_error_check(cudaError_t err , const char *msg )
{
    if(err != cudaSuccess)
    {
      printf("The error is %s, %s \n", cudaGetErrorString(err), msg );
      exit(1);
    }
}

__device__ void __gpu_sync(int goalVal)
{
    //thread ID in a block
    int tid_in_block = threadIdx.x * blockDim.y + threadIdx.y;
    // only thread 0 is used for synchronization
    if (tid_in_block == 0)
    {
            atomicAdd((int *)&g_mutex, 1);

            while(g_mutex != goalVal)
            {
            //Do nothing here
            }        
    }
    
    //__threadfence();
    __syncthreads();
}

__global__ void matrix(int *d_a)
{

    int row = (blockIdx.y * blockDim.y) + threadIdx.y ; 
    int col = (blockIdx.x * blockDim.x) + threadIdx.x ;
    //int L1  = blockDim.x ;
    //int L2  = blockDim.y ;
      
    int offset = col + row * blockDim.x * gridDim.x ;

    int blockid = blockIdx.y * gridDim.x + blockIdx.x ;
    int off = col + row * gridDim.x  - blockDim.x*blockid ;

    off += (blockDim.x - gridDim.x ) * threadIdx.y;

    int off_x = off - (threadIdx.y * blockDim.x) ;
    int off_y = (off - threadIdx.x) / blockDim.x ;

    //while(blockid!=g_mutex ){}

    if(blockid==0)
            d_a[offset] = g_mutex ;  


    __gpu_sync(4);

    if(blockid!=0)
            d_a[offset] = g_mutex ;  
/*
    if(blockid == g_mutex && off==0)
    {   atomicAdd((int *)&g_mutex, 1);
    }*/

}

int main(int argc , char **argv)
{

    int a[(GRID_SIZE*GRID_SIZE)*BLOCK_SIZE*BLOCK_SIZE];
    int i,j;
    //int m,n;
     
    for(j=0;j<(GRID_SIZE*GRID_SIZE)*BLOCK_SIZE*BLOCK_SIZE;j++)
    {
        a[j]=0;
    }

    int *d_a ;

    cuda_error_check(cudaSetDevice(0) , "cudaSetDevice failed!" );

    cuda_error_check(cudaMalloc((void **)&d_a , (GRID_SIZE*GRID_SIZE)*BLOCK_SIZE*BLOCK_SIZE* sizeof(int)),"cudaMalloc Failed!");
    cuda_error_check(cudaMemcpy(d_a , a , (GRID_SIZE*GRID_SIZE)*BLOCK_SIZE*BLOCK_SIZE * sizeof(int) , cudaMemcpyHostToDevice),"cudaMemcpy H-D failed!");

    dim3 dimBlock(BLOCK_SIZE , BLOCK_SIZE );
    dim3 dimGrid(GRID_SIZE , GRID_SIZE);

    matrix<<< dimGrid , dimBlock >>>(d_a);

    cuda_error_check(cudaMemcpy(a , d_a , (GRID_SIZE*GRID_SIZE)*BLOCK_SIZE*BLOCK_SIZE * sizeof(int) , cudaMemcpyDeviceToHost),"cudaMemcpy D-H failed!");


for(j=0;j<(GRID_SIZE*GRID_SIZE)*BLOCK_SIZE*BLOCK_SIZE;j++)
{
    if(j%(BLOCK_SIZE*GRID_SIZE) == 0)
        cout<<endl;
    if(j%(GRID_SIZE*BLOCK_SIZE*BLOCK_SIZE)==0)
        cout<<endl;
    if(j%BLOCK_SIZE ==0)
        cout<<"  ";
        
    cout.width( 3 );
    cout<<a[j]<<" ";
}
            printf("\n\n");

 /*for(m=0;m<2;m++)
 {
    for(n=0;n<2;n++)
    {*/
   /*
         for(i=0;i<2*BLOCK_SIZE;i++)
         {
            for(j=0;j<2*BLOCK_SIZE;j++)
            {                  // [m*2+n]
                //printf("%d\t",a[i*BLOCK_SIZE + j]);
                cout.width( 3 );
                //cout<<a[(m*2+n)*16 + i*BLOCK_SIZE + j]<<" ";
                cout<<a[i*(2*BLOCK_SIZE) + j]<<" ";

            }
            printf("\n");
         }

         */
        /*    printf("\n");
    }
}
*/
return 0;
}





 /*      
if(blockIdx.x == 0 && blockIdx.y == 1)
    d_a[offset] = offset;
else
    d_a[offset] = -1;

if(row==0)
        d_a[offset] = -1 * offset;
if(col==0)
        d_a[offset] = -1 * offset/(BLOCK_SIZE*BLOCK_SIZE);
*/