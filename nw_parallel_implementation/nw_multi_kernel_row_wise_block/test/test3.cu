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


// threadfence();

using namespace std;

__device__ volatile int Arrayin[100];
__device__ volatile int Arrayout[100];


void cuda_error_check(cudaError_t err , const char *msg )
{
    if(err != cudaSuccess)
    {
      printf("The error is %s, %s \n", cudaGetErrorString(err), msg );
      exit(1);
    }
}

__device__ void __gpu_sync(int goalVal , volatile int *Arrayin , volatile int *Arrayout)
{
    // thread ID in a block
    int tid_in_blk = threadIdx.x * blockDim.y  + threadIdx.y;
    int nBlockNum = gridDim.x * gridDim.y;
    int bid = blockIdx.x * gridDim.y + blockIdx.y;


    // only thread 0 is used for synchronization
    if (tid_in_blk == 0)
    {
        Arrayin[bid] = goalVal;
    }
    
    if (bid == 1)
    {
        if (tid_in_blk < nBlockNum)
        {
            while (Arrayin[tid_in_blk] != goalVal){
            //Do nothing here
            }
        }

        __syncthreads();
        __threadfence();

        if (tid_in_blk < nBlockNum)
        {
            Arrayout[tid_in_blk] = goalVal;
        }
    }
    if (tid_in_blk == 0)
    {
            while (Arrayout[bid] != goalVal){
            //Do nothing here
            }
    }

    __syncthreads();
    __threadfence();

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
    __gpu_sync(4, Arrayin , Arrayout) ;

    d_a[offset] = 0 ;  

    /*if(blockid == g_mutex && off==0)
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

    cudaDeviceReset();

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

