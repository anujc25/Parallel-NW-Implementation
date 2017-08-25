#include <cuda.h>
#include <cuda_runtime_api.h>
#include <stdio.h>
#include <iostream>
#include <string.h>
#include <algorithm>
#include <stdlib.h>

#define N         40
#define GAP       -2
#define MATCH      1
#define MISMATCH  -1

//#include "kernels.h"

using namespace std ;


__device__ volatile int g_mutex;
__device__ volatile int g_mutex_sync;

__device__ void __gpu_sync(int goalVal)
{
    //thread ID in a block
    int tid_in_block = threadIdx.x * blockDim.y + threadIdx.y;
    // only thread 0 is used for synchronization
    if (tid_in_block == 0)
    {
            atomicAdd((int *)&g_mutex_sync, 1);

            while(g_mutex_sync != goalVal)
            {
            //Do nothing here
            }        
            //g_mutex = 0 ;
    }
    
    //__threadfence();
    __syncthreads();

}

__global__ void calc_matrix(int *d_mat , int l1 , int l2)
{  
    int row = (blockIdx.y * blockDim.y) + threadIdx.y ; 

    int col = (blockIdx.x * blockDim.x) + threadIdx.x ;

    int offset = col + row * blockDim.x * gridDim.x ;

    __gpu_sync(1);

    d_mat[offset]       =  gridDim.x     ;

}



void cuda_error_check(cudaError_t err , const char *msg );

int main()
{
        const int NN = 64 ;
        const int  L1 = NN ;
        const int  L2 = NN ;
           

      //  int  mat[L1*L2];        // Dynamic Prog. Matrix

        int *mat = new int[L1*L2] ;
        

        //printf("\nSize : %d\n",sizeof(mat));
// #################################  Parallel #####################################################        
        int *d_mat ;
        
        cuda_error_check(cudaSetDevice(0) , "cudaSetDevice failed!" );


        cuda_error_check(cudaMalloc((void **)&d_mat       , L1 * L2 * sizeof(int)),"cudaMalloc Failed!");

        int BLOCK_SIZE = 16;
        int GRID_SIZE  = 1  ;

        dim3 dimBlock(BLOCK_SIZE , BLOCK_SIZE );
        dim3 dimGrid(GRID_SIZE , 1);

        calc_matrix<<< dimBlock , dimGrid >>> ( d_mat , L1 , L2 ) ;

        printf("\nKernel Complete Size : %d %d : %d ", L1 ,L2 , L1*L2 );

        cuda_error_check(cudaMemcpy(mat , d_mat ,L1* L2 * sizeof(int) , cudaMemcpyDeviceToHost),"cudaMemcpy D-H failed! 1");
        printf("\nCudamemcpy D-H Complete");

        cudaFree(d_mat);

       // cudaDeviceReset();

// #################################################################################################        

        printf("\nMatrix: \n");

        cout << "        ";
        for( int j = 0; j < L1; j++ )
        {
                //cout << seq_2[ j ] << "   ";
        }
        cout << "\n";

        for( int i = 0; i < L2; i++ )
        {
                if( i > 0 )
                {
                        //cout << seq_1[ i-1 ] << " ";
                }
                for( int j = 0; j < L1; j++ )
                {
                        cout.width(3);
                        cout << mat[i*L1 + j] << " ";
                }
                cout << endl;
        }

// #################################################################################################
        
}

void cuda_error_check(cudaError_t err , const char *msg )
{
    if(err != cudaSuccess)
    {
      printf("The error is %s, %s \n", cudaGetErrorString(err), msg );
      exit(1);
    }
}
