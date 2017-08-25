#include <cuda.h>
#include <cuda_runtime_api.h>
#include <stdio.h>
#include <iostream>
#include <string.h>
#include <algorithm>
#include <stdlib.h>

#define N 4
#define BLOCK_SIZE 4
#define GRID_SIZE 1

void cuda_error_check(cudaError_t err , const char *msg )
{
    if(err != cudaSuccess)
    {
      printf("The error is %s, %s \n", cudaGetErrorString(err), msg );
      exit(1);
    }
}

__global__ void matrix(int *d_a)
{

int row = (blockIdx.y * blockDim.y) + threadIdx.y ; 
int col = (blockIdx.x * blockDim.x) + threadIdx.x ;
int L1  = blockDim.x ;
int L2  = blockDim.y ;
       
__shared__ int s_a[N*N];


s_a[row*L1 + col] =  d_a[row*L1 + col] ;

__syncthreads();


d_a[row*L1 + col] = s_a[row*L1 + col] ;

/*if(row == 0)
{
    d_a[col] = col * -1 ;
}else if(col == 0)
{
    d_a[row*L1] = row * -1 ;
}else
{
    d_a[row*L1 + col] = 0 ;
}

__syncthreads();

        int z ;
        for( z = 2; z <= L1 + L2 - 1; z++ )  
        {
            
            int mx = ( 1>z-L1 ? 1  : z-L1);
            int mn = ( L1<z-1 ? L1 : z-1 );
            if(z>L1)    mx++ ;


                if(row >=mx && row<=mn)
                {
                    d_a[row*L1 + (z-row)] = d_a[(row-1)*L1 + (z-row)] + d_a[row*L1 + (z-row-1)]   ;
                }

        __syncthreads();
        
        }
*/

}

int main(int argc , char **argv)
{

int a[N*N];
int i,j;


 


int *d_a ;


    cuda_error_check(cudaSetDevice(0) , "cudaSetDevice failed!" );

    cuda_error_check(cudaMalloc((void **)&d_a , N*N* sizeof(int)),"cudaMalloc Failed!");
    cuda_error_check(cudaMemcpy(d_a , a , N*N * sizeof(int) , cudaMemcpyHostToDevice),"cudaMemcpy H-D failed!");

    dim3 dimBlock(BLOCK_SIZE , BLOCK_SIZE );
    dim3 dimGrid(GRID_SIZE , GRID_SIZE);

    matrix<<< dimGrid , dimBlock >>>(d_a);

    cuda_error_check(cudaMemcpy(a , d_a , N*N * sizeof(int) , cudaMemcpyDeviceToHost),"cudaMemcpy D-H failed!");
 
printf("\n\n");

 for(i=0;i<N;i++)
 {
    for(j=0;j<N;j++)
    {
        //printf("%d\t",a[i*N + j]);
        printf("%d\t",a[i*N + j]);
    }
    printf("\n");
 }
return 0;
}
