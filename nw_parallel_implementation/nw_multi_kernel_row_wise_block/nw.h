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

//__device__ int d_mat[1000*1000] ;
//__device__ int d_traceback[1000*1000];

__device__ int maximum(int f1 ,int f2 ,int f3,int *ptr)
{
      int  max = 0 ;

        if( f1 >= f2 && f1 >= f3 )  
        {
                max = f1 ;
                *ptr = -1 ;
        }
        else if( f2 > f3 )              
        {
                max = f2 ;
                *ptr = 1 ;
        }
        else
        {
                max = f3 ;
                *ptr = 0 ;
        }

    return max;
}

__device__ volatile int g_mutex;
__device__ volatile int g_mutex_sync;
__device__ volatile int Arrayin[100];
__device__ volatile int Arrayout[100];



//GPU lock-free synchronization function
__device__ void __gpu_sync_free(int goalVal)
{
    // thread ID in a block
    int tid_in_blk = threadIdx.x * blockDim.y + threadIdx.y;
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

}


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
            g_mutex = 0 ;
    }
    
    //__threadfence();
    __syncthreads();

}

__global__ void reset_g_mutex()
{
    g_mutex = 0 ;
}

__global__ void calc_matrix(char *dseq1 , char *dseq2 , int *d_mat , int *d_traceback , int l1 , int l2 , int loop)
{

    int row = (blockIdx.y * blockDim.y) + threadIdx.y ; 
    int col = (blockIdx.x * blockDim.x) + threadIdx.x ;
    const int L1  = blockDim.x ;
    const int L2  = blockDim.y ;
    const int G1  = gridDim.x  ;
    int row_back  = row ;  

    int per_row = L1 * L1 * gridDim.x ;

    extern __shared__ int s_mat[];      // Will Work for both s_mat[First l1*l2 elements] & s_traceback [Last l1*l2 elements] 

  // ##############################################################################################################

//int loop =0 ;

//for(int loop = 0 ; loop < gridDim.x  ;loop++)
//{
    //row = row_back + 4 ;

    int offset = col + row * blockDim.x * gridDim.x ;

    int new_offset = offset + per_row * loop ;  // New Offset as Global Index in Matrix

    int blockid = blockIdx.y * gridDim.x + blockIdx.x ; 
    int off = col + row * gridDim.x - blockDim.x*blockid ;

    off += (blockDim.x - gridDim.x ) * threadIdx.y;

    int off_x = off - (threadIdx.y * blockDim.x) ;
    int off_y = (off - threadIdx.x) / blockDim.x ;


    int flag=0 , flag2 = 0;

 
    //s_mat[off] = -9999  ;

    __syncthreads();
    
// Block Syncronization 
    while(blockid!=g_mutex ){}

    if(blockid==0)
    {
            g_mutex_sync =0;
    }

// Find First Row & Col of Matrix 

    if(row == 0 && loop ==0)
    {
            s_mat[off]       = GAP * col ;
            s_mat[L1*L2 + off] = -1;
            flag2 =1;

    }
    else if(col == 0)
    {
            s_mat[off]       = GAP * (row+ blockDim.y *loop) ;
            s_mat[L1*L2 +off] = 1;
            flag2 =1;
    }
    else
    {
            s_mat[off]       = -9999 ;
            s_mat[L1*L2 +off] = -9 ;
    }

    d_mat[new_offset] = s_mat[off] ;

    __syncthreads();


    if(flag2==0)
    {
        // Calculation for First row & col for blocks which are not part of 1st global row and col .
        // #############################
            if(row == 0 && col%L1 == 0 )
            {
                        int left , top , dia ;

                        left = d_mat[new_offset -1]   + GAP ;
                        top  = d_mat[new_offset -G1*L1]  + GAP ;
                        
                        if(dseq1[row-1] == dseq2[col-1])
                            dia  = d_mat[new_offset -G1*L1-1] + MATCH ;
                        else
                            dia  = d_mat[new_offset -G1*L1-1]  MISMATCH ;

                        s_mat[off] = maximum(left,top,dia,&s_mat[L1*L2 +off]) ;
                        flag = 1;

                        d_mat[new_offset] = s_mat[off] ;
            }


            for(int i=1 ; i<L1 ;i++)
            {
               
                    if(row == 0  && off_x==i)
                    {

                            int left , top , dia ;

                            left = d_mat[new_offset -1]   + GAP ;
                            top  = d_mat[new_offset -G1*L1]  + GAP ;
                           
                            if(dseq1[row-1] == dseq2[col-1])
                                dia  = d_mat[new_offset -G1*L1-1] + MATCH ;
                            else
                                dia  = d_mat[new_offset -G1*L1-1]  MISMATCH ;

                            s_mat[off] = maximum(left,top,dia,&s_mat[L1*L2 +off]) ;


                            flag = 1;

                            d_mat[new_offset] = s_mat[off] ;

                    }
                    if(col%L1 == 0 && col!=0 && off_y ==i)
                    {

                            int left , top , dia ;

                            left = d_mat[new_offset -1]   + GAP ;
                            top  = d_mat[new_offset -G1*L1]  + GAP ;
                            
                            if(dseq1[row-1] == dseq2[col-1])
                                dia  = d_mat[new_offset -G1*L1-1] + MATCH ;
                            else
                                dia  = d_mat[new_offset -G1*L1-1]  +MISMATCH ;

                            s_mat[off] =  maximum(left,top,dia,&s_mat[L1*L2 +off])  ;
                            
                            flag = 1;
                    
                            d_mat[new_offset] = s_mat[off] ;

                    }
                
              
                    __syncthreads();
            }
            // #############################################

            // Calculation of matrix using Skewing Transformation Technique 

            if (flag==0)
            {
                row = off_x ;

                int z ;
                    for( z = 2; z <= L1 + L2 - 1; z++ )  
                    {
                        
                        int mx = ( 1>z-L1 ? 1  : z-L1);
                        int mn = ( L1<z-1 ? L1 : z-1 );
                        if(z>L1)    mx++ ;


                            if(row >=mx && row<=mn)
                            {
                                int left , top , dia ;

                                left = s_mat[row*L1 + (z-row-1)]     + GAP ;
                                top  = s_mat[(row-1)*L1 + (z-row)]   + GAP ;

                                
                                    if(dseq1[row-1 + L1 * loop] == dseq2[z-row-1+ L1 * blockIdx.x])
                                        dia  = s_mat[(row-1)*L1 + (z-row-1)] + MATCH ;
                                    else
                                        dia  = s_mat[(row-1)*L1 + (z-row-1)] + MISMATCH ;

                                    s_mat[row*L1 + (z-row)] = maximum(left,top,dia,&s_mat[L1*L2 + row*L1 + (z-row)])  ;
                            
                            }

                    __syncthreads();
                    
                    }

            }       
    }         
   
    // Write Back the value to global Matrix 
    //d_mat[offset]       =  g_mutex ;
    //d_mat[new_offset]       =  blockid + loop * gridDim.x ;
    d_mat[new_offset]       =  s_mat[off] ;
    //d_traceback[new_offset] =  0 ; 
    d_traceback[new_offset] =  s_mat[L1*L2 + off] ;
    
  // Increment the Value for g_mutex used for block sync.
    if(blockid == g_mutex && off==0)
    {   atomicAdd((int *)&g_mutex, 1);
    }
    

//    __gpu_sync_free((loop+1));

     //g_mutex=0;

//    __gpu_sync_free((loop+1)*2+1);

    //__gpu_sync(gridDim.x);

    //row = row_back ;
    
//}

// ##############################################################################################################


}


void display_mat(int *mat , int *traceback , int L1 ,int L2 );
void cuda_error_check(cudaError_t err , const char *msg );
char *strrev(char *str);

void nw(                                                          
        char*       seq_1,          /*  Needleman-Wunsch   */
        char*       seq_2,          /*  algorithm for      */
        char*       seq_1_al,       /*  global alignment   */
        char*       seq_2_al        /*  of nt sequence.    */
      )
{
        const int  L1 = strlen(seq_1) + 1 ;
        const int  L2 = strlen(seq_2) + 1 ;
           

        //int  mat[L1*L2];        // Dynamic Prog. Matrix
        int *mat = new int[L1*L2] ;
        int *traceback = new int[L1*L2] ;
        //int  traceback[L1*L2];   // TraceBack Direction matrix

// ################################# Sequential  ##################################################     
  /*     
        for(int i=0 ;i < L1 ; i++){
            mat[i] = i * GAP ;
            traceback[i] = -1 ;
        }
        for(int i=0 ;i < L2 ; i++){
            mat[i*L1] = i * GAP ;
            traceback[i*L1] = -1 ;
        }
  */
// #################################  Parallel #####################################################        
        int *d_mat, *d_traceback;
        char *dseq1 , *dseq2;

        cuda_error_check(cudaSetDevice(0) , "cudaSetDevice failed!" );


        cuda_error_check(cudaMalloc((void **)&d_mat       , L1 * L2 * sizeof(int)),"cudaMalloc Failed! 1 ");
        cuda_error_check(cudaMalloc((void **)&d_traceback , L1 * L2 * sizeof(int)),"cudaMalloc Failed! 2 ");


        cuda_error_check(cudaMalloc((void **)&dseq1 , L1 * sizeof(char)),"cudaMalloc Failed!");
        cuda_error_check(cudaMemcpy(dseq1 , seq_1 , L1 * sizeof(char) , cudaMemcpyHostToDevice),"cudaMemcpy H-D failed!");
        cuda_error_check(cudaMalloc((void **)&dseq2 , L2 * sizeof(char)),"cudaMalloc Failed!");
        cuda_error_check(cudaMemcpy(dseq2 , seq_2 , L2 * sizeof(char) , cudaMemcpyHostToDevice),"cudaMemcpy H-D failed!");

        //init<<< L2 , L1>>>();
        //for(int i=1;i<L1+L2-2;i++)
        int BLOCK_SIZE = 30;
        int GRID_SIZE  = L1/30  ;

        dim3 dimBlock(BLOCK_SIZE , BLOCK_SIZE );
        dim3 dimGrid(GRID_SIZE , 1);

        for(int loop = 0 ; loop < GRID_SIZE  ; loop++)
            {
            calc_matrix<<< dimGrid , dimBlock , 2 * BLOCK_SIZE * BLOCK_SIZE * sizeof(int) >>> ( dseq1 ,dseq2 , d_mat , d_traceback , L1 , L2 , loop) ;
            reset_g_mutex <<<1,1>>>();
            }
        //calc_matrix<<< dimGrid , dimBlock >>> ( dseq1 ,dseq2 , d_mat , d_traceback , L1 , L2) ;

        printf("\nKernel Complete Size : %d %d : %d ", L1 ,L2 , L1*L2 );

        cuda_error_check(cudaMemcpy(mat       , d_mat       , L1 * L2 * sizeof(int) , cudaMemcpyDeviceToHost),"cudaMemcpy D-H failed! 1");
        cuda_error_check(cudaMemcpy(traceback , d_traceback , L1 * L2 * sizeof(int) , cudaMemcpyDeviceToHost),"cudaMemcpy D-H failed! 2");
        printf("\nCudamemcpy D-H Complete");
       
        //cuda_error_check(cudaMemcpyFromSymbol(&mat , d_mat , sizeof(int)*L1*L2 , 0 , cudaMemcpyDeviceToHost),"cudaMemcpy D-H failed!");
        //cuda_error_check(cudaMemcpyFromSymbol(&traceback , d_traceback , sizeof(int)*L1*L2 , 0 , cudaMemcpyDeviceToHost),"cudaMemcpy D-H failed!");

        cudaFree(dseq1);
        cudaFree(dseq2);
        cudaFree(d_mat);
        cudaFree(d_traceback);

        //cudaDeviceReset();

// #################################################################################################        
/*
        printf("\nMatrix: \n");

        cout << "        ";
        for( int j = 0; j < L1; j++ )
        {
                cout << seq_2[ j ] << "   ";
        }
        cout << "\n  ";

        for( int i = 0; i < L2; i++ )
        {
                if( i > 0 )
                {
                        cout << seq_1[ i-1 ] << " ";
                }
                for( int j = 0; j < L1; j++ )
                {
                        cout.width(3);
                        cout << mat[i*L1 + j] << " ";
                }
                cout << endl;
        }
*/
/*
        printf("\nDirections: \n");
        cout << "      ";
        for( int j = 0; j < L1; j++ )
        {
                cout.width( 3 );
                cout << seq_1[ j ] << " ";
        }
        cout << "\n  ";

        for( int i = 0; i < L2; i++ )
        {
                if( i > 0 )
                {
                        cout << seq_2[ i-1 ] << " ";
                }
                for( int j = 0; j < L1; j++ )
                {
                        cout.width( 3 );
                        cout << traceback[i*L1 +j] << " ";
                }
                cout << endl;
        }*/

// #################################################################################################
        
        int i=L1-1 , j = L2-1;


        char  *temp;
        temp = (char *)malloc((L1+100) * sizeof(char));


        while( i > 0 || j > 0 )
        {
                switch( traceback[ i *L1 + j ] )
                {
                        case 1 :        strcat(seq_2_al,"-");
                                        sprintf(temp,"%c",seq_1[i-1]);
                                        strcat(seq_1_al,temp);
                                        i-- ;
                                        break ;

                        case 0 :        sprintf(temp,"%c",seq_2[j-1]);
                                        strcat(seq_2_al,temp);
                                        sprintf(temp,"%c",seq_1[i-1]);                                            
                                        strcat(seq_1_al,temp);
                                        i-- ;  j-- ;
                                        break ;

                        case -1 :       sprintf(temp,"%c",seq_2[j-1]);
                                        strcat(seq_2_al,temp);
                                        strcat(seq_1_al,"-");
                                        j-- ;
                }
                //k++ ;
        }

        strcpy(temp,strrev(seq_1_al));
        strcpy(seq_1_al,temp);
        strcpy(temp,strrev(seq_2_al));
        strcpy(seq_2_al,temp);

        //return 0 ;

        //traceback[0]=0;
        //display_mat(mat,L1,L2);
  
}

void display_mat(int *mat , int L1 ,int L2 )
{
    printf("\nMatrix: \n");
    for(int i=0 ;i < L2 ; i++){
        for(int j=0 ;j<L1 ;j++)
            printf("%d\t",mat[i*L1 + j]);
        printf("\n");}

/*  printf("\nDirections: \n");
    for(int i=0 ;i < L2 ; i++){
        for(int j=0 ;j<L1 ;j++)
            printf("%d\t",traceback[i*L1 +j]);
        printf("\n");}*/
}

void cuda_error_check(cudaError_t err , const char *msg )
{
    if(err != cudaSuccess)
    {
      printf("The error is %s, %s \n", cudaGetErrorString(err), msg );
      exit(1);
    }
}

char *strrev(char *str)
{
      char *p1, *p2;

      if (! str || ! *str)
            return str;
      for (p1 = str, p2 = str + strlen(str) - 1; p2 > p1; ++p1, --p2)
      {
            *p1 ^= *p2;
            *p2 ^= *p1;
            *p1 ^= *p2;
      }
      return str;
}