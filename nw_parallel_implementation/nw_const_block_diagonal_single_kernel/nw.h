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




__global__ void calc_matrix( char *dseq1 ,char *dseq2 ,int *d_mat ,int *d_traceback , int L1 ,int L2 ,int BLK_P_ROW)
{

    int thread_id =  threadIdx.x + blockIdx.x * blockDim.x   ;
    int size      =  blockDim.x * gridDim.x ;

    int grd_size  = gridDim.x  ;
    int blk_size  = blockDim.x ; 
    int blk_id    = blockIdx.x ;
/*    __shared__ int count ;
    if(blk_id==0 && threadIdx.x==0)
        count=0;*/

    for(int ii = 0 ; ii < grd_size ; ii++)  // Block Sync 1st N diagonals Blocks
    {
        thread_id = threadIdx.x ;
        int blk_r = blk_id ;
        int blk_c = ii - blk_id ;

        if(blk_id<=ii)
        {   
             int offset = blk_c*blk_size + blk_r*size*blk_size + threadIdx.x  ;

                
                // First Row and First Column Calculation
                if(offset<size)
                    d_mat[offset] = offset * GAP ;
                if (offset%size == 0)
                {   
                    for(int xx =0 ; xx < blk_size ; xx++)
                        d_mat[offset+size*xx] = ((offset+size*xx)/size) * GAP ;                
                }

          
                // Inside Block Sync 1st N diagonals
                    for(int jj=0 ; jj< blk_size ;jj++)  
                    {
                        //int offset = blk_c*blk_size + blk_r*size*blk_size + size*jj + threadIdx.x  ;

                        int r = thread_id + blk_size * blk_r ;
                        int c = jj - thread_id + blk_size * blk_c ;

                        if(thread_id<=jj && r!=0 && c!=0 )
                        {
                            //d_mat[ r * size + c] =   blk_id ; // r * size + c ;

                            int left , top , dia ;

                            left = d_mat[r * size + c -1]   + GAP ;
                            top  = d_mat[(r-1) * size + c]  + GAP ;
                            
                            if(dseq1[r-1] == dseq2[c-1])
                                dia  = d_mat[(r-1) * size + c -1] + MATCH ;
                            else
                                dia  = d_mat[(r-1) * size + c -1] + MISMATCH ;

                            d_mat[ r * size + c] =   maximum(left,top,dia , &d_traceback[r*size+c] ) ;
                        }
                        __syncthreads();
                    } // Inside Block Sync 1st N diagonals

                    thread_id = thread_id + 1 ;

                // Inside Block Sync Last N-1 diagonals
                    for(int j = 1 ; j < blk_size ; j++)
                    {
                        int r = blk_size-1-thread_id+j +  blk_size * blk_r ;
                        int c = thread_id +  blk_size * blk_c  ;
                
                        if(thread_id < blk_size )
                        {
                            //d_mat[ r * size + c] =  blk_id ; //r * size + c ;

                            int left , top , dia ;

                            left = d_mat[r * size + c -1]   + GAP ;
                            top  = d_mat[(r-1) * size + c]  + GAP ;
                            
                            if(dseq1[r-1] == dseq2[c-1])
                                dia  = d_mat[(r-1) * size + c -1] + MATCH ;
                            else
                                dia  = d_mat[(r-1) * size + c -1] + MISMATCH ;

                            d_mat[ r * size + c] =   maximum(left,top,dia , &d_traceback[r*size+c] ) ;

                        }
                        thread_id = thread_id + 1 ;
                         __syncthreads();

                    }   // Inside Block Sync Last N-1 diagonals
                
        } // Block selection Condition 

        __gpu_sync_free(ii+1);
    
    } // Block Sync 1st N diagonals Blocks

    blk_id = blk_id + 1 ;

    for(int pp = 1 ; pp < grd_size ; pp++)
    {
        thread_id = threadIdx.x ;
        int blk_r = grd_size -1- blk_id + pp ;
        int blk_c = blk_id   ;

        if( blk_id < grd_size ) // Block selection Condition 
        {
                //int offset = blk_c*blk_size + blk_r*size*blk_size + threadIdx.x  ;

                //d_mat[offset] = jj ;

                // Inside Block Sync 1st N diagonals
                    for(int jj=0 ; jj< blk_size ;jj++)  
                    {
                        //int offset = blk_c*blk_size + blk_r*size*blk_size + size*jj + threadIdx.x  ;

                        int r = thread_id + blk_size * blk_r ;
                        int c = jj - thread_id + blk_size * blk_c ;

                        if(thread_id<=jj && r!=0 && c!=0 )
                        {
                            //d_mat[ r * size + c] =   blk_id ; // r * size + c ;

                            int left , top , dia ;

                            left = d_mat[r * size + c -1]   + GAP ;
                            top  = d_mat[(r-1) * size + c]  + GAP ;
                            
                            if(dseq1[r-1] == dseq2[c-1])
                                dia  = d_mat[(r-1) * size + c -1] + MATCH ;
                            else
                                dia  = d_mat[(r-1) * size + c -1] + MISMATCH ;

                            d_mat[ r * size + c] =   maximum(left,top,dia , &d_traceback[r*size+c] ) ;
                        }
                        __syncthreads();
                    } // Inside Block Sync 1st N diagonals

                    thread_id = thread_id + 1 ;

                // Inside Block Sync Last N-1 diagonals
                    for(int j = 1 ; j < blk_size ; j++)
                    {
                        int r = blk_size-1-thread_id+j +  blk_size * blk_r ;
                        int c = thread_id +  blk_size * blk_c  ;
                
                        if(thread_id < blk_size )
                        {
                            //d_mat[ r * size + c] =  blk_id ; //r * size + c ;

                            int left , top , dia ;

                            left = d_mat[r * size + c -1]   + GAP ;
                            top  = d_mat[(r-1) * size + c]  + GAP ;
                            
                            if(dseq1[r-1] == dseq2[c-1])
                                dia  = d_mat[(r-1) * size + c -1] + MATCH ;
                            else
                                dia  = d_mat[(r-1) * size + c -1] + MISMATCH ;

                            d_mat[ r * size + c] =   maximum(left,top,dia , &d_traceback[r*size+c] ) ;

                        }
                        thread_id = thread_id + 1 ;
                         __syncthreads();

                    }   // Inside Block Sync Last N-1 diagonals

        } // Block selection Condition 

        blk_id = blk_id + 1 ;

        __gpu_sync_free(pp+100);

    }


}


//###################################### Threads in Block Logic ##############################################

                /*  d_mat[thread_id] = thread_id * GAP ; 
                    d_traceback[thread_id] = -1 ;

                    d_mat[thread_id*(size)] = thread_id * GAP ; 
                    d_traceback[thread_id*(size)] = 1 ; 
                


                    for(int i =0 ; i < size ; i++)
                    {
                        int r = thread_id ;
                        int c = i - thread_id  ;

                        if(thread_id<=i && r!=0 && c!=0 )
                        {
                            d_mat[ r * size + c] =   blk_id ; // r * size + c ;

                            int left , top , dia ;

                            left = d_mat[r * size + c -1]   + GAP ;
                            top  = d_mat[(r-1) * size + c]  + GAP ;
                            
                            if(dseq1[r-1] == dseq2[c-1])
                                dia  = d_mat[(r-1) * size + c -1] + MATCH ;
                            else
                                dia  = d_mat[(r-1) * size + c -1] + MISMATCH ;

                            d_mat[ r * size + c] =   maximum(left,top,dia , &d_traceback[r*size+c] ) ;
                        }
                        __syncthreads();

                    }

                    //thread_id = thread_id + 1 ;

                    for(int j = 1 ; j < size ; j++)
                    {
                        int r = size-1-thread_id+j ;
                        int c = thread_id   ;


                
                        if(thread_id < size )
                        {
                            d_mat[ r * size + c] =  blk_id ; //r * size + c ;

                            /*int left , top , dia ;

                            left = d_mat[r * size + c -1]   + GAP ;
                            top  = d_mat[(r-1) * size + c]  + GAP ;
                            
                            if(dseq1[r-1] == dseq2[c-1])
                                dia  = d_mat[(r-1) * size + c -1] + MATCH ;
                            else
                                dia  = d_mat[(r-1) * size + c -1] + MISMATCH ;

                            d_mat[ r * size + c] =   maximum(left,top,dia , &d_traceback[r*size+c] ) ;

                        }
                        thread_id = thread_id + 1 ;
                         __syncthreads();

                    }*/

        
    



void display_mat(int *mat , int *traceback , int L1 ,int L2 );
void cuda_error_check(cudaError_t err , const char *msg );
char *strrev(char *str);

void nw(                                                          
        char*       seq_1,          /*  Needleman-Wunsch   */
        char*       seq_2,          /*  algorithm for      */
        char*       seq_1_al,       /*  global alignment   */
        char*       seq_2_al,        /*  of nt sequence.    */
        int         given_block_size
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


        cuda_error_check(cudaMalloc((void **)&d_mat       , L1 * L2 * sizeof(int)),"cudaMalloc Failed!");
        cuda_error_check(cudaMalloc((void **)&d_traceback , L1 * L2 * sizeof(int)),"cudaMalloc Failed!");


        cuda_error_check(cudaMalloc((void **)&dseq1 , L1 * sizeof(char)),"cudaMalloc Failed!");
        cuda_error_check(cudaMemcpy(dseq1 , seq_1 , L1 * sizeof(char) , cudaMemcpyHostToDevice),"cudaMemcpy H-D failed!");
        cuda_error_check(cudaMalloc((void **)&dseq2 , L2 * sizeof(char)),"cudaMalloc Failed!");
        cuda_error_check(cudaMemcpy(dseq2 , seq_2 , L2 * sizeof(char) , cudaMemcpyHostToDevice),"cudaMemcpy H-D failed!");

        //init<<< L2 , L1>>>();
        //for(int i=1;i<L1+L2-2;i++)


        int BLOCK_SIZE = given_block_size ; 
        int BLK_P_ROW  = L1/BLOCK_SIZE  ;
        int GRID_SIZE  = L1/BLOCK_SIZE  ;

        dim3 dimBlock(BLOCK_SIZE , 1);
        dim3 dimGrid(GRID_SIZE , 1);

        calc_matrix<<< GRID_SIZE , BLOCK_SIZE  >>> (dseq1 ,dseq2 , d_mat , d_traceback , L1 , L2 , BLK_P_ROW) ;
        //calc_matrix<<< dimGrid , dimBlock , 2 * BLOCK_SIZE * BLOCK_SIZE * sizeof(int) >>> ( dseq1 ,dseq2 , d_mat , d_traceback , L1 , L2 , BLK_P_ROW) ;
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
        }*/

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
        }
*/
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





















/*
    int row = (blockIdx.y * blockDim.y) + threadIdx.y ; 
    int col = (blockIdx.x * blockDim.x) + threadIdx.x ;
    const int L1  = blockDim.x ;
    const int L2  = blockDim.y ;
    const int G1  = BLK_P_ROW  ;
    //int row_back  = row ;  
    int col_back  = col ;
    int per_row = L1 * L1 * G1 ;

    extern __shared__ int s_mat[];      // Will Work for both s_mat[First l1*l2 elements] & s_traceback [Last l1*l2 elements] 

  // ##############################################################################################################
    int loop=0,loopin=0;

                        int offset = col + row * blockDim.x * G1 ;
                        int new_offset = offset + per_row * loop ;  // New Offset as Global Index in Matrix



                        int blockid = blockIdx.y * G1 + blockIdx.x + loopin ; 
                        int off = col + row * G1 - blockDim.x*blockid ;

                        off += (blockDim.x - G1 ) * threadIdx.y;

                        int off_x = off - (threadIdx.y * blockDim.x) ;
                        int off_y = (off - threadIdx.x) / blockDim.x ;


                        int flag=0 , flag2 = 0;

                     
                        //s_mat[off] = -9999  ;

                        __syncthreads();
                        
                    //////////////////////////////////////////////////////////////////////////////////////////

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

                                
                                if (flag==0)   // Calculation of matrix using Skewing Transformation Technique 
                                {
                                   int row = off_y ;

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

                                }       // Calculation of matrix using Skewing Transformation Technique  
                        }   // Calculation for First row & col for blocks which are not part of 1st global row and col .

                       
                        // Write Back the value to global Matrix 
                        //d_mat[offset]         =  g_mutex ;
                        //d_mat[new_offset]     =  off_x ;
                        d_mat[new_offset]       =  s_mat[off] ;
                        //d_traceback[new_offset] =  0 ; 
                        d_traceback[new_offset] =  s_mat[L1*L2 + off] ;
                        
                  */