#include <cuda.h>
#include <cuda_runtime_api.h>
#include <stdio.h>
#include <iostream>
#include <string.h>
#include <algorithm>
#include <stdlib.h>

#define N         10
#define GAP       -2
#define MATCH      1
#define MISMATCH  -1

//#include "kernels.h"

using namespace std ;

//__device__ int d_mat[2000*2000] ;
//__device__ int d_traceback[2000*2000];


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
/*
__global__ void init()
{
    int id = ( blockIdx.x * blockDim.x ) + threadIdx.x ;
    
    if(blockIdx.x == 0)
    {
        d_mat[id] = GAP * threadIdx.x ;
        d_traceback[id] = -1;
    }
    else if(threadIdx.x == 0)
    {
        d_mat[id] = GAP * blockIdx.x ;
        d_traceback[id] = 1;
    }  
    else
    {
        d_mat[id] =  -999 ;
        d_traceback[id] = 0;
    }
}
*/

__global__ void calc_matrix(int Num ,char *dseq1 , char *dseq2 , int *d_mat , int *d_traceback ,int L1 , int L2)
{
 
    int row = (blockIdx.y * blockDim.y) + threadIdx.y ; 
    int col = (blockIdx.x * blockDim.x) + threadIdx.x ;

    int blockid = blockIdx.y * gridDim.x + blockIdx.x ; 

    //int tx = threadIdx.x ;
    //int ty = threadIdx.y ;
    //int bx = blockIdx.x  ;
    int by = blockIdx.y  ;
    //int bd = blockDim.x  ;
    int gdx = gridDim.x   ;

    int offset = col + row * blockDim.x * gridDim.x ;   // Global Index

    int off = col + row * gridDim.x - blockDim.x*blockid ;  // index in current Block

    off += (blockDim.x - gridDim.x ) * threadIdx.y;

    int off_x = off - (threadIdx.y * blockDim.x) ;
    int off_y = (off - threadIdx.x) / blockDim.x ;


    extern __shared__ int s_mat[];      // Will Work for both s_mat[First l1*l2 elements] & s_traceback [Last l1*l2 elements] 


    //d_mat[offset] = off ;

        int mx = ( 1 > Num-gdx ? 1  : Num-gdx);
        int mn = ( gdx< Num-1  ? gdx : Num-1 );

// ############################################################################################ //

        if(by >=(mx-1) && by<=(mn-1))       // Condition for Block Range 
        {         
            if(blockid == (by*gdx + Num-by -2))     // Condition to enter the block End                             //d_mat[by*gdx*2*2 + (Num-by-2)*2  ] = Num ;
            {
                int loop = Num - 2 ;  
            // ################################################################################ //     
                      d_mat[offset] = Num ;
                      /*
                    int flag=0 , flag2 = 0;
                    //__syncthreads();
                    
                    // Block Syncronization 
                        //while(blockid!=g_mutex ){}

                    // Find First Row & Col of Matrix 
                    
                        if(row == 0)
                        {
                                s_mat[off]       = GAP * col ;
                                s_mat[L1*L2 + off] = -1;
                                flag2 =1;

                        }
                        else if(col == 0)
                        {
                                s_mat[off]       = GAP * row ;
                                s_mat[L1*L2 +off] = 1;
                                flag2 =1;
                        }
                        else
                        {
                                s_mat[off]       = -9999 ;
                                s_mat[L1*L2 +off] = -9 ;
                        }

                        d_mat[offset] = s_mat[off] ;
                        //d_mat[offset] = offset;
                        //__syncthreads();

                        if(flag2==0)
                        {
                            // Calculation for First row & col for blocks which are not part of 1st global row and col .
                            // #############################
                                if(row == 0 && col%L1 == 0 )
                                {
                                            int left , top , dia ;

                                            left = d_mat[offset -1]   + GAP ;
                                            top  = d_mat[offset -gdx*L1]  + GAP ;
                                            
                                            if(dseq1[row-1] == dseq2[col-1])
                                                dia  = d_mat[offset -gdx*L1-1] + MATCH ;
                                            else
                                                dia  = d_mat[offset -gdx*L1-1] + MISMATCH ;

                                            s_mat[off] = maximum(left,top,dia,&s_mat[L1*L2 +off]) ;
                                            flag = 1;

                                            d_mat[offset] = s_mat[off] ;
                                }


                                for(int i=1 ; i<L1 ;i++)
                                {
                                   
                                     if(row != 0  && row%L1==0 && col%L1==0)
                                        {

                                                int left , top , dia ;

                                                left = d_mat[offset -1]   + GAP ;
                                                top  = d_mat[offset -gdx*L1]  + GAP ;
                                               
                                                if(dseq1[row-1] == dseq2[col-1])
                                                    dia  = d_mat[offset -gdx*L1-1] + MATCH ;
                                                else
                                                    dia  = d_mat[offset -gdx*L1-1] + MISMATCH ;

                                                s_mat[off] = maximum(left,top,dia,&s_mat[L1*L2 +off]) ;


                                                flag = 1;

                                                d_mat[offset] = s_mat[off] ;

                                        }

                                        if(row == 0  && off_x==i)
                                        {

                                                int left , top , dia ;

                                                left = d_mat[offset -1]   + GAP ;
                                                top  = d_mat[offset -gdx*L1]  + GAP ;
                                               
                                                if(dseq1[row-1] == dseq2[col-1])
                                                    dia  = d_mat[offset -gdx*L1-1] + MATCH ;
                                                else
                                                    dia  = d_mat[offset -gdx*L1-1] + MISMATCH ;

                                                s_mat[off] = maximum(left,top,dia,&s_mat[L1*L2 +off]) ;


                                                flag = 1;

                                                d_mat[offset] = s_mat[off] ;

                                        }

                                        if(row%L1 == 0 && row!=0 && off_x ==i)
                                        {

                                                int left , top , dia ;

                                                left = d_mat[offset -1]   + GAP ;
                                                top  = d_mat[offset -gdx*L1]  + GAP ;
                                                
                                                if(dseq1[row-1] == dseq2[col-1])
                                                    dia  = d_mat[offset -gdx*L1-1] + MATCH ;
                                                else
                                                    dia  = d_mat[offset -gdx*L1-1]  +MISMATCH ;

                                                s_mat[off] =  maximum(left,top,dia,&s_mat[L1*L2 +off])  ;
                                                
                                                flag = 1;
                                        
                                                d_mat[offset] = s_mat[off] ;

                                        }

                                        if(col%L1 == 0 && col!=0 && off_y ==i)
                                        {

                                                int left , top , dia ;

                                                left = d_mat[offset -1]   + GAP ;
                                                top  = d_mat[offset -gdx*L1]  + GAP ;
                                                
                                                if(dseq1[row-1] == dseq2[col-1])
                                                    dia  = d_mat[offset -gdx*L1-1] + MATCH ;
                                                else
                                                    dia  = d_mat[offset -gdx*L1-1]  +MISMATCH ;

                                                s_mat[off] =  maximum(left,top,dia,&s_mat[L1*L2 +off])  ;
                                                
                                                flag = 1;
                                        
                                                d_mat[offset] = s_mat[off] ;

                                        }
                                    
                                  
                                       // __syncthreads();
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
                        
                        
                        d_mat[offset]       =  s_mat[off] ;
                        
                        d_traceback[offset] =  s_mat[L1*L2 + off] ;*/

                        d_traceback[offset] = 0 ;  

            }       // Condition to enter the block End
        }   // Condition for Block Range End



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
	
		int *mat = new int[L1*L2] ;
        int *traceback = new int[L1*L2] ;

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

        cuda_error_check(cudaMalloc((void **)&dseq1 , L1 * sizeof(char)),"cudaMalloc Failed! 1");
        cuda_error_check(cudaMemcpy(dseq1 , seq_1 , L1 * sizeof(char) , cudaMemcpyHostToDevice),"cudaMemcpy H-D failed!");
        cuda_error_check(cudaMalloc((void **)&dseq2 , L2 * sizeof(char)),"cudaMalloc Failed! 2");
        cuda_error_check(cudaMemcpy(dseq2 , seq_2 , L2 * sizeof(char) , cudaMemcpyHostToDevice),"cudaMemcpy H-D failed!");

    	//init<<< L2 , L1>>>();

        int BLOCK_SIZE = 30;
        int GRID_SIZE  = L1/30  ;

        dim3 dimBlock(BLOCK_SIZE , BLOCK_SIZE );
        dim3 dimGrid(GRID_SIZE , GRID_SIZE);
        
       //int i=5;
        for(int i=2;i<=L1+L1;i++)
        {
            calc_matrix<<< dimGrid , dimBlock , 2 * BLOCK_SIZE * BLOCK_SIZE * sizeof(int) >>> (i , dseq1 ,dseq2,  d_mat , d_traceback ,BLOCK_SIZE , BLOCK_SIZE) ;
        }
        
        cuda_error_check(cudaMemcpy(mat       , d_mat       , L1 * L2 * sizeof(int) , cudaMemcpyDeviceToHost),"cudaMemcpy D-H failed! 1");
        cuda_error_check(cudaMemcpy(traceback , d_traceback , L1 * L2 * sizeof(int) , cudaMemcpyDeviceToHost),"cudaMemcpy D-H failed! 2");
        printf("\nCudamemcpy D-H Complete");
       

        //cuda_error_check(cudaMemcpyFromSymbol(&mat , d_mat , sizeof(int)*L1*L2 , 0 , cudaMemcpyDeviceToHost),"cudaMemcpy D-H failed!");
        //cuda_error_check(cudaMemcpyFromSymbol(&traceback , d_traceback , sizeof(int)*L1*L2 , 0 , cudaMemcpyDeviceToHost),"cudaMemcpy D-H failed!");

        cudaFree(dseq1);
        cudaFree(dseq2);
        cudaFree(d_mat);
        cudaFree(d_traceback);

        cudaDeviceReset();
    
// #################################################################################################        
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
        //char temp[4000];

        char  *temp;
        temp = (char *)malloc((L1+1000) * sizeof(char));


        while( i > 0 || j > 0 )
        {
                switch( traceback[ i *L1 + j ] )
                {
                        case 1 :        strcat(seq_1_al,"-");
                                        sprintf(temp,"%c",seq_2[i-1]);
                                        strcat(seq_2_al,temp);
                                        i-- ;
                                        break ;

                        case 0 :        sprintf(temp,"%c",seq_1[j-1]);
                                        strcat(seq_1_al,temp);
                                        sprintf(temp,"%c",seq_2[i-1]);                                            
                                        strcat(seq_2_al,temp);
                                        i-- ;  j-- ;
                                        break ;

                        case -1 :       sprintf(temp,"%c",seq_1[j-1]);
                                        strcat(seq_1_al,temp);
                                        strcat(seq_2_al,"-");
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

/*	printf("\nDirections: \n");
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