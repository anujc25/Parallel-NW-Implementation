
#include "nw.h"
#include <omp.h>


using namespace std;

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

int nw(                                                          
        char*       seq_1,          /*  Needleman-Wunsch   */
        char*       seq_2,          /*  algorithm for      */
        char*       seq_1_al,       /*  global alignment   */
        char*       seq_2_al,       /*  of nt sequence.    */
        bool         prm
      )
{
        //printf("CCC");

        int  d = 2 ;                 /* gap  */

        int  L1 = strlen(seq_1);
        int  L2 = strlen(seq_2);

        // Dynamic programming matrix
        int ** F = new int * [ L2+1 ];
        for( int i = 0; i <= L2; i++ )  F[ i ] = new int [ L1 ];

        // Traceback matrix
        char ** traceback = new char * [ L2+1 ];
        for( int i = 0; i <= L2; i++ )  traceback[ i ] = new char [ L1 ];

        // Initialize traceback and F matrix (fill in first row and column)
        dpm_init( F, traceback, L1, L2, d );

        printf("CCC");

        // Create alignment
        nw_align( F, traceback, seq_1, seq_2, seq_1_al, seq_2_al, d );

        #if DEBUG
            int  L_al = strlen(seq_1_al);
            cout << "Length after alignment: " << L_al << endl;
        #endif

        //if( prm )
        //{
                cout << "\nDynamic programming matrix: " << "\n\n";
                //print_matrix( F, seq_1, seq_2 );

                cout << "\nTraceback matrix: " << "\n\n";
                //print_traceback( traceback, seq_1, seq_2 );

                cout << endl;
        //}

      /*  for( int i = 0; i <= L2; i++ )  delete F[ i ];
        delete[] F;
        for( int i = 0; i <= L2; i++ )  delete traceback[ i ];
        delete[] traceback;*/

        return  0 ;
}


void  dpm_init( int ** F, char ** traceback, int L1, int L2, int d )
{
        F[ 0 ][ 0 ] =  0 ;
        traceback[ 0 ][ 0 ] = 'n' ;

        int i=0, j=0;

        for( j = 1; j <= L1; j++ )
        {
                F[ 0 ][ j ] =  -j * d ;
                traceback[ 0 ][ j ] =  '-' ;
        }
        for( i = 1; i <= L2; i++ )
        {
                F[ i ][ 0 ] =  -i * d ;
                traceback[ i ][ 0 ] =  '|' ;
        }
}


int nw_align(                  // Needleman-Wunsch algorithm
              int **     F,
              char **    traceback,
              char*     seq_1,
              char*     seq_2,
              char*    seq_1_al,
              char*    seq_2_al,
              int        d         // Gap penalty
            )
{

        int        k = 0, x = 0, y = 0;
        int        fU, fD, fL ;
        char       ptr, nuc ;
        int        i = 0, j = 0;
        //char       temp[4000];

        const int  a =  1;   // Match
        const int  b = -1;   // Mismatch

        const int  s[ 4 ][ 4 ] = { { a, b, b, b },    /* substitution matrix */
                                   { b, a, b, b },
                                   { b, b, a, b },
                                   { b, b, b, a } } ;

        int  L1 = strlen(seq_1);
        int  L2 = strlen(seq_2);

        char  *temp;
        temp = (char *)malloc((L1+100) * sizeof(char));

        int z ;
        for( z = 2; z <= L1 + L2; z++ )
        {
            
            int mx = ( 1>z-L1 ? 1  : z-L1);
            int mn = ( L1<z-1 ? L1 : z-1 );

            #pragma omp parallel for
                for( i = mx; i <= mn ; i++ )
                {
                      nuc = seq_1[ z-i-1 ] ;

                        switch( nuc )
                        {
                                case 'A':  x = 0 ;  break ;
                                case 'C':  x = 1 ;  break ;
                                case 'G':  x = 2 ;  break ;
                                case 'T':  x = 3 ;
                        }

                        nuc = seq_2[ i-1 ] ;

                        switch( nuc )
                        {
                                case 'A':  y = 0 ;  break ;
                                case 'C':  y = 1 ;  break ;
                                case 'G':  y = 2 ;  break ;
                                case 'T':  y = 3 ;
                        }

                        fU = F[ i-1 ][ z-i ] - d ;
                        fD = F[ i-1 ][ z-i-1 ] + s[ x ][ y ] ;
                        fL = F[ i ][ z-i-1 ] - d ;

                        F[ i ][ z-i ] = max( fU, fD, fL, &ptr ) ;

                        traceback[ i ][ z-i ] =  ptr ;
                }
        }



/*
      for( i = 1; i <= L2; i++ )
        {
                for( j = 1; j <= L1; j++ )
                {
                        nuc = seq_1[ j-1 ] ;

                        switch( nuc )
                        {
                                case 'A':  x = 0 ;  break ;
                                case 'C':  x = 1 ;  break ;
                                case 'G':  x = 2 ;  break ;
                                case 'T':  x = 3 ;
                        }

                        nuc = seq_2[ i-1 ] ;

                        switch( nuc )
                        {
                                case 'A':  y = 0 ;  break ;
                                case 'C':  y = 1 ;  break ;
                                case 'G':  y = 2 ;  break ;
                                case 'T':  y = 3 ;
                        }

                        fU = F[ i-1 ][ j ] - d ;
                        fD = F[ i-1 ][ j-1 ] + s[ x ][ y ] ;
                        fL = F[ i ][ j-1 ] - d ;

                        F[ i ][ j ] = max( fU, fD, fL, &ptr ) ;

                        traceback[ i ][ j ] =  ptr ;
                }
        }
*/
        i=L2;
        j=L1;

        //i-- ; j-- ;

        while( i > 0 || j > 0 )
        {
                switch( traceback[ i ][ j ] )
                {
                        case '|' :      strcat(seq_1_al,"-");
                                        sprintf(temp,"%c",seq_2[i-1]);
                                        strcat(seq_2_al,temp);
                                        i-- ;
                                        break ;

                        case '\\':      sprintf(temp,"%c",seq_1[j-1]);
                                        strcat(seq_1_al,temp);
                                        sprintf(temp,"%c",seq_2[i-1]);                                            
                                        strcat(seq_2_al,temp);
                                        i-- ;  j-- ;
                                        break ;

                        case '-' :      sprintf(temp,"%c",seq_1[j-1]);
                                        strcat(seq_1_al,temp);
                                        strcat(seq_2_al,"-");
                                        j-- ;
                }
                k++ ;
        }
        //char temp[20];
        strcpy(temp,strrev(seq_1_al));
        strcpy(seq_1_al,temp);
        strcpy(temp,strrev(seq_2_al));
        strcpy(seq_2_al,temp);
        //reverse( seq_1_al.begin(), seq_1_al.end() );
        //reverse( seq_2_al.begin(), seq_2_al.end() );

        return  0 ;
}


int  max( int f1, int f2, int f3, char * ptr )
{
        int  max = 0 ;

        if( f1 >= f2 && f1 >= f3 )  
        {
                max = f1 ;
                *ptr = '|' ;
        }
        else if( f2 > f3 )              
        {
                max = f2 ;
                *ptr = '\\' ;
        }
        else
        {
                max = f3 ;
                *ptr = '-' ;
        }
        
        return  max ;   
}


void  print_matrix( int ** F, char* seq_1, char* seq_2 )
{
        int  L1 = strlen(seq_1);
        int  L2 = strlen(seq_2);

        cout << "        ";
        for( int j = 0; j < L1; j++ )
        {
                cout << seq_1[ j ] << "   ";
        }
        cout << "\n  ";

        for( int i = 0; i <= L2; i++ )
        {
                if( i > 0 )
                {
                        cout << seq_2[ i-1 ] << " ";
                }
                for( int j = 0; j <= L1; j++ )
                {
                        cout.width( 3 );
                        cout << F[ i ][ j ] << " ";
                }
                cout << endl;
        }
}


void  print_traceback( char ** traceback, char* seq_1, char* seq_2 )
{
        int  L1 = strlen(seq_1);
        int  L2 = strlen(seq_2);

        cout << "    ";
        for( int j = 0; j < L1; j++ )
        {
                cout << seq_1[ j ] << " ";
        }
        cout << "\n  ";

        for( int i = 0; i <= L2; i++ )
        {
                if( i > 0 )
                {
                        cout << seq_2[ i-1 ] << " ";
                }
                for( int j = 0; j <= L1; j++ )
                {
                        cout << traceback[ i ][ j ] << " ";
                }
                cout << endl;
        }
}


void  print_al( char* seq_1_al, char* seq_2_al )
{
        cout << seq_1_al << endl;
        cout << seq_2_al << endl;
}


