#include <fstream>
#include <time.h>
#include "nw.h"

using namespace std;


int  main( int argc, char ** argv )
{
    
        bool    prm = false;

/*        if( argc < 3 )
        {
                cerr << "\n   Usage:   " << "File" << " seq_1 seq_2 [-p]\n";
                cerr << "\n   -p:       Print matrices\n";
                cerr << "\n   Output:   alignment\n\n";
                exit( 1 ) ;
        }                       
*/            

        size_t len;
        // Sequences to be aligned
        char  *seq_1;  // =  argv[ 1 ] ;
        char  *seq_2;  // =  argv[ 2 ] ;

        //FILE  *file1 , *file2 ;
        int   size1 , size2 ;
        ifstream file1("../s1.txt");
        if (file1==NULL)
             perror ("Error opening file 1");
        else
        {
            file1.seekg (0, ios::end);
            size1 = file1.tellg();
            file1.seekg (0, ios::beg);

            seq_1 = (char *) malloc( 2 * size1 * sizeof(char));
            
            file1 >> seq_1 ;
            
            file1.close();
            printf("Seq 1: %s size1: %d\n",seq_1,size1);
        }
        

        ifstream file2("../s2.txt");
        if (file2==NULL)
             perror ("Error opening file 2");
        else
        {
            file2.seekg (0, ios::end);
            size2 = file2.tellg();
            file2.seekg (0, ios::beg);
            seq_2 = (char *) malloc(  2 * size2 * sizeof(char));
            
            file2 >> seq_2 ;
            
            file2.close();
            printf("Seq 2: %s size2: %d\n",seq_2,size2);
        }
        
        printf("AAAAAAAAA  ");


	//strcpy(seq_1,argv[1]);
	//strcpy(seq_2,argv[2]);

    
        //if( argc == 4 )
        //{
        //        string  pr_arg  =  argv[ 3 ] ;
        //        if( pr_arg == "-p" )  
        prm = true;   // Print matrices
        //}                       

        // Aligned sequences
        char  *seq_1_al;
        char  *seq_2_al;

	seq_1_al = (char *) malloc( 2 * size1 * sizeof(char));
	seq_2_al = (char *) malloc( 2 * size2 * sizeof(char));

    strcpy(seq_1_al,"");
    strcpy(seq_2_al,"");



    struct timespec t1,t2; double dt1;
    clock_gettime(CLOCK_REALTIME,  &t1);

        // Get alignment
        nw( seq_1, seq_2, seq_1_al, seq_2_al, prm ) ;   


    clock_gettime(CLOCK_REALTIME,  &t2);
    dt1 = (t2.tv_sec - t1.tv_sec)  + (double) (t2.tv_nsec - t1.tv_nsec) * 1e-9  ;
    double time=dt1*1000 ;

    printf("\n%s \n%s \n",seq_1_al,seq_2_al);



    printf("\n\n%10f filltightstring kernel Time elapsed \n", time);    


       return  0 ;
}

























