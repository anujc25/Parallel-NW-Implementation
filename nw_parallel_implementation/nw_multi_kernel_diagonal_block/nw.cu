/* Needleman Wunsch Algorithm Parallel Implementation */ 

#include "nw.h"
#include <fstream>


int main(int argc , char **argv)
{


    //size_t len;
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

            seq_1 = (char *) malloc( 2 *size1 * sizeof(char));
            
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
            seq_2 = (char *) malloc( 2* size2 * sizeof(char));
            
            file2 >> seq_2 ;
            
            file2.close();
            printf("Seq 2: %s size2: %d\n",seq_2,size2);
        }
        

        // Aligned sequences
        char  *seq_1_al;
        char  *seq_2_al;

    seq_1_al = (char *) malloc( 2 * size1 * sizeof(char));
    seq_2_al = (char *) malloc( 2 * size2 * sizeof(char));


    strcpy(seq_1_al,"");
    strcpy(seq_2_al,"");




/*        if( argc < 3 )
        {
                cout << "\n   Usage:   " << "File" << " seq_1 seq_2 [-p]\n";
                cout << "\n   -p:       Print matrices\n";
                cout << "\n   Output:   alignment\n\n";
                exit( 1 ) ;
        }                       
        
        // Sequences to be aligned
        char  *seq_1;  // =  argv[ 1 ] ;
        char  *seq_2;  // =  argv[ 2 ] ;

	seq_1 = (char *) malloc( strlen(argv[1]) * sizeof(char));
	seq_2 = (char *) malloc( strlen(argv[2]) * sizeof(char));

	strcpy(seq_1,argv[1]);
	strcpy(seq_2,argv[2]);

        // Aligned sequences
        char  *seq_1_al;
        char  *seq_2_al;

	seq_1_al = (char *) malloc( 3 * strlen(argv[1]) * sizeof(char));
	seq_2_al = (char *) malloc( 3 * strlen(argv[2]) * sizeof(char));

    strcpy(seq_1_al,"");
    strcpy(seq_2_al,"");
*/

    struct timespec t1,t2; double dt1;
    clock_gettime(CLOCK_REALTIME,  &t1);


    // Get alignment
    nw( seq_1, seq_2, seq_1_al, seq_2_al) ;   


    clock_gettime(CLOCK_REALTIME,  &t2);
    dt1 = (t2.tv_sec - t1.tv_sec) + (double) (t2.tv_nsec - t1.tv_nsec) * 1e-9;
    double time=dt1*1000 ;

    printf("\nOriginal Sequences:");
    printf("\n\t\t\t%s \n\t\t\t%s\n",seq_1,seq_2);

    printf("\nAfter Alignment:");
    printf("\nSeq1: \t\t\t%s \n\nSeq2:\t\t\t%s\n\n",seq_1_al,seq_2_al);

    printf("\n%10f  kernel Time elapsed with only threads\n", time);



return 0;
}
