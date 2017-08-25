#include <iostream>
#include <string.h>
#include <algorithm>
#include <stdlib.h>
#include <stdio.h>

#define DEBUG 0

using namespace std;

extern int  nw( 
                 char*, char*, 
                 char*, char*,
                 bool
              );

extern int  nw_align( 
                      int **, char **,
                      char*, char*, 
                      char*, char*,
                      int 
                    );

extern void  dpm_init        ( int **, char **, int, int, int );
extern void  print_al        ( char*, char* );
extern void  print_matrix    ( int ** const, char*, char* );
extern void  print_traceback ( char ** const, char*, char* );
extern int   max             ( int, int, int, char * );


