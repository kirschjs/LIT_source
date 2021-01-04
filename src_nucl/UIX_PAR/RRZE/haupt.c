#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include <sys/stat.h>
#include <unistd.h>
#include <errno.h>
#include <time.h>

#include "comm.h"

// Steuerprogramm

#ifdef PC_FORTRAN
#define calcu calcu_
#define drquaak drquaak_
#endif

#ifdef CRAY_FORTRAN
#define calcu CALCU
#define drquaak DRQUAAK
#endif

#ifdef SUN_FORTRAN
#define calcu calcu_
#define drquaak drquaak_
#endif

void drquaak();
void calcu();


int main( int argc, char *argv[])
{
  int  myid, namelen, maxlen;
  char *processor_name;

  init_comm( argc, argv );
  machine_id( &myid );
  processor_name = machine_name(&namelen);


  /* setup fortran interface hitachi 
  if(hf_fint("-Fport(econv,prcntl)")!=0)
    {
      fprintf(stderr,"Error with hf_fint()\n");
      exit(1);
    }
    */


  if( myid==0) { // Chef!    
    printf("Chef auf %s.\n", processor_name);

    drquaak();

  } else {
    printf( "Sklave %d auf %s\n", myid, processor_name );
    //   nice( 10 );

    calcu();
  }
  comm_barrier();
  end_comm();
  
  exit(0);
  /* reset fortran interface hitachi
  if(hf_fend()!=0) 
    {
      fprintf(stderr,"Error with hf_fend()\n");
      exit(1);
    }*/

}

