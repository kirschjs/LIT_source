#include <time.h>

#ifdef PC_FORTRAN
#define zeit zeit_
#endif

#ifdef SUN_FORTRAN
#define zeit zeit_
#endif

#ifdef CRAY_FORTRAN
#define zeit ZEIT
#endif

/* Wall clock time seit zeit() das letzte Mal aufgerufen wurde */
void zeit( double *cpu )
{
  static double cpu_zeit=0;
  static double zeitdiff=0;

  if( 0==cpu_zeit ) cpu_zeit=1.0*time( NULL );

  if( time( NULL ) < cpu_zeit ){
    /* zeitdiff unveraendert lassen, Ueberlauf hat stattgefunden */
  } else {
    zeitdiff = fabs( 1.0*time( NULL ) - cpu_zeit);
  }

  cpu_zeit = time( NULL );

  *cpu=zeitdiff;
}
 
