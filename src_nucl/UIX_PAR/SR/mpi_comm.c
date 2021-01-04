#include <stdio.h>
#include "comm.h"
#include <mpi.h>


MPI_Status *status=NULL;

#ifdef DEBUG_INT
FILE *FD;
#endif

/* Stellt die vom gen. Algorithmus benötigten Kommunikations-
   formen über MPI dar.

   MPI KANN DIE VIRTUELLE MASCHINE NICHT ZUR LAUFZEIT AENDERN.
 */

/* Send data
   Zum Empfang muss data schon vorher alloziert sein.
   dest==0 soll immer der MASTER sein !!!!
 */
void send_double( double *data, int *len, int *dest){
#ifdef DEBUG_DOUBLE
  int i;
  fprintf( FD, "Diuble send to %d: %d values\n", *dest, *len );
  for( i=0; i<*len; i++ ) fprintf( FD, "%f ", data[i] );
#endif
  if (NULL == status) status=(MPI_Status*)malloc(sizeof(MPI_Status));
  MPI_Ssend( data, *len, MPI_DOUBLE, *dest, 0, MPI_COMM_WORLD); 
}

void send_int( int *data, int *len, int *dest){
#ifdef DEBUG_INT
  int i,j;
  fprintf( FD, "Send to %d: %d Bytes\n", *dest, *len );
  j=0;
  for( i=0; i<*len; i++ ){
    j++;
    if( j%25 == 0 ){
      fprintf( FD, "\n" );
      j=0;
    }
    fprintf( FD, "%d ", data[i] );
  }
  fprintf( FD, "\n" );
#endif
  if (NULL == status) status=(MPI_Status*)malloc(sizeof(MPI_Status));
      MPI_Ssend( data, *len, MPI_INT, *dest, 0, MPI_COMM_WORLD); 
}

void send_long( long *data, int *len, int *dest){
  if (NULL == status) status=(MPI_Status*)malloc(sizeof(MPI_Status));
      MPI_Ssend( data, *len, MPI_LONG, *dest, 0, MPI_COMM_WORLD); 
}

/* Receive data */
void receive_double( double *data, int *len, int *source, int *flag){
#ifdef DEBUG_DOUBLE
  int i;
#endif
  if (NULL == status) status=(MPI_Status*)malloc(sizeof(MPI_Status));
  if ( 0==*flag ) {
    MPI_Recv( data, *len, MPI_DOUBLE,*source,0,MPI_COMM_WORLD,status);
#ifdef DEBUG_DOUBLE
  fprintf( FD, "Double rec. from %d: %d values\n", *source, *len );
  for( i=0; i<*len; i++ ) fprintf( FD, "%f ", data[i] );
#endif
    return;
  }
  // Auf jeden hoeren, source entsprechend setzen
  if( 1==*flag ){
    MPI_Recv( data, *len, MPI_DOUBLE, MPI_ANY_SOURCE, MPI_ANY_TAG,
	      MPI_COMM_WORLD,status);
    *source = status->MPI_SOURCE;
    return;
  }
}

void receive_int( int *data, int *len, int *source, int *flag){
#ifdef DEBUG_INT
  int i,j;
#endif
  if (NULL == status) status=(MPI_Status*)malloc(sizeof(MPI_Status));
  if ( 0==*flag ) {
    MPI_Recv( data, *len, MPI_INT,*source,0,MPI_COMM_WORLD,status);
#ifdef DEBUG_INT
  fprintf( FD , "Receive from %d: %d Bytes\n", *source, *len );
  j=0;
  for( i=0; i<*len; i++ ){
    j++;
    if( j%25 == 0 ){
      fprintf( FD, "\n" );
      j=0;
    }
    fprintf( FD, "%d ", data[i] );
  }
  fprintf( FD, "\n" );
#endif
    return;
  } else {
    MPI_Recv( data, *len, MPI_INT, MPI_ANY_SOURCE, 0, MPI_COMM_WORLD, status );
    *source = status->MPI_SOURCE;
  }
}

void receive_long( long *data, int *len, int *source, int *flag){
  if (NULL == status) status=(MPI_Status*)malloc(sizeof(MPI_Status));
  if ( 0==*flag ) {
    MPI_Recv( data, *len, MPI_LONG, *source, 0, MPI_COMM_WORLD, status );
    return;
  }
}

/* Get information on virtual machine */
void machine_size( int *id ){
  MPI_Comm_size( MPI_COMM_WORLD, id );
}

void machine_id( int *myid ){
  MPI_Comm_rank( MPI_COMM_WORLD, myid );
}

char* machine_name(int *len){
  char *processor_name;
  processor_name = (char*)calloc(MPI_MAX_PROCESSOR_NAME+1,sizeof(char));
  MPI_Get_processor_name(processor_name,len);
  return processor_name;

}

/* Configure virtual machine
   (init, abort, shutdown) */
void init_comm( int argc, char *argv[]){
  char *name;
  int myid;
  name = (char*)calloc(128,sizeof(char));
  
  MPI_Init(&argc,&argv);
  MPI_Comm_rank( MPI_COMM_WORLD, &myid );
  
  strcpy( name, "INT_DUMP." );
  sprintf( name+9, "%04d", myid );
  
  #ifdef DEBUG_INT
  FD = fopen( name, "w" );
  #endif
}

void abort_comm( int code ){
  MPI_Abort( MPI_COMM_WORLD, code );
}

void end_comm(){
  MPI_Finalize();
#ifdef DEBUG_INT  
  fclose(FD);
#endif
}

void comm_barrier(){
  MPI_Barrier( MPI_COMM_WORLD );
}
