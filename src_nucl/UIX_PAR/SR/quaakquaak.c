#include <stdio.h>
/*---------------------------------------------------------------
  Programm QUAAKQUAAK
  liest ein Quaout des parallelisierten TNIs
  und dreht die Byte-Ordnung um.

  QUAOUT Format:
  
  read(13) nrec,mfl,mfr,mkc,nzbvl,nzbvr
  do jrec=1,nrec
  read(13) NUML,NUMR,IK1,JK1,(((DM(K,L,I),K=1,ik1),L=1,jk1), I=1,2)
  end do
  --------------------------------------------------------------*/

/* --------------------------------------------------
   Kann neuerdings auch das Binaerformat der quafs_v3.f verarbeiten.
   Erkennt selbstaendig, ob es sich um NN oder NNN QUAOUT handelt
   
   write(nband1) nrec,mfl,mfr,mkc
   do 512 jrec=1,nrec
   512   WRITE(NBAND1) mfl,mfr,mkc,NUML,NUMR,IK1,JK1,(((DM(K,L,I),K=M1,M2),L=N1,N2), I=1,2)
   -------------------------------------------------- */

typedef union {
  char   c[8];
  int    i[4];
  long   h[2];
  double x;
} swap;

int TRANSSEX;
int TEST;
int DEBUG = 0;
int ZWEITEILCHEN = 0;

int rw_long( FILE *in, FILE *out, long *wert )
{
  int i;
  swap l, umgedreht;

  for (i=0;i<4;i++){
    l.c[i] = (char) fgetc(in);
    if( feof(in)!=0 || ferror(in)!=0 ){
      if( DEBUG>2 ) perror("rw_long");
      exit(9);
    }
  }
  for (i=0;i<4;i++) umgedreht.c[i]=l.c[ 3-i ];
  if( 0==TEST ) for (i=0;i<4;i++) fputc( umgedreht.c[i], out);

  if( 0==TRANSSEX ) *wert = l.h[0];
  else *wert = umgedreht.h[0];

  return 0;
}

int rw_double( FILE *in, FILE *out, double *wert )
{
  int i;
  swap l, umgedreht;
  
  for (i=0;i<8;i++){
    l.c[i] = (char) fgetc(in);
    if( feof(in)!=0 || ferror(in)!=0 ){
      if( DEBUG>2 ) perror("rw_double");
      exit(8);
    }
  }
  for (i=0;i<8;i++) umgedreht.c[i]=l.c[ 7-i ];
  if( 0==TEST ) for (i=0;i<8;i++) fputc( umgedreht.c[i], out );
  
  if( 0==TRANSSEX ) *wert = l.x;
  else *wert = umgedreht.x;

  return 0;
}

int main(int argc, char *argv[])
{
  long n, i=0, j, k, nrec;
  long numl, numr, ik1, jk1;
  int opts;
  double wert;
  swap l, umgedreht;
  FILE *in, *out;
  
  in=stdin;
  out=stdout;

  /* Kommandozeile */
  i=1;
  opts=0;
  TRANSSEX = 0;
  TEST = 0;
  while (i < argc) {
    if (strcmp(argv[i],"-i") == 0) {
      opts = opts | 1;
      i++;
      if ((in = fopen(argv[i],"r")) == NULL) {
	exit(1);
      }
    } else if (strcmp(argv[i],"-o") == 0) {
      opts = opts|2;
      i++;
      if ((out = fopen(argv[i],"w")) == NULL) {
	exit(1);
      }
    } else if (strcmp(argv[i],"-d") == 0) {
      i++;
      sscanf( argv[i], "%d", &DEBUG );
      printf( "Debug auf %d gesetzt.\n", DEBUG );
    } else if (strcmp(argv[i],"-t") == 0) {
      opts = opts|2;
      TEST = 1;
    } 
    i++;
  }

  /* Bedienungsanleitung */
  if( (opts != 3) ){
    printf("QUAAK QUAAK!\n");
    printf("Wandelt QUAF-QUAOUT von einer Architektur auf andere\n");
    printf("Kommandozeile: [-i <Original-Datei>] [-o <Ziel-Datei>]\n");

    exit(2);
  }
  
  /*
    Blocklaenge lesen, sollte 24 = 6*long sein.
    
    read(13) nrec,mfl,mfr,mkc,nzbvl,nzbvr
  */
  for (i=0;i<4;i++){
    l.c[i] = (char) fgetc(in);
    if( feof(in)!=0 || ferror(in)!=0 ) exit(7);
  }

  for (n=0;n<4;n++)  umgedreht.c[n]=l.c[ 3-n ];
  
  if( 24==l.h[0] ){
    TRANSSEX = 0;
    if( DEBUG>0 ) printf("Input-Datei von dieser Architektur, ");
  } else  if( 24==umgedreht.h[0] ){
    TRANSSEX = 1;
    if( DEBUG>0 ) printf("Input-Datei NICHT von dieser Architektur, ");
  } else if( 16==l.h[0] ){
    ZWEITEILCHEN = 1;
    TRANSSEX = 0;
    if( DEBUG>0 ) printf("NN-QUAOUT  ");
    if( DEBUG>0 ) printf("Input-Datei von dieser Architektur, ");
  } else  if( 16==umgedreht.h[0] ){
    ZWEITEILCHEN = 1;
    TRANSSEX = 1;
    if( DEBUG>0 ) printf("NN-QUAOUT  ");
    if( DEBUG>0 ) printf("Input-Datei NICHT von dieser Architektur, ");
  } else if( DEBUG>0 ){
    printf("Ich kenn mich nicht aus: %ld %ld\n",
	   l.h[0], umgedreht.h[0]);
    exit(9);
  }

  if( 0==TEST ) for(i=0;i<4;i++) fputc( umgedreht.c[i], out);

  rw_long( in, out, &nrec );
  if( DEBUG>1 ) printf( "hat %ld records\n", nrec);
  
  for(i=0;i<(ZWEITEILCHEN==0?6:4);i++){
    rw_long( in, out, &j );
    if( DEBUG>2 ) printf("%ld.",j);
  }
  if( DEBUG>0 ) printf("\n");

  /* Records */
  for( n=0; n<nrec; n++ ){
    rw_long( in, out, &k );
    if( DEBUG>1 ) printf( "Record %ld hat %ld Bytes, 4*long und %ld*doubkle\n",
			  n, k, ZWEITEILCHEN==0?(k-4*4)/8:(k-7*4)/8 );

    if( k<32 ) exit(5);
    
    if( ZWEITEILCHEN>0 ){
      rw_long( in, out, &ik1 ); /* mfl */
      rw_long( in, out, &ik1 ); /* mfr */
      rw_long( in, out, &ik1 ); /* mkc */
    }
    rw_long( in, out, &numl );
    rw_long( in, out, &numr );
    rw_long( in, out, &ik1 );
    rw_long( in, out, &jk1 );
    if( DEBUG>2 ) printf("numl %ld. numr %ld, ik1 %ld jk1 %ld\n", numl, numr, ik1, jk1);

    for (i=0;i<(ZWEITEILCHEN==0?(k-4*4)/8:(k-7*4)/8);i++){
      rw_double( in, out, &wert );
      if( DEBUG>4 ) printf("%ld %#16.15g\n", i, wert);
    }
    rw_long( in, out, &k );
  }
  
  fclose(in);
  fclose(out);

  return 0;
}
