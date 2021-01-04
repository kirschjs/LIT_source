      PROGRAM QUAL
      IMPLICIT double precision (A-H,O-Z)
      include 'mpif.h'
C
C
C     FUER ELEKTROMAGNETISCHE OPERATOREN OHNE LANGWELLE
C
C                KEINE POLYNOME
C
C     VERSION MIT LOOP'S BASISVEKTOR*RADIALPARAMETER ZUM VERARBEITEN
C     GROSSER MATRIZEN
C
C     BEI DIESER VERSION KOENNEN MIT NBAND5#0 TEILERGEBNISSE
C     ANDERER QUAL-RECHNUNGEN UEBERNOMMEN WERDEN.
C
C
C     VERSION APRIL 88
C     AENDERUNG 27.06.88: KONTROLLE FUER NDIM2 EINGEFUEGT UND DIMENSIONIERUNG
C                         GEAENDERT M.U.
C     DIMENSIONIERUNG VON NDIM1 GEAENDERT 06.07.88 M.U.
C     AENDERUNG 15.09.88: OPEN-STATEMENTS EINGEFUEGT, VARIABELE BEIM
C                         SCHREIBEN/LESEN VON NBAND5 UND NBAND7 GEAEN-
C                         DERT UND KORRIGIERT  M.U.
C     AENDERUNG 13.12.88: FALSCHE VORBESETZUNG DES FELDES IND KORRIGIERT M.U.
C     AENDERUNG 08.06.89: FALSCHE DIMENSIONIERUNG VON MVK IN MAT KORRIGIERT
C                         M U.
C     AENDERUNG 25.10.90: DIMENSIONIERUNG GEAENDERT  M.U.
C     LETZTE AENDERUNG 25.10.90  M.U.
C
C     INTEGERS IM FORMAT 24I3, REALS IM FORMAT 6E12.4
C
C
C     WENN NBAND5=0, WIRD ALLES  NEU GERECHNET; WENN NBAND5=12, WIRD VON
C     DIESEM EIN TEIL UBERNOMMEN
C
C
      INCLUDE 'par/QUAL'
C     NZOPER: ANZAHL DER OPERATOREN IN QUAL
C     NZOPOB:   "     "       "      " OBEM
C     NZOPLU:   "     "       "      " LUISE
C     NZTMAX: MAXIMALE ANZAHL DER TEILCHEN
C     NZFMAX:    "       "     "  ZERLEGUNGEN
C     NZCMAX:    "       "     "  CLUSTER
C     MZGMAX:    "       "     "  SPINFUNKTIONEN
C     NZLWMA:    "       "     "  DREHIMPULSSTRUKTUREN
C     NZRHOM:    "       "     "  BASISVEKTOREN EINER ZERLEGUNG
C     MZPARM:    "       "     "  RADIALPARAMETER
C     NZPARM:    "       "     "  SAETZE INNERER WEITEN
C     NZBMAX:    "       "     "  BASISVEKTOREN UEBER ALLE ZERLEGUNGEN
C     NZPOMA:    "       "     "  POLYNOMSTRUKTUREN
C     NDIM  :    "       "    BASISVEKTOREN * RADIALPARAMETER
C             NDIM .LE. NZBMAX * MZPARM
C
      integer i, info, nproc, nhost
      INTEGER iStatus( MPI_STATUS_SIZE )
      integer myid
      integer speed
      character*18 nodename, host
      character*8 arch
      character*180 name, fnumber, qname, oname, cquaf
      character*180 cquaout, coutput,slaveout
      CHARACTER*1 dummy(80)
      character*1 cdummy
      integer nteil

      INTEGER DATE_TIME (8)
      CHARACTER (LEN = 10) BIG_BEN (3)
C
      parameter(lbeg = 1 + 1 + 1 + 1,
     * lstore = nzparm*nzrhom*nzfmax + 2*nztmax*nztmax +
     *         mzgmax*mzgmax*npdc +
     *         nztmax*(nzcmax) + nzvmax*(2*nzcmax-1) +
     *         mzparm*nzfmax + (nzavma)*nzparm*nzfmax,
     * lluri = ndim1 + ndim1 + ndim1 + ndim1*nziqma + 1 + 1 + 1 +
     *         2*2*(nzavma) + 1 + 1,
     * llurr = ndim1 + ndim4*ndim4,
     * limp = ndim4*ndim4 + 2*nzbmax*nzbmax + 1 + 1 + 1 + 1 + 1 +1 +
     *                                      1 + 1 + 1 + 1,
     * lkomi = 4*nzrhom*nzfmax + nzfmax + nzfmax + nzfmax + nzfmax +
     *        nztmax*npdc + 2*nzcmax*nzfmax +2*ndim5 + npdc + 5,
     * lindy = nzrhom*nzrhom,
     * lbig = ndim*ndim*ndim2 + ndim*ndim,
     * lcsh = ndim5,
     * llup2 = 3+ndim1+ndim1+ndim1*nzavma,
     * lfltw = ndim4*ndim4*ndim2)
C
      PARAMETER (LIND=(NZFMAX*(NZFMAX+1))/2*NZOPER)
      DIMENSION INDEX(0:LIND)
C
      COMMON /STORE/ COF(NZPARM,NZRHOM,NZFMAX),RVEC(NZTMAX,NZTMAX),
     *    SVEC(NZTMAX,NZTMAX), U(MZGMAX,MZGMAX,NPDC),
     *    VW(NZTMAX,NZCMAX),
     *    QFCL(NZVMAX,2*NZCMAX-1),
     *    RPAR(MZPARM,NZFMAX),CPAR(NZAVMA,NZPARM,NZFMAX)
      dimension rstore(lstore)
      equivalence (rstore(1), cof)
C
      COMMON /BEG/ IZLURE, NBAUS, JWIED, NZAVZU
      dimension ibeg(lbeg)
      equivalence (ibeg(1), izlure)
C  commoners ----------------------------------------------
      COMMON /LURI/ KVK(NDIM1), JQ(NDIM1), INDPO(NDIM1),
     *             MVK(NZIQMA,NDIM1),  NDIM3,
     *             IZ, JZ, KZHV(2,2*NZAVMA), LUPAUS, IQM
      dimension iluri(lluri)
      equivalence (iluri(1), kvk)
C
      COMMON /LURR/ EPO(NDIM1),WERTL(NDIM4,NDIM4)
      dimension rlurr(llurr)
      equivalence (rlurr(1), epo)
C
      COMMON /IMP/ MKC,LAMB(NDIM4*NDIM4),LAHI(NZBMAX*NZBMAX),
     *       NUHI(NZBMAX*NZBMAX),
     *       MD1, MFL, NCC, NREB, NZAV, NZV1,
     *       LPARR, NZV,MC1
      dimension iimp(limp)
      equivalence (iimp(1), mkc)
C
      COMMON /KOMIX/ NUM(4,NZRHOM,NZFMAX), NZC(NZFMAX),
     *          NZPAR(NZFMAX), MZPAR(NZFMAX),NZRHO(NZFMAX),
     *          LC(NZTMAX,NPDC),NGRV(2,NZCMAX,NZFMAX),
     *          NSH1(NDIM5), NSH2(NDIM5), NT1(NPDC),
     *          ITV2, NTE, NZT, MFR, NSH
      dimension ikomi(lkomi)
      equivalence (ikomi(1), num)
C
      COMMON /LUR/ IENT(NZOPLU), NDIM1H, NBAND3
C
      COMMON /INDY/ IND(NZRHOM,NZRHOM)
      dimension iindy(lindy)
      equivalence (iindy(1),ind)

      COMMON /BIG/ F(NDIM,NDIM),DM(NDIM,NDIM,NDIM2) 
      dimension rbig(lbig)
      equivalence (rbig(1),f)

C
C  MULPOL common blocks -----------------------------------
      COMMON /LUP2/ NDIM4H,MUL,LINDEX(NDIM1),LINT(NZAVMA,NDIM1),
     *              LSUM(NDIM1), NDIM2H
      dimension ilup2(llup2)
      equivalence (ilup2(1),NDIM4H)
C
      COMMON /CSH/ SH(NDIM5)
      dimension rcsh(lcsh)
      equivalence (rcsh(1),sh)     
C
      COMMON /FLTW/ WERT(NDIM4,NDIM4,NDIM2)
      dimension rfltw(lfltw)
      equivalence (rfltw(1),wert)      
C ---------------------------------------------------------
      DIMENSION NREG(NZOPOB),
     *          MZG(NZFMAX),NOL(NZFMAX),
     *          MS(MZGMAX,NZFMAX)
      DIMENSION LREG(NZOPER), KOM(NZOPER,NZFMAX,NZFMAX)
      DIMENSION kom2(NZFMAX,NZFMAX)
      DIMENSION MMASSE(2,MZGMAX,NZFMAX), MLAD(2,MZGMAX,NZFMAX),
     *          MSS(2,MZGMAX,NZFMAX)
      DIMENSION NZLW(NZFMAX), LW(2*NZCMAX-3,NZLWMA,NZFMAX),
     *          LZWERT(5,NZLWMA,NZFMAX), 
     *          NZPO(NZFMAX),KP(NZCMAX-1,NZPOMA,NZFMAX),
     *          ITPO(NZFMAX),
     *          NZC1(NZFMAX)
      DIMENSION LREG1(NZOPER), NZRHO1(NZFMAX)
      DIMENSION IMV(NZMAT),UM(NZMAT,NZMAT),
     *          LCALT(NZTMAX),NKSU(NZFMAX)
C
      CHARACTER*80 INFIL  ,INFILU
C     
      character*(MPI_MAX_PROCESSOR_NAME) my_name
      integer my_len , error

      integer jdummy(1 : 80)

C     MPI initialisieren
      call MPI_INIT(ierr)
      call MPI_COMM_RANK(MPI_COMM_WORLD,myid,ierr)
      call MPI_GET_PROCESSOR_NAME( my_name , my_len ,error )
      
C     ######################################################################
C     Der gesamte Master wird nur fuer myid .eq. 0 durchlaufen,
C     sonst muss in den slave gesprungen werden
C     ######################################################################   
      if (myid .eq. 0) then

C     Anzahl der Prozessoren
         call MPI_COMM_SIZE(MPI_COMM_WORLD,nproc,ierr)


      OPEN(UNIT=3,STATUS='SCRATCH',FORM='UNFORMATTED')
      OPEN(UNIT=5,FILE='INQUA',STATUS='OLD')
      OPEN(UNIT=6,FILE='OUTPUT')
      OPEN(UNIT=8,FILE='OBOUT',STATUS='OLD',FORM='UNFORMATTED')
      OPEN(UNIT=9,FILE='LUOUT',STATUS='OLD',FORM='UNFORMATTED')
      OPEN(UNIT=10,FILE='QUAOUT',STATUS='UNKNOWN',FORM='UNFORMATTED')
      OPEN(UNIT=12,FILE='QUAALT',STATUS='UNKNOWN',FORM='UNFORMATTED')
C
      INPUT=5
      NOUT=6
      write(NOUT,*) 'Anzahl der Prozessoren: ',nproc
      WRITE (NOUT,1000)
1000  FORMAT('1',//,1X,'PROGRAMM QUAL_parallel',/,1X,13('='),////,
     *       1X,'VERSION VOM 27.05.20')
C  
      NDIM1H=NDIM1
      NDIM2H=NDIM2
      NDIM3=NZLWMA
      NDIM4H=NDIM4
      PHI=3.1415926536
      PHI2 =  SQRT(PHI)
      NQ = NZBMAX * NZBMAX
      DO 10 I = 1,NQ
   10 LAHI(I) = 1
      READ (INPUT,1002) NBAND1,NBAND2,NBAND3,NBAND4,NBAND5,
     1 naufset,NAUS,MOBAUS,LUPAUS,NBAUS,NLES
C
      if(naufset.eq.1) stop 'naufset=1'
      open(unit=18,file='FINDEX',status='unknown',form='unformatted')
      if(naufset.eq.0) then
         naufset=1
         nband5=0
         do i=1,lind
           index(i)=-1
         enddo
      else
         nband5=12
            if(nles.eq.0) then
      read(18) nger,(index(i),i=1,nger)
            else
                READ(INPUT,6666) NLES
6666       FORMAT(I6)
                read(18) nger,(index(i),i=1,nles)
           write(nout,*) ' Von ',nger,' gerechneten operatoren werden',
     *          nles,' eingelesen'
                nger=nles
            endif
      write(nout,*) 'nger',nger,index(nger)
c     write(6,*) 'nger',nger,(index(i),i=1,nger)
      do 59 i=1,nger
      if(index(i).lt.naufset) goto 59
      write(6,*) ' Fuer ',i,' ist index ',index(i),
     *   'nicht groesser naufset'
      stop 'naufset falsch'
59    continue
      do i=nger+1,lind
           index(i)=-1
      enddo
      endif
C
C     NBAND1  TAPE FUER ERGEBNISSE
C     NBAND2  TAPE VON OBEM
C     NBAND3  TAPE VON LUISE
C     NBAND4  TAPE ZUM ZWISCHENSPEICHERN DER SPIN-ISOSPIN-FUNKT.
C     NBAND5  TAPE ZUM KOPIEREN VON TEILERGEBNISSEN
C     NAUS    STEUERT ZUSAETZLICHEN AUSDRUCK IN QUAL
C     MOBAUS  ERGEBNISSE AUS OBEM WERDEN AUSGEDRUCKT
C     LUPAUS       "      "  LUISE WERDEN AUSGREDRUCKT; BEI LUPAUS
C             > 2 ZUSAETZLICHER AUSDRUCK IN MAT
C     NBAUS   ZUSAETZLICHER AUSDRUCK IN BEGRI

C
      READ(INPUT,1002)(LREG(K),K=1,NZOPER)
C
C
C     EINLESEN DER FUNKTIONSEIGENSCHAFTEN UND PARAMETER
C
      REWIND  NBAND2
      READ(NBAND2) NZF,NZT,NZV,(NREG(K),K=1,NZOPOB), INDEPO
      IF(INDEPO.NE.0) WRITE(NOUT,*)' OBERMATRIXELEMENTE VON EINZELFILES'
C
      do kl=1,nzfmax
         do kr=1,kl
            kom2(kl,kr)=1
         enddo
      enddo
C
      WRITE(NOUT,6100)
      DO 94 KL=1,NZF
      DO 94 KR=1,KL
      DO 94 L=1,NZOPER
      I=I+1
      INDEX(I)=-1
   94 KOM(L,KL,KR)=0
   93 CONTINUE
      NZV1 = NZV - 1
      DO 20  K = 1,NZF
      READ (NBAND2) NZC(K),MZG(K),NOL(K)
C     NZC     ANZAHL DER CLUSTER
C     MZG       "     "  SPIN-ISOSPIN-FUNKTIONEN
C     NOL     ANZAHL DER CLUSTER IM 1. FRAGMENT; WIRD WEITER UNTEN
C             UM EINS ERHOEHT!!!!
C
      M = MZG(K)
      IF(M.LE.MZGMAX) GOTO 6002
      N1=1
      N2=M
      N3=MZGMAX
      call fehler(6040)
6002  NOL(K)=NOL(K)+1
      READ (NBAND2) ((MMASSE(N,L,K),MLAD(N,L,K),MSS(N,L,K),N=1,2),
     1     MS(L,K),L=1, M)
C     MMASSE(I,.,.)  MASSE DES I-TEN FRAGMENTS IN NUKLEONENZAHLEN
C     MLAD(I,.,.)    LADUNG DES I-TEN FRAGMENTS IN ELEMENTARLADUNGEN
C     MSS(I,.,.)     SPIN DES I-TEN FRAGMENTS
C     MS(.,.)        GESAMTSPIN
C 
      M = NZC(K)
      READ (NBAND2) ((NGRV(N,L,K),L=1,M),N=1,2)
C     NGRV(1,L,K)  ANZAHL DER TEILCHEN IN DEN CLUSTERN 1 BIS L-1
C                  PLUS 1
C     NGRV(2,L,K)  ANZAHL DER TEILCHEN IM CLUSTER L
C
   20 CONTINUE
      REWIND NBAND3
      READ(NBAND3) NZF1,MUL,(NZLW(K),NZC1(K),K=1,NZF1),
     1(IENT(K),K=1,NZOPLU), (NZPO(MH),MH=1,NZF1),INDELU
      IF(INDELU.NE.0) WRITE(NOUT,*)' LUDWMATRIXELEMENTE VON EINZELFILES'     
      IF(NAUS.GT.0) WRITE(NOUT,248) NZF1,MUL,(IENT(J),J=1,NZOPLU),
     1(NZLW(K),NZC1(K),NZPO(K),K=1,NZF1)
248   FORMAT(' VON LUDI,ZAHL DER ZERLEGUNGEN',I3,'  MULTIPOLARITAET',
     1I3,'     OPERATOREN',7I4,/,('     ZAHL DER DREHIMPULSE',I4,
     2'     ZAHL DER CLUSTER',I4,'     ZAHL DER POLYNOME',I4,/))
      IF(NZF.EQ.NZF1) GOTO 251
       write (nout,*) 'nzf, nzf1',nzf,nzf1
      call fehler(251)
  250 call fehler(250)      
C
  251 DO 252 K=1,NZF
      ITPO(K)=0
      IF(NZPO(K).NE.0) ITPO(K)=2
      NZPO(K)=MAX0(1,NZPO(K))
      IF(NZC(K).NE.NZC1(K)) call fehler( 252)
      IF(NZPO(K).GT.NZPOMA) call fehler( 252)
  252 CONTINUE
C
C     UEBERPRUEFUNG, OB DATEN AUS LUISE(IENT) UND OBEM(NREG) ZUSAMMENPASSEN
      DO 254, K=1,NZOPER
      KH=INT(K/2)+1
      IF (LREG(K).LE.0) GOTO 254
      write(nout,*) LREG(K),NREG(K),IENT(KH)
      IF (LREG(K).GT.NREG(K)) call fehler( 254)
      IF (LREG(K).GT.IENT(KH)) call fehler( 254)
254   CONTINUE
C     ENDE UEBERPRUEFUNG
C
C     KONSTRUKTION DER BAHNDREHIMPULSE
C
C     BESTIMMUNG DER LZWERTE(.,.,.) FUER ALLE DREHIMPULSSTRUKTUREN
C     LZWERT(1,.,.): BAHNDREHIMPULS DES 1. FRAGMENTS
C     LZWERT(2,.,.):       "         "  2.     "
C     LZWERT(3,.,.): AUS LZWERT(1,.,.) UND LZWERT(2,.,.) GEKOPPELTER DREHIMPULS
C     LZWERT(4,.,.): RELATIVDREHIMPULS ZWISCHEN DEN FRAGMENTEN
C     LZWERT(5,.,.): GESAMTBAHNDREHIMPULS
C
      DO 264 K=1,NZF
      M1=2*NZC(K)-3
      M2=NZLW(K)   
      IF(M2.LE.NDIM3) GOTO 6000
      N1=2
      N2=M2
      N3=NDIM3
      call fehler(6040)
 6000 CONTINUE
      M3=NZC(K)-1
      MPO=NZPO(K)
      READ(NBAND3) ((LW(M,L,K),M=1,M1),L=1,M2),((KP(M,KH,K),M=1,M3),
     1 KH=1,MPO)
      DO 265 L=1,M2
       DO 266 M=1,3
  266  LZWERT(M,L,K)=0
C
       LZWERT(4,L,K)=LW(M3,L,K)
       LZWERT(5,L,K)=LW(M1,L,K)
C
       IF(M3-2) 265,267,268
C
267    CONTINUE
C      3-CLUSTER-FALL
       LZWERT(3,L,K)=LW(1,L,K)
       IF(NOL(K).GE.NZC(K)) GOTO 271
       LZWERT(2,L,K)=LW(1,L,K)
       GO TO 265
  271  LZWERT(1,L,K)=LW(1,L,K)
       GO TO 265
C
268    CONTINUE
C      4-CLUSTER-FALL
       IF(NOL(K).LT.NZC(K)) GOTO 272
       LZWERT(3,L,K)=LW(4,L,K)
       LZWERT(1,L,K)=LZWERT(3,L,K)
       GO TO 265
  272  IF(NOL(K).GT.2) GOTO 275
       LZWERT(2,L,K)=LW(4,L,K)
       LZWERT(3,L,K)=LZWERT(2,L,K)
       GO TO 265
  275  LZWERT(3,L,K)=LW(4,L,K)
       LZWERT(1,L,K)=LW(1,L,K)
       LZWERT(2,L,K)=LW(2,L,K)
C
  265  CONTINUE
C
  264 CONTINUE
C     ENDE DER KONSTRUKTION DER BAHNDREHIMPULSE
C
C
C
C     KONSTRUKTION DER BASISVEKTOREN
C
      I=0
      WRITE (NOUT,1006)
 1006 FORMAT(1H0)
      DO 22  K = 1,NZF
       NZPAR(K) = 0
       MZPAR(K) = 0
C
       READ(INPUT,1002) NZRHO(K)
C      ANZAHL DER BASISVEKTOREN DER ZERLEGUNG
C
       KK=NZRHO(K)
       IF(KK.LT.NZRHOM) GOTO 388
       N1=8
       N2=KK
       N3=NZRHOM
       call fehler(6040)
388    CONTINUE
       IF(KK.EQ.0) GOTO 22
       M = 2*NZC(K) - 2
       WRITE (NOUT,1001)  K
 1001  FORMAT(//25H WEITEPARAMETER FUER DIE ,  I3,
     1        47H TE CLUSTERSTRUKTUR IN REZIPROKEN FERMIQUADRAT     )
C 
       READ (INPUT,1002) NZPAR(K), MZPAR(K)
c        write(nout,*) K,'/',NZF
C      NZPAR: ANZAHL DER SAETZTE INNERER WEITEN
C      MZPAR: ANZAHL DER RADIALWEITEN
C
       LM = NZPAR(K)
       KM=MZPAR(K)
       IF(LM.LE.NZPARM) GOTO 390
       N1=7
       N2=LM
       N3=NZPARM
       call fehler(6040)
390    DO 24  L = 1,LM
        READ (INPUT,1003)  (CPAR(N,L,K),N=1,M)
C       CPAR(N,L,K): N-TE INNERE WEITE DES L-TEN SATZES AN INNEREN
C                    WEITEN DER K-TEN ZERLEGUNG
C
 1003  FORMAT(6E12.4)
        WRITE (NOUT,1004)L,(  CPAR(N,L,K),N=1,M)
24      CONTINUE
 1004  FORMAT(/I3,25H TER SATZ INNERER WEITEN         /8F12.6)
C
       READ(INPUT,1003)  (RPAR(L,K),L=1,KM)
C      RPAR(L,K): L-TER RADIALPARAMETER DER K-TEN ZERLEGUNG
C
       WRITE(NOUT,1005)   (RPAR(L,K),L=1,KM)
 1005  FORMAT(/22H SATZ RADIALPARAMETER   /3(8F12.6/))
C
       DO 240 N=1,KK
        READ(INPUT,1002) NUM(1,N,K),NUM(2,N,K),NUM(4,N,K)
        NUM(3,N,K)=I+N
C       NUM(3,.,.)  NUMMER DES BASISVEKTORS
C       NUM(1,.,.)    "    DER SPINFUNKTION DES BASISVEKTORS
C       NUM(2,.,.)    "     "  DREHIMPULSFUNKTION DES BASISVEKTORS
C       NUM(4,.,.)    "     "  POLYNOMFUNKTION DES BASISVEKTORS
C
        READ(INPUT,1003) (COF(L,N,K),L=1,LM)
C       COF(L,N,K)  KOEFFIZIENT, MIT DER DER L-TE SATZ INNERER WEITEN
C                   ZUM N-TEN BASISVEKTOR DER K-TEN ZERLEGUNG BETRAEGT
C
        WRITE(NOUT,1050) NUM(3,N,K),K,NUM(1,N,K),NUM(2,N,K),
     *                  NUM(4,N,K),(COF(L,N,K),L=1,LM)
 1050  FORMAT(//15H DEFINITION DER,I3,22H TEN BASISFUNKTION IST/
     1  I3,13H TE ZERLEGUNG,I3,23H TE SPINISOSPINFUNKTION/I3,
     2  26H TE BAHNDREHIMPULSFUNKTION,I3,' TES POLYNOM',/,
     3  ' SUMMATION UEBER INNERE WEITEN MIT KOEFF',/,(1P20E12.4))
  240  CONTINUE
C
       I=I+KK
   22 CONTINUE
C     ENDE DEFINITION BASISVEKTOREN
C
C
      IF(I.GT.NZBMAX) THEN
        WRITE(NOUT,1985)
        IF(I.GT.NZBMAX) call fehler (22)
      ELSE
        WRITE(NOUT,*)'ANZ BASISVEKTOREN = ',I
      ENDIF

 1985 FORMAT(' NZBMAX ZU KLEIN')
C
C
C     AUSSCHREIBEN DER ZU BERECHNENDEN OPERATOREN
C
2011  FORMAT(////18H BERECHNET WERDEN )
      WRITE (NOUT,2011)
      write(nout,1002)(LREG(K),K=1,NZOPER)
      IF(LREG(1).GT.0) WRITE(NOUT,1011)
      IF(LREG(2).GT.0) WRITE(NOUT,1012) MUL
      IF(LREG(3).GT.0) WRITE(NOUT,1013) MUL
      IF(LREG(4).GT.0) WRITE(NOUT,1014) MUL
      IF(LREG(5).GT.0) WRITE(NOUT,1015) MUL
      IF(LREG(6).GT.0) WRITE(NOUT,1016) MUL
      IF(LREG(7).GT.0) WRITE(NOUT,1017) MUL
      IF(LREG(8).GT.0) WRITE(NOUT,1018) MUL
      IF(LREG(9).GT.0) WRITE(NOUT,1019) MUL
      IF(LREG(10).GT.0) WRITE(NOUT,1020) MUL
      IF(LREG(11).GT.0) WRITE(NOUT,1021) MUL
      IF(LREG(12).GT.0) WRITE(NOUT,1022) MUL
      IF(LREG(13).GT.0) WRITE(NOUT,1023) MUL
C     ENDE AUSSCHREIBEN
C
      REWIND NBAND1
C      IF(NBAND5.LE.0) GOTO 800
CC
CC     DER FOLGENDE PROGRAMMTEIL WIRD NUR ZUM KOPIEREN VON NBAND5
CC     DURCHLAUFEN
C      READ(NBAND5) NZF1,MUL,(LREG1(K),K=1,NZOPER),I1,
C     * (NZRHO1(K),K=1,NZF1)
C      MFEH1 = 1
C      MFEH2 = NZF
C      MFEH3 = NZF1
C      IF(NZF.LT.NZF1) GOTO 890
C      DO 804 K=1,NZOPER
C      MFEH1 = 3+ K
C      MFEH2 = LREG(K)
C      MFEH3 = LREG1(K)
C      IF(LREG(K)-LREG1(K))890,804,806
C  806 DO 807 ML=1,NZF1
C      DO 807 MR=ML,NZF1
C  807 KOM(K,ML,MR)=0
C  804 CONTINUE
C      DO 805 K=1,NZF1
C      MFEH1=11 + K
C      MFEH2=NZRHO(K) 
C      MFEH3 = NZRHO1(K)
C      IF (MFEH2 .LT. MFEH3) GOTO 890
C  805 CONTINUE
C      DO 6020 ML=1,NZF
C      DO 6020 MR=1,ML
C      IF(ML.LE.NZF1) GOTO 6020
C      DO 6022 MC=1,NZOPER
C 6022 KOM(MC,ML,MR)=0
C 6020 CONTINUE
C      GO TO 800
C  890 WRITE(NOUT,891) MFEH1,MFEH2,MFEH3
C  891 FORMAT( 17H VERGLEICHSFEHLER ,3I4)
C      STOP
C
C
  800 CONTINUE
C
C     AUSSCHREIBEN DER BISHERIGEN WERTE AUS NBAND1
      WRITE(NBAND1) NZF,MUL,(LREG(K),K=1,NZOPER),I,(NZRHO(K),K=1,NZF)
      DO 950  K = 1,NZF
      M=MZG(K)
      N3=MZPAR(K)
      MC1=NZC(K)-1
      KK=NZRHO(K)
      IF(KK.LE.0) GOTO 950
      DO 4536 N=1,KK
      N1=NUM(1,N,K)
      N2=NUM(2,N,K)
       N4=NUM(4,N,K)
      WRITE(NBAND1) N3,MMASSE(1,N1,K),MMASSE(2,N1,K),MLAD(1,N1,K),
     1 MLAD(2,N1,K),MSS(1,N1,K),MSS(2,N1,K),MS(N1,K),
     2 (LZWERT(L,N2,K),L=1,5),(RPAR(L,K),L=1,N3),KP(MC1,N4,K)
 4536 CONTINUE
      IF(NBAND5.LE.0) GOTO 810
       IF(K.GT.NZF1)GOTO 810
       KK1=NZRHO1(K)
      DO 813 N=1,KK1
813    READ(NBAND5)
810   CONTINUE
  950 CONTINUE
C     ENDE AUSSCHREIBEN
         cquaf='DMOUT.'
         write(fnumber,*) naufset
         do i=1, 255
            if(fnumber(i:i).ne.' ') goto 370
         end do
 370      do j=i, 255
            if(fnumber(j:j).eq.' ') goto 380
         end do  
 380      do k=1, 80
            if(cquaf(k:k).eq.' ') goto 490
         end do
 490      oname = cquaf(1:k-1) // fnumber(i:j)
          open(unit=27,file=oname,status='replace',form='unformatted')

      nteil = 1
      narbeit = 0
      nsend= 1
 1933 format(a72)
      
      IF(INDEPO.NE.0) THEN                             
         MEFI=23                                           
         CLOSE(UNIT=8,STATUS='KEEP')                               
         NBAND2=MEFI                                         
      ENDIF   
         IF(INDELU.NE.0) THEN
         MEFLU=25
         CLOSE(UNIT=9,STATUS='KEEP')
         NBAND3=MEFLU
         ENDIF

      CLOSE(UNIT=10,STATUS='KEEP')
      
        icount=0
C
C     BESTIMMUNG DER ORTSMATRIXELEMENTE VOR EINSETZEN DER DREHIMPULSE
C
C
      DO 40  MFL = 1,NZF
C     LOOP ZERLEGUNGEN LINKS
C
      IZLW=NZLW(MFL)
      IZPW=NZPO(MFL)
      IRHO=NZRHO(MFL)
      IK1=MZPAR(MFL)
      MC = NZC(MFL)
      MC1= MC - 1
      MC2 = MC + MC1
      MC3 = MC2 - 1
      NC = NZPAR(MFL)
      NCC = MZPAR(MFL)
C
      MD2 = MD + MD1
      MD3 = MD2 - 1
      ND = NZPAR(MFR)
      NDD = MZPAR(MFR)
C 
C
      IF(INDEPO.EQ.0) THEN            
C        SVEC   JAKOBIMATRIX
C        RVEK   TRANPONIERTE JAKOBIMATRIX
         READ(NBAND2)   ((RVEC(M,N),M=1,NZV),N=1,NZT)
         READ(NBAND2)   ((SVEC(N,M),M=1,NZV),N=1,NZT)
         READ(NBAND2)   ((QFCL(M,K),M=1,NZV),K=1,MC2)   
      ENDIF
C
      DO 40 MFR=1,MFL
C     LOOP ZERLEGUNGEN RECHTS
C
      JZLW=NZLW(MFR)
       JZPW=NZPO(MFR)
      JRHO=NZRHO(MFR)
      JK1=MZPAR(MFR)
      WRITE (NOUT,3009) MFL, MFR
3009  FORMAT(////19H ZWISCHEN ZERLEGUNG,I3,15H UND ZERLEGUNG ,I3,
     124H WIRD ZUR ZEIT GERECHNET )
      MD=NZC(MFR)
      MD1= MD - 1
c
      IF(INDEPO.NE.0) THEN                            
         READ(INPUT,378) INFIL                      
 378     FORMAT(A80)               
         write(NOUT,*) 'Oeffnen? MFL,MFR,KOM2(MFL,MFR)',
     $        MFL,MFR,KOM2(MFL,MFR),INFIL
         if(KOM2(MFL,MFR).EQ.0) goto 9010
         write(NOUT,*) 'Oeffne ',INFIL
c        CLOSE(UNIT=6,STATUS='KEEP')
#ifdef  PC_FORTRAN
c        OPEN(UNIT=6,FILE='OUTPUT', ACCESS='APPEND')
#endif
#ifdef  CRAY_FORTRAN
c        OPEN(UNIT=6,FILE='OUTPUT', POSITION='APPEND')
#endif
         OPEN(UNIT=MEFI,FILE=INFIL,STATUS='OLD',FORM='UNFORMATTED')    
         REWIND MEFI                                      
         READ(NBAND2)   ((RVEC(M,N),M=1,NZV),N=1,NZT)                   
         READ(NBAND2)   ((SVEC(N,M),M=1,NZV),N=1,NZT)                
         READ(NBAND2)   ((QFCL(M,K),M=1,NZV),K=1,MC2)  
      ENDIF  
9010  IF(INDELU.NE.0) THEN
         READ(INPUT,378) INFILU
         if(KOM2(MFL,MFR).EQ.0) goto 9011
         OPEN(UNIT=MEFLU,FILE=INFILU,STATUS='OLD',FORM='UNFORMATTED')
          REWIND MEFLU
      ENDIF

      if(KOM2(MFL,MFR).EQ.0) goto 9011
      READ(NBAND2)    ((VW(M,K),M=1,NZT),K=1,MD1)
 9011 CONTINUE      
C
      DO 42 MKC =1, NZOPER
C     LOOP OPERATOREN
C
C      REWIND NBAND4
      NREB = 0
      IZ=(IZPW-1)*NDIM3+IZLW
      JZ=(JZPW-1)*NDIM3+JZLW
      NZAV=MC1+MD1
      IZLURE=1
      KWIED=0
      NZAVZU=0
C             1   2   3   4   5   6   7   8   9 SiP SiN      
      GOTO (401,402,403,404,405,406,407,408,409,410,411,412,413,414),MKC
401   IF (LREG(MKC).EQ.1) WRITE(NOUT,1011)
      GOTO 420
402   IF (LREG(MKC).EQ.1) WRITE(NOUT,1012) MUL
      GOTO 420
403   IF (LREG(MKC).EQ.1) WRITE(NOUT,1013) MUL
      GOTO 420
404   IF (LREG(MKC).EQ.1) WRITE (NOUT,1014) MUL
      IZLURE=2
      KWIED=1
      GOTO 420
405   IF (LREG(MKC).EQ.1) WRITE (NOUT,1015) MUL
      IZLURE=2
      KWIED=1
      GOTO 420
406   IF (LREG(MKC).EQ.1) WRITE (NOUT,1016) MUL
      NREB=1
      IZLURE=2*NZAV
      GOTO 420
407   IF (LREG(MKC).EQ.1) WRITE (NOUT,1017) MUL
      NREB=1
      IZLURE=2*NZAV
      GOTO 420
408   IF (LREG(MKC).EQ.1) WRITE (NOUT,1018) MUL
      NREB=1
      IZLURE=NZAV
      GOTO 420
409   IF (LREG(MKC).EQ.1) WRITE (NOUT,1019) MUL
      NREB=1
      IZLURE=NZAV
      GOTO 420
410   IF (LREG(MKC).EQ.1) WRITE (NOUT,1020) MUL
      GOTO 420
411   IF (LREG(MKC).EQ.1) WRITE (NOUT,1021) MUL
      GOTO 420
412   IF (LREG(MKC).EQ.1) WRITE (NOUT,1022) MUL
      IZLURE=NZAV+1
      NREB=1
      NZAVZU=1
      GOTO 420
413   IF (LREG(MKC).EQ.1) WRITE (NOUT,1023) MUL
      IZLURE=NZAV+1
      NREB=1
      NZAVZU=1
      GOTO 420
414   call fehler(414)
420   CONTINUE
C
C     EINLESEN DER SPIN-ISOSPINMATRIXELEMENTE AUS OBER
      NL = MZG(MFL)
      NR = MZG(MFR)
      IF(NREG(MKC).LE.0) GOTO 440
      READ(NBAND2) NTE,NTE,NTE
      NTEH=NTE
      IF (LREG(MKC).EQ.0) GOTO 425
      IF (KOM(MKC,MFL,MFR)) 421,422,423
  421 WRITE(NOUT,6034)
      GO TO  425
  422 WRITE(NOUT,6035)
      GO TO  425
  423 WRITE(NOUT,6036)
  425 CONTINUE
      IF(NTE.LE.0) GOTO 440
      DO 430   MTE = 1,NTE
C       READ(NBAND2)    NT1,(LC(K),K=1,NZT),
C     1   ((U(NFL,NFR),NFL=  1,NL),NFR=1,NR)        
c -------------------------------------------------------------------
C 1) READ        
       READ(NBAND2)   NT1ALT,(LCALT(K),K=1,NZT),
     1   ((UM(NFL,NFR),NFL=1,NL),NFR=1,NR)

C 2) UEBERPRUEFEN, OB SPIN-ISOSPIN-MATRIXELEMENT
C    GLEICH 0
       IMQ=0
       IM=-1
       DO 427, NFL=1, NL
       DO 427, NFR=1, NR
        IM = IM + 1
        NTEH=NTEH-1
        IF (ABS(UM(NFL,NFR)).LT.1.E-10) GOTO 427
        NTEH=NTEH+1
        IMQ = IMQ + 1
        IMV(IMQ) = IM
       
427    CONTINUE
429   NT1(MTE)=NT1ALT
      DO 27 I=1,NZT
27     LC(I,MTE)=LCALT(I)
       DO 28, NFL=1, NL
        DO 28, NFR=1, NR
         U(NFL,NFR,MTE)=UM(NFL,NFR)
28    CONTINUE
430   CONTINUE
      NTE=NTEH

C 3) reset indices
      DO 11 J=1,NR
      DO 11 I=1,NL
      U(I,J,MTE)= 0.
11    CONTINUE
      DO 12 IK=1,IMQ
      I= IMV(IK)/(NR*MKCM)
      IM=MOD(IMV(IK),NR*MKCM)
      J=IM/MKCM
      M=MOD(IM,MKCM)
12    U(I+1,J+1,M+1,MTE)=UM(IK)
      IF(IZT.EQ.NZT) THEN
      DO 7 I=1,IZT
7     LC(I,MTE)=LCALT(I)
      GOTO 921
      ENDIF
      DO 177 I=1,NZT
177    LC(I,MTE)=LC(I,MTE-1)
921   CONTINUE
c -------------------------------------------------------------------
       IF(MOBAUS.EQ.0) GOTO 426
       IF(LREG(MKC).EQ.0) GOTO 426
       WRITE (NOUT,838) NT1,(LC(K),K=1,NZT)
 838   FORMAT (/,1X,'VON OBEM',20I3)
       IF(MOBAUS.EQ.0) GOTO 426
       WRITE (NOUT,839) ((U(NFL,NFR),NFL=1,NL),NFR=1,NR)
 839   FORMAT (1P10E12.4)
C
426    CONTINUE

C     ENDE EINLESEN DER SPIN-ISOSPIN-MATRIXELEMENTE
C
440   KH=INT(MKC/2)+1
      IF (IENT(KH).LE.0) GOTO 42
C            1   2           5   6             SiP SiN
      GOTO(450,450,460,450,460,450,460,450,460,450,460,450,460,445),MKC
445   STOP 8
450   CALL LURE (MKC,ITV2)
C     EINLESEN DER ELEMENTE AUS LUISE
C
460   CONTINUE
      IF(KOM(MKC,MFL,MFR).EQ.0) GOTO 95
      if(mfl .ne. mfr) then
         nnn = nzrho(mfl)*nzrho(mfr)
      else
         nnn = (nzrho(mfl)*(nzrho(mfl)+1))/2
      endif
      if(lreg(mkc)*nnn .le. 0) goto 42
      
      IF(KOM(MKC,MFL,MFR).EQ.1) GOTO 942
      WRITE(NOUT,6037)
6037  FORMAT(' FUER DIESEN OPERATOR WIRD ERGAENZT')

95    IF (LREG(MKC).LE.0) GOTO 42
      IF(IRHO*JRHO.LE.0) GOTO 42
C
C     Nichts tun, koennen wir auch selber
      if( nte*itv2.eq.0 ) then


         write( nout,* ) 'Nichts zu tun.',nteil
         index(nteil)=0
         icount=icount+1

         WRITE (27) nteil,icount,index(nteil)
         WRITE (nout,*)'nteil,icount,index(nteil):',
     *                  nteil,icount,index(nteil)         


         goto 942
      endif

c schon alle beschaeftigt?
      nfertig=0
      write(NOUT,*) 'narbeit, nproc, nteil',narbeit,nproc,nteil
      if(narbeit .ge. (nproc-1) ) then
c        alle beschaeftigt, erstmal auf Ergebnisse warten
         call MPI_RECV(ngerechnet,1,MPI_INTEGER,MPI_ANY_SOURCE,
     *                 MPI_ANY_TAG,MPI_COMM_WORLD,istatus,ierr)
         call MPI_RECV(nteilalt,1,MPI_INTEGER,ngerechnet,0,
     *        MPI_COMM_WORLD,istatus,ierr)

         call MPI_RECV(iindy,lindy,MPI_INTEGER,ngerechnet,0,
     *        MPI_COMM_WORLD,istatus,ierr)
         call MPI_RECV(rbig,lbig,MPI_DOUBLE_PRECISION,ngerechnet,0,
     *        MPI_COMM_WORLD,istatus,ierr)


         narbeit = narbeit - 1
         write(NOUT,*)  'Rueckmeldung von ', ngerechnet
         index(nteilalt)=naufset
         icount=icount+1

         WRITE (27) nteilalt,icount,index(nteilalt)
         WRITE (nout,*)'nteilalt,icount,index(nteilalt):', 
     *                  nteilalt,icount,index(nteilalt)

         write(18) nteil,(index(i),i=1,nteil)

         write(27) (iindy(lb),lb=1,lindy)
         write(27) (rbig(lb),lb=1,lbig)

         nfertig=nfertig+1
C     Der Sklave darf gleich noch einmal, wenn noch Arbeit da ist
         nstarte=ngerechnet
      endif
            
c     noch einen starten
      if (nsend.lt. nproc) then
         if (nfertig.eq.0) then
            nstarte=nsend
         endif 
      endif


c     CLOSE(UNIT=6,STATUS='KEEP')
#ifdef  PC_FORTRAN
c        OPEN(UNIT=6,FILE='OUTPUT', ACCESS='APPEND')
#endif
#ifdef  CRAY_FORTRAN
c     OPEN(UNIT=6,FILE='OUTPUT', POSITION='APPEND')
#endif

c und mit Arbeit versorgen
      idummy=0
      call MPI_SEND(idummy,1,MPI_INTEGER,nstarte,0,
     *     MPI_COMM_WORLD,ierr)
      
      call MPI_SEND(nteil,1,MPI_INTEGER,nstarte,0,
     *     MPI_COMM_WORLD,ierr)


      read (coutput,'(80A1)') dummy

      do 9002 n=1,80
 9002    jdummy(n)=ichar(dummy(n))
      call MPI_SEND(jdummy(1),80,MPI_INTEGER,nstarte,0,
     *        MPI_COMM_WORLD,ierr)

      call MPI_SEND(rstore,lstore,MPI_DOUBLE_PRECISION,nstarte,0,
     *     MPI_COMM_WORLD,ierr)

      call MPI_SEND(ibeg,lbeg,MPI_INTEGER,nstarte,0,
     *     MPI_COMM_WORLD,ierr)

      call MPI_SEND(iluri,lluri,MPI_INTEGER,nstarte,0,
     *     MPI_COMM_WORLD,ierr)

      call MPI_SEND(rlurr,llurr,MPI_DOUBLE_PRECISION,nstarte,0,
     *     MPI_COMM_WORLD,ierr)

      call MPI_SEND(iimp,limp,MPI_INTEGER,nstarte,0,
     *     MPI_COMM_WORLD,ierr)

      call MPI_SEND(ikomi,lkomi,MPI_INTEGER,nstarte,0,
     *     MPI_COMM_WORLD,ierr)

      call MPI_SEND(rcsh,lcsh,MPI_DOUBLE_PRECISION,nstarte,0,
     *     MPI_COMM_WORLD,ierr)
      call MPI_SEND(rfltw,lfltw,MPI_DOUBLE_PRECISION,nstarte,0,
     *     MPI_COMM_WORLD,ierr)
      call MPI_SEND(ilup2,llup2,MPI_INTEGER,nstarte,0,
     *     MPI_COMM_WORLD,ierr)      
   
c     ... mit Arbeit versorgt
      
C     CALL DATE_AND_TIME (BIG_BEN (1), BIG_BEN (2), 
C    $        BIG_BEN (3), DATE_TIME)
C      write(NOUT,*) 'Rechner ',nstarte," um :: ",date_time,
C    $     ' abgefertigt'

      narbeit = narbeit + 1
      nsend= nsend +1
            
 942  nteil = nteil + 1

   42 CONTINUE
C        LOOP OPERATOREN
c     CLOSE(UNIT=6,STATUS='KEEP')
#ifdef  PC_FORTRAN
c        OPEN(UNIT=6,FILE='OUTPUT', ACCESS='APPEND')
#endif
#ifdef  CRAY_FORTRAN
c     OPEN(UNIT=6,FILE='OUTPUT', POSITION='APPEND')
#endif
   40 CONTINUE
      
c     hier muss die Schleife zum Aufsammeln der restlichen Ergebnisse hin
      do i=1, narbeit
         call MPI_RECV(ngerechnet,1,MPI_INTEGER,MPI_ANY_SOURCE,
     *        MPI_ANY_TAG,MPI_COMM_WORLD,istatus,ierr)
         call MPI_RECV(nteilalt,1,MPI_INTEGER,ngerechnet,0,
     *        MPI_COMM_WORLD,istatus,ierr)

      call MPI_RECV(iindy,lindy,MPI_INTEGER,ngerechnet,0,
     *        MPI_COMM_WORLD,istatus,ierr)
      call MPI_RECV(rbig,lbig,MPI_DOUBLE_PRECISION,ngerechnet,0,
     *     MPI_COMM_WORLD,istatus,ierr)

         index(nteilalt)=naufset
         icount=icount+1

         WRITE(27) nteilalt,icount,index(nteilalt)
         WRITE(nout,*) 'nteilalt,icount,index(nteilalt),rbig(1)',
     *          nteilalt,icount,index(nteilalt),(iindy(nb),nb=1,10)

         write(27) (iindy(lb),lb=1,lindy)
         write(27) (rbig(lb),lb=1,lbig)
         write(18) nteil,(index(iy),iy=1,nteil)

         istop=66
         call MPI_SEND(istop,1,MPI_INTEGER,istatus(MPI_SOURCE),
     *        0,MPI_COMM_WORLD,ierr)
         write(NOUT,*)  'Rueckmeldung VON ', ngerechnet
      enddo
c      
      istop=666
      
      write (NOUT,*) 'Aufsammeln'

      call MPI_BCAST(istop,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierror)

      write (NOUT,*) 'BROADCAST'
      write (NOUT,*) 'AUFSAMMELN MUSS NOCH MANUELL GEMACHT WERDEN!'
      rewind 18
         write(18) nteil,(index(i),i=1,nteil)
         write(NOUT,*) nteil,(index(i),i=1,nteil)
      close(unit=18,status='keep')


c     CLOSE(UNIT=6,STATUS='KEEP')
#ifdef  PC_FORTRAN
c        OPEN(UNIT=6,FILE='OUTPUT', ACCESS='APPEND')
#endif
#ifdef  CRAY_FORTRAN
c     OPEN(UNIT=6,FILE='OUTPUT', POSITION='APPEND')
#endif
      
      goto 666

C     Das Aufsammeln, wie bisher, wird nicht mehr durchlaufen, sondern
C     ist ein extra (serielles!) Programm.
C     Viel Spass beim Warten.


CC    RETURN
C     FUER SMIN AKTIVIEREN UND STOP DEAKTIVIEREN
C     Hier ist der Master fertig, laso das endif zum myid .eq. 0
      endif

C     ######################################################################
C     Jetzt kommt der quaf Sklave
C     ######################################################################

      write(slaveout,'(a7,i1)') "OUTPUT.",myid
      open(12, file=slaveout )
      write(slaveout,'(a7,i1)') "tempor.",myid
      open(4, file=slaveout ,FORM='UNFORMATTED')

      if (myid .ne. 0) then


 1234    continue

         call MPI_RECV(idummy,1,MPI_INTEGER,0,0,
     *        MPI_COMM_WORLD,istatus,ierr)
         if(idummy.eq.66) then
            call MPI_BCAST(istop,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierror)
            if(istop.eq.666) goto 666
         endif
      
         call MPI_RECV(nteil,1,MPI_INTEGER,0,0,
     *        MPI_COMM_WORLD,istatus,ierr)
      
         call MPI_RECV(jdummy(1),80,MPI_INTEGER,0,0,
     *           MPI_COMM_WORLD,istatus,ierr)
         do 9003 n=1,80
 9003       coutput(n:n)=char(jdummy(n))

C         CALL DATE_AND_TIME (BIG_BEN (1), BIG_BEN (2), 
C     $           BIG_BEN (3), DATE_TIME)
C         write(12,*) 'STARTE Empfang  um :: ',date_time
   
         call MPI_RECV(rstore,lstore,MPI_DOUBLE_PRECISION,0,0,
     *        MPI_COMM_WORLD,istatus,ierr)

C         CALL DATE_AND_TIME (BIG_BEN (1), BIG_BEN (2), 
C     $        BIG_BEN (3), DATE_TIME)
C         write(12,*) 'Empfangen lstore um :: ',date_time,lstore

         call MPI_RECV(ibeg,lbeg,MPI_INTEGER,0,0,
     *        MPI_COMM_WORLD,istatus,ierr)

C         CALL DATE_AND_TIME (BIG_BEN (1), BIG_BEN (2), 
C     $        BIG_BEN (3), DATE_TIME)
C         write(12,*) 'Empfangen lbeg um :: ',date_time,lbeg

         call MPI_RECV(iluri,lluri,MPI_INTEGER,0,0,
     *        MPI_COMM_WORLD,istatus,ierr)

         call MPI_RECV(rlurr,llurr,MPI_DOUBLE_PRECISION,0,0,
     *        MPI_COMM_WORLD,istatus,ierr)

         call MPI_RECV(iimp,limp,MPI_INTEGER,0,0,
     *        MPI_COMM_WORLD,istatus,ierr)

         call MPI_RECV(ikomi,lkomi,MPI_INTEGER,0,0,
     *        MPI_COMM_WORLD,istatus,ierr)
         call MPI_RECV(rcsh,lcsh,MPI_DOUBLE_PRECISION,0,0,
     *        MPI_COMM_WORLD,istatus,ierr)
         call MPI_RECV(rfltw,lfltw,MPI_DOUBLE_PRECISION,0,0,
     *        MPI_COMM_WORLD,istatus,ierr)
         call MPI_RECV(ilup2,llup2,MPI_INTEGER,0,0,
     *        MPI_COMM_WORLD,istatus,ierr)

c Namen fuer Slave bestimmen
C         CALL DATE_AND_TIME (BIG_BEN (1), BIG_BEN (2), 
C     $        BIG_BEN (3), DATE_TIME)
C         write(12,*) 'Empfangen Alles um :: ',date_time


         write(fnumber,*) nteil
         do i=1, 255
            if(fnumber(i:i).ne.' ') goto 70
         end do
 70      do j=i, 255
            if(fnumber(j:j).eq.' ') goto 80
         end do
 80      do k=1, 80
            if(coutput(k:k).eq.' ') goto 90
         end do
 90      oname = coutput(1:k-1) // fnumber(i:j)

      

C         TIME0=MCLOCK()/100.
C         RT0=timef()/1000.


c        write(15, *) 'Programm hat Bisher ',time0,'s gerechnet.'
C         CALL DATE_AND_TIME (BIG_BEN (1), BIG_BEN (2), 
C     $        BIG_BEN (3), DATE_TIME)
C         write(NOUT,*) date_time
C         write(NOUT,*) date_time,my_name
c          write(12,*)'NTE,ITV2',NTE,ITV2
C     Ein Aufruf
         CALL CALCU
C     ... und fertig
         
C         RT1=timef()/1000.

C         TIME1=MCLOCK()/100.
C         TIMED=MCLOCK()/100.-TIME0
c        write(15, *) 'Ich bin Nummer ',myid
c        
c        write(15, *) 'Programm hat bisher ',time0,'s gerechnet.'
c        write(15, *) 'Gerechnet: ',timed,' gebraucht: ',RT1,'s'
c        write(15, *) 'insgesamt jetzt ',time1,'s.'

C         CALL DATE_AND_TIME (BIG_BEN (1), BIG_BEN (2), 
C     $        BIG_BEN (3), DATE_TIME)
c        write(15,*) 'Beendet: ', date_time
         
C     Zurueckmelden
         
c        close(unit=15)
         call MPI_SEND(myid,1,MPI_INTEGER,0,0,MPI_COMM_WORLD,ierr)

         call MPI_SEND(nteil,1,MPI_INTEGER,0,0,
     *        MPI_COMM_WORLD,ierr)

      call MPI_SEND(iindy,lindy,MPI_INTEGER,0,0,
     *     MPI_COMM_WORLD,ierr)
      call MPI_SEND(rbig,lbig,MPI_DOUBLE_PRECISION,0,0,
     *     MPI_COMM_WORLD,ierr)


C     Sklave geht zurueck zum Anfang. Vielleicht kommt ja noch was nach
      
      goto 1234

      endif

C     Ende Hauptprogramm
      close(unit=6)
 666  continue
      close(unit=27,status='keep')
      ifehler = 42
      call MPI_FINALIZE(ierror)
      STOP
C       
6034  FORMAT(34H FUER DIESEN OPERATOR WIRD ERSETZT /)
6035  FORMAT(40H FUER DIESEN OPERATOR WIRD NEU GERECHNET /)
6036  FORMAT(34H FUER DIESEN OPERATOR WIRD KOPIERT /)
6100  FORMAT(/31H ES WIRD VOELLIG NEU GERECHNET /)
1002  FORMAT(20I3)
1011  FORMAT(/,1X,'NORMOPERATOR')
1012  FORMAT(/,1X,'E',I2,' SPINOPERATOR FUER PROTONEN')
1013  FORMAT(/,1X,'E',I2,' SPINOPERATOR FUER NEUTRONEN')
1014  FORMAT(/,1X,'M',I2,' SPINOPERATOR FUER PROTONEN')
1015  FORMAT(/,1X,'M',I2,' SPINOPERATOR FUER NEUTRONEN')
1016  FORMAT(/,1X,'E',I2,' BAHNOPERATOR FUER PROTONEN')
1017  FORMAT(/,1X,'E',I2,' BAHNOPERATOR FUER NEUTRONEN')
1018  FORMAT(/,1X,'M',I2,' BAHNOPERATOR FUER PROTONEN')
1019  FORMAT(/,1X,'M',I2,' BAHNOPERATOR FUER NEUTRONEN')
1020  FORMAT(/,1X,'E',I2,' COULOMB-/SIEGERT-OPERATOR FUER PROTONEN')
1021  FORMAT(/,1X,'E',I2,' COULOMB-/SIEGERT-OPERATOR FUER NEUTRONEN')
1022  FORMAT(/,1X,'E',I2,' KORREKTUROPERATOR FUER PROTONEN')
1023  FORMAT(/,1X,'E',I2,' KORREKTUROPERATOR FUER NEUTRONEN')

      END
C >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

      SUBROUTINE CALCU
      IMPLICIT double precision (A-H,O-Z)
C
      include 'par/QUAL'
C
      PARAMETER (NOUT=12, NBAND1=10, NBAND7=4)
C
      COMMON /STORE/ COF(NZPARM,NZRHOM,NZFMAX),RVEC(NZTMAX,NZTMAX),
     *    SVEC(NZTMAX,NZTMAX), U(MZGMAX,MZGMAX,NPDC),
     *    VW(NZTMAX,NZCMAX),
     *    QFCL(NZVMAX,2*NZCMAX-1),
     *    RPAR(MZPARM,NZFMAX),CPAR(NZAVMA,NZPARM,NZFMAX)
C
      COMMON /BEG/ IZLURE, NBAUS, JWIED, NZAVZU
C
      COMMON /LURI/ KVK(NDIM1), JQ(NDIM1), INDPO(NDIM1),
     *             MVK(NZIQMA,NDIM1),  NDIM3,
     *             IZ, JZ, KZHV(2,2*NZAVMA), LUPAUS, IQM
C
      COMMON /LURR/ EPO(NDIM1),WERTL(NDIM4,NDIM4)
C
      COMMON /IMP/ MKC,LAMB(NDIM4*NDIM4),LAHI(NZBMAX*NZBMAX),
     *       NUHI(NZBMAX*NZBMAX),
     *       MD1, MFL, NCC, NREB, NZAV, NZV1,
     *       LPARR, NZV,MC1
C
      COMMON /KOMIX/ NUM(4,NZRHOM,NZFMAX), NZC(NZFMAX),
     *          NZPAR(NZFMAX), MZPAR(NZFMAX),NZRHO(NZFMAX),
     *          LC(NZTMAX,NPDC),NGRV(2,NZCMAX,NZFMAX),
     *          NSH1(NDIM5), NSH2(NDIM5), NT1(NPDC),
     *          ITV2, NTE, NZT, MFR, NSH
C
      COMMON /LUR/ IENT(NZOPLU), NDIM1H, NBAND3
C
      COMMON /INDY/ IND(NZRHOM,NZRHOM)
      COMMON /BIG/ F(NDIM,NDIM),DM(NDIM,NDIM,NDIM2)
C
C  MULPOL common blocks -----------------------------------
      COMMON /LUP2/ NDIM4H,MUL,LINDEX(NDIM1),LINT(NZAVMA,NDIM1),
     *              LSUM(NDIM1), NDIM2H
C
      COMMON /CSH/ SH(NDIM5)
C
      COMMON /FLTW/ WERT(NDIM4,NDIM4,NDIM2)
C
      COMMON /CALBE/ QO(NZVMAX,NZVMAX), QOR(NZVMAX,NZVMAX),
     *        VR(NZTMAX,NZCMAX-1), WR(NZTMAX),
     *        VZ(NZVMAX,NZVMAX)
C
      DIMENSION QFCR(NZVMAX,NZVMAX,2*NZCMAX-1), LC1(NZTMAX),
     *          POR(NZVMAX,NZVMAX),
     *          POL(NZVMAX,NZVMAX)
C
      DIMENSION POLL(NZVMAX,NZVMAX), PORR(NZVMAX,NZVMAX)
C
      call MPI_COMM_RANK(MPI_COMM_WORLD,myid,ierr)
      write(nout,*)myid,' rechnet in CALCU!'

      PHI=3.1415926535897932384D0
      IRHO=NZRHO(MFL)
C      IK1=K1VEC(MFL)
      IK1=MZPAR(MFL)
      MC = NZC(MFL)
      MC1= MC - 1
      MC2 = MC + MC1
      MC3 = MC2 - 1
      NC= NZPAR(MFL)
      NCC = MZPAR(MFL)
      DO 490 M=1,NZT
  490 LC1(M)=0
      JRHO=NZRHO(MFR)
C      JK1=K1VEC(MFR)
      JK1=MZPAR(MFR)
      MD=NZC(MFR)
      MD1= MD - 1
      MD2 = MD + MD1
      MD3 = MD2 - 1
      ND = NZPAR(MFR)
      NDD = MZPAR(MFR)

464   IF (NTE*ITV2.LE.0) GOTO 900
C >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
      DO 44 MTE=1,NTE
C     LOOP DC'S
C
C      READ(NBAND4)    NT1,(LC(K),K=1,NZT),
C     1     ((U(NFL,NFR),NFL=1,NL),NFR=1,NR)
      MM=0
      DO 491 M=1,NZT
  491 MM=MM+IABS(LC(M)-LC1(M))
      IF(MM.LE.0) GOTO 492
      DO 494 M=1,NZT
  494 LC1(M)=LC(M)
      DO 50   M = 1,NZV
      DO 50   N = 1,M
      DO 51   K= 1,MD2
   51 QFCR(M,N,K) = .0
      DO 50 K=1,MD
   50 VR(M,K) = .0
      DO 290  K = 1,MD
      I1 = NGRV(1,K,MFR)
      I2 = NGRV(2,K,MFR) + I1 - 1
      IF(NGRV(2,K,MFR).LE.1)GOTO 290
      FZZ = 1./ FLOAT(NGRV(2,K,MFR))
      I3=I2-1
      DO 52  L1 = I1,I3
      I4=L1+1
      DO 52  L2 =I4,I2
      L3 = LC(L1)
      L4 = LC(L2)
      DO 52  M = 1,NZV
      DO 52  N = 1,NZV
      MM=MIN0(M,N)
      NN=MAX0(N,M)
   52 QFCR(NN,MM,K) = QFCR(NN,MM,K) + (RVEC(M,L3)-RVEC(M,L4))*
     1(RVEC(N,L3)-RVEC(N,L4)) *FZZ
  290 CONTINUE
      DO 53  K = 2,MD
      KI = MD1 + K
      DO 100   M = 1,NZT
      MM = LC(M)
      DO 100  N = 1,NZV
  100 VR(N,K-1) = VR(N,K-1) + VW(M,K-1) *RVEC(N,MM)
      DO 53   M = 1,NZV
      DO 53   N = 1,NZV
      MM=MIN0(M,N)
      NN=MAX0(N,M)
      QFCR(NN,MM,KI ) = QFCR(NN,MM,KI ) + VR(M,K-1)*VR(N,K-1)
   53 CONTINUE
  492 CONTINUE
      DO 54  M = 1,NZV
54    WR(M)=RVEC(M,NT1)
C      WRITE(NOUT,1100) WR(NZV),NZV
C 1100 FORMAT(//' J(A-1,I) = ',E16.8,I3)
C
C
  700 DO 60   KPAR = 1,NC
C     LOOP INNERE WEITEN LINKS
C
      DO 62   M = 1,NZV
      DO 62   N = 1,M
   62 POL (M,N) = .0
      DO 64   M = 1,NZV
      IF(MC3.LE.0) GOTO 64
      DO 63   K = 1,MC3
63    POL (M,M)=POL (M,M) + CPAR(K,KPAR,MFL)*QFCL(M,K)
   64 CONTINUE
C

C
      DO 70 LPAR=1,ND
C     LOOP INNERE WEITEN RECHTS
C
      JWIED=KWIED
C     KWIED/JWIED REGELT ZWEIFACHEN DURCHLAUF BEI MAGN. SPINOPERATOR
C     1. DURCHLAUF: JWIED=1  MUL+1
C     2. DURCHLAUF: JWIED=0  MUL-1
C
708   CONTINUE
C     WIEDERHOLPUNKT FUER MAGN. SPINOPERATOR
C
C
C
710   CONTINUE
      REWIND NBAND7
      DO 730, MM=1, IRHO
       DO 720, NN=1, JRHO
        IND(MM,NN)=0
720    CONTINUE
730   CONTINUE
C      IF (NBAND5.NE.0) READ (NBAND5) ((IND(MM,NN),NN=1,JRHO),
C     *                                            MM=1,IRHO)
      IZRS=((IRHO-1)/(NDIM/IK1))+1
      JZRS=((JRHO-1)/(NDIM/JK1))+1
      DO 880, KZRS=1,IZRS
C     LOOP BASISVEKTOREN*RADIALPARAMETER LINKS
C
      MKANA=(KZRS-1)*(NDIM/IK1)+1
      MKANE=MIN0(KZRS*(NDIM/IK1),IRHO)
      MKANZ=MKANE-MKANA+1
C
C
      DO 870, LZRS=1,JZRS
C     LOOP BASISVEKTOREN*RADIALPARAMETER RECHTS
C
      NKANA=(LZRS-1)*(NDIM/JK1)+1
      NKANE=MIN0(LZRS*(NDIM/JK1),JRHO)
      NKANZ=NKANE-NKANA+1
      IF (MFR.EQ.MFL.AND.MKANE.LT.NKANA) GOTO 870
C
C     FOLGENDER PROGRAMMTEIL DIENT NUR ZUM KOPIEREN BEI NBAND5#0
C      IF(KOM(MKC,MFL,MFR).EQ.0) GOTO 95
C      DO 540, M=1, MKANZ
C      MM=M-1+MKANA
CC
C      DO 530, N=1, NKANZ
C      NN=N-1+NKANA
CC
C      IF (IND(MM,NN).EQ.0) GOTO 530
C      READ (NBAND5) NUML,NUMR,IK1H,JK2H,LL1,
C     *              ((F(K,L),(JHI,DM(K,L,J),J=1,LL1),L=1,JK1),
C     *                                               K=1,IK1)
C      IF(KOM(MKC,MFL,MFR).LE.0) GOTO 530
C      WRITE (NBAND7) NUML,NUMR,IK1H,JK2H,LL1,
C     *               ((F(K,L),(DM(K,L,J),J=1,LL1),L=1,JK1),
C     *                                            K=1,IK1)
C530   CONTINUE
C  540 CONTINUE
C      IF(KOM(MKC,MFL,MFR).GT.0) GOTO 870
C95    CONTINUE
C
C
      DO 150 K1 = 1,NDIM
      DO 150 K2 = 1,NDIM
      DO 150 K3 = 1,NDIM2
      DM(K1,K2,K3) = .0
  150 CONTINUE
      I=0
C
C
      DO 101 M=1,MKANZ
C     LOOP BASISVEKTOREN LINKS
C
      MM = M-1+MKANA
      KSL = NUM(1,MM,MFL)
      KLL = NUM(2,MM,MFL)
      KPL=NUM(4,MM,MFL)
C      JPL=(KPL-1)*NDIM3+KLL
      JPL=KPL*NDIM3+KLL      
      NUML = NUM(3,MM,MFL)
      TS = COF(KPAR,MM,MFL)
      IF (ABS(TS).LT.1.E-10) GOTO 101
      MADL = (M-1) * IK1
C
C
      DO 105 N=1,NKANZ
C     LOOP BASISVEKTOREN RECHTS
C
      NN = N-1+NKANA
      KSR = NUM(1,NN,MFR)
      KLR = NUM(2,NN,MFR)
      KPR=NUM(4,NN,MFR)
C      JPR=(KPR-1)*NDIM3+KLR
      JPR=KPR*NDIM3+KLR      
      NUMR = NUM(3,NN,MFR)
      IF (NUML.LT.NUMR) GOTO 105
      IF(WERTL(JPL,JPR).NE.1.) GOTO 105
      A = TS * U(KSL,KSR)*COF(LPAR,NN,MFR)
      IF(ABS(A).LT.1.E-20) GOTO 105
      MADR = (N-1) * JK1
      IF (I.LT.NDIM5) GOTO 110
      WRITE(NOUT,111) I
 111  FORMAT (1X,'NDIM5 ZU KLEIN', I4)
      STOP
 110  I = I + 1
      NSH1(I) = NDIM * (MADR-1) + MADL
      NSH2(I) = NDIM4 * (JPR-1) + JPL
      SH(I) = A
      NUHI(I) = NZBMAX*(NUML-1) + NUMR
C
C     ENDE LOOP BASISVEKTOREN RECHTS
 105  CONTINUE
C
C     ENDE LOOP BASISVEKTOREN LINKS
 101  CONTINUE
C
      NSH = I
      IF (NSH.EQ.0) GOTO 80
      DO 72   M = 1,NZV
      DO 72   N = 1,NZV
   72 POR(M,N) = .0
      DO 74  M = 1,NZV
      DO 74  N = M,NZV
      IF(MD3.LE.0) GOTO 74
      DO 73  K = 1,MD3
   73 POR(N,M)=POR(N,M)+CPAR(K,LPAR,MFR)*QFCR(N,M,K)
   74 CONTINUE
      DO 67 M=1,NZV
      DO 67  N = M,NZV
67    POLL(N,M)=POL(N,M)
C
C
      DO 76 LPARR=1,NDD
C     LOOP RADIALWEITEN RECHTS
C
      DO 78 M=1,NZV
      DO 78 N=1,NZV
      QOR(N,M)=0.0
78    QO(N,M)=0.0
      DO 77 M=1,NZV
      DO 77 N=M,NZV
      PORR(N,M)=POR(N,M)+RPAR(LPARR,MFR)*QFCR(N,M,MD2)
      QOR(N,M)=PORR(N,M)
77    QO(N,M)=POLL(N,M)+PORR(N,M)
      CALL BEGRI
      
C
C     ENDE LOOP RADIALWEITEN RECHTS
   76 CONTINUE
C
C
   80 DO 480 M=1,MKANZ
C     LOOP BASISVEKTOREN LINKS
C
      MM=M-1+MKANA
      NUML=NUM(3,MM,MFL)
C
C
      DO 481 N=1,NKANZ
C     LOOP BASISVEKTOREN RECHTS
C
      NN=N-1+NKANA
      NUMR=NUM(3,NN,MFR)
      IF(NUML.LT.NUMR) GOTO 481
      NUHILF = NZBMAX*(NUML-1) + NUMR
      LL1 = LAHI(NUHILF)
      M1=(M-1)*IK1+1
      M2=M*IK1
      N1=(N-1)*JK1+1
      N2=N*JK1
      A = 0.
      DO 510 K=M1,M2
      DO 510 L=N1,N2
      DO 510 J=1,LL1
C      IF(NAUS.LT.7) GOTO 511
C      WRITE(nout,516) DM(K,L,J)
C  516 FORMAT(1X,' QUADI: DM = ',E16.8/)
  511 A = A + ABS(DM(K,L,J))
  510 CONTINUE
      IF (A.NE.0.) IND(MM,NN)=1
      IF(IND(MM,NN).EQ.0) GOTO 481
C
      WRITE(NOUT,112) MTE,KPAR,LPAR
  112 FORMAT(//I4,'TE DC, ',I3,'TE INNERE WEITE LINKS, ',
     $ I3,'TE INNERE WEITE RECHTS')
      WRITE(NOUT,1002)NUML,NUMR,IK1,JK1,LL1,IND(MM,NN)
C      
      WRITE (NBAND7) NUML, NUMR, IK1, JK1, LL1,
     *               ((F(K,L),(DM(K,L,J),J=1,LL1),L=N1,N2)
     *                           ,K=M1,M2)
C      IF(NAUS.LT.5) GOTO 481
c                    DO 520 K=M1,M2
c                    DO 520 L=N1,N2
c              520   WRITE(NOUT,1051) F(K,L),(J-1,DM(K,L,J),J=1,LL1)
c              1051  FORMAT(' F= ',E12.5,4(' IQ= ',I3,' DM= ',E12.5))
C
C     ENDE LOOP BASISVEKTOREN RECHTS
  481 CONTINUE
C
C     ENDE LOOP BASISVEKTOREN LINKS
  480 CONTINUE
C
C     ENDE LOOP BASISVEKTOREN*RADIALPARAMETER RECHTS
870   CONTINUE
C
C     ENDE LOOP BASISVEKTOREN*RADIALPARAMETER LINKS
880   CONTINUE
C
C     AUSSCHREIBEN DER ERGEBNISSE AUF NBAND1
C
C      WRITE (NBAND1) ((IND(MM,NN), NN=1, JRHO), MM=1, IRHO)
      REWIND NBAND7
C
      DO 920, MM=1, IRHO
C     LOOP BASISVEKTOREN LINKS
C
      DO 910, NN=1, JRHO
C     LOOP BASISVEKTOREN RECHTS

C
      IF (IND(MM,NN).EQ.0) GOTO 910
      READ (NBAND7) NUML, NUMR, IK1H, JK1H, LL1,
     *              ((F(K,L),(DM(K,L,J),J=1,LL1),L=1,JK1),
     *                                           K=1,IK1)
      WRITE(NOUT,1002)NUML,NUMR,IK1H,JK1H,LL1,IND(MM,NN)
      WRITE(NOUT,*)   NUML,NUMR,IK1H,JK1H,LL1,
     *          ((F(K,L),(J-1,DM(K,L,J),J=1,LL1),L=1,JK1H),
     *                                           K=1,IK1H)
C      WRITE (NBAND1) NUML, NUMR, IK1H, JK1H, LL1,
C     *               ((F(K,L),(J-1,DM(K,L,J),J=1,LL1),L=1,JK1),
C     *                                                K=1,IK1)
C
C     ENDE LOOP BASISVEKTOREN RECHTS
910   CONTINUE
C
C     ENDE LOOP BASISVEKTOREN LINKS
920   CONTINUE
C
C
      IF (JWIED.EQ.0) GOTO 70
      JWIED=0
      GOTO 708
C
C     ENDE LOOP INNERE WEITEN RECHTS
   70 CONTINUE
C
C     ENDE LOOP INNERE WEITEN LINKS
   60 CONTINUE
C
C     ENDE LOOP DC'S
   44 CONTINUE
C
  900 CONTINUE
C
1002  FORMAT(20I3)
      RETURN
      END

C >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>      
      SUBROUTINE LURE (MKC,KZAHL)
      IMPLICIT double precision (A-H,O-Z)
C     LURE LIEST DIE RECORDS AUS LUISE EIN
C
      INCLUDE 'par/QUAL'
C
C  commoners ----------------------------------------------
      COMMON /LURI/ KVK(NDIM1), JQ(NDIM1), INDPO(NDIM1),
     *             MVK(NZIQMA,NDIM1),  NDIM3,
     *             IZ, JZ, KZHV(2,2*NZAVMA), LUPAUS, IQM
C
      COMMON /LURR/ EPO(NDIM1),WERTL(NDIM4,NDIM4)
C
      COMMON /LUR/ IENT(NZOPLU), NDIM1H, NBAND3
C
      COMMON /BIG/ F(NDIM,NDIM),DM(NDIM,NDIM,NDIM2)
C
C  MULPOL common blocks -----------------------------------
      COMMON /LUP2/ NDIM4H,MUL,LINDEX(NDIM1),LINT(NZAVMA,NDIM1),
     *              LSUM(NDIM1), NDIM2H
C ---------------------------------------------------------
C
      COMMON /BEG/ IZLURE, NBAUS, JWIED, NZAVZU
C
      NOUT=6
      GOTO(20,20,100,20,100,20,100,20,100,20,100,20,100,10),MKC
10    STOP 40
20    DO 30 I=1,IZ
       DO 30 J=1,JZ
        WERTL(I,J)=0.
30    CONTINUE
C
C
      KZAHL=0
      IF (IZLURE.GT.2*NZAVMA) STOP 41
      DO 40 MIZ=1,IZLURE
      READ(NBAND3) KZAHL1,IQM
      IF(LUPAUS.GT.0) WRITE (NOUT,1000) MIZ,KZAHL1,IQM
      IQM=MAX0(IQM,1)
      IF (IQM.GT.NZIQMA) STOP 42
      KZHV(1,MIZ)=KZAHL
      KZHV(2,MIZ)=KZAHL1
      KZAHL2=KZAHL+1
      KZAHL=KZAHL+KZAHL1
      IF(KZAHL.GT.NDIM1) GOTO 110
      IF(KZAHL1.EQ.0) GOTO 40
      READ(NBAND3) (JQ(KW),INDPO(KW),EPO(KW),(MVK(IW,KW),IW=1,IQM),
     *             LINDEX(KW),KW=KZAHL2,KZAHL)
C
      CALL FAPOR(KZAHL2,KZAHL)
40    CONTINUE
C
C
100   RETURN
C
110   PRINT 1001,KZAHL1,KZAHL,NDIM1
      STOP 43
1000  FORMAT(1X,I2,'. RECORD ZUM OPERATOR AUS LUISE   ',
     *       'RECORDLAENGE:',I6,'  MAXIMALE ANZAHL DER ',
     *       'SIGMAFAKTOREN:',I4)
1001  FORMAT(1X,'GESAMTRECORD AUS LUISE ZU LANG',3I10)
      END
      SUBROUTINE FAPOR(IA,IE)
      IMPLICIT double precision (A-H,O-Z)
C     FAPOR REPRODUZIERT DIE DREHIMPULS- UND POLYNOMSTRUKTUR
C     DES MATRIXELEMENTS AUS INDPO UND DIE INNEREN DREHIMPULSE
C     AUS INDEX
C
      INCLUDE 'par/QUAL'
C
      COMMON /LURI/ KVK(NDIM1), JQ(NDIM1), INDPO(NDIM1),
     *             MVK(NZIQMA,NDIM1),  NDIM3,
     *             IZ, JZ, KZHV(2,2*NZAVMA), LUPAUS, IQM
C
      COMMON /LURR/ EPO(NDIM1),WERTL(NDIM4,NDIM4)
C
      COMMON /LUP2/ NDIM4H,MUL,LINDEX(NDIM1),LINT(NZAVMA,NDIM1),
     *              LSUM(NDIM1), NDIM2H
C 
      IF (NDIM4H.NE.NDIM4) STOP 50
C
C      IF(LUPAUS.EQ.0) GOTO 20
C      PRINT 1000,IA,IE
C      IF(LUPAUS.LT.2) GOTO 20
C      DO 10 I=IA,IE
C      PRINT 1001,JQ(I),INDPO(I),EPO(I),LINDEX(I),
C     1 (MVK(J,I),J=1,IQM)
C10    CONTINUE
C
20    DO 40, I=IA,IE
C
C      REPRODUKTION DER EXTERNEN DREHIMPULSE UND POLYNOME
C
       KPL=INDPO(I)/100000+1
       KPR=MOD(INDPO(I)/10000,10)+1
       LL=MOD(INDPO(I)/100,100)+1
       LR=MOD(INDPO(I),100)+1
       LPL=(KPL-1)*NDIM3+LL
       LPR=(KPR-1)*NDIM3+LR
       WERTL(LPL,LPR) = 1.
       L=NDIM4*(LPR-1)+LPL
       KVK(I)=L
C
C      REPRODUKTION DER INNEREN DREHIMPULSE
C
       INDEXH=LINDEX(I)
       LSUMH=0
       DO 30, KK=1, NZAVMA
        LINT(KK,I)=MOD(INDEXH,10)
        INDEXH=INT(INDEXH/10)
        LSUMH=LSUMH+LINT(KK,I)
30     CONTINUE
       IF (INT(LSUMH/2)+1.GT.NDIM2H) THEN
        PRINT 1002, I, INT(LSUMH/2)+1, NDIM2H
        call fehler(51)
       ENDIF
       LSUM(I)=LSUMH
40    CONTINUE
      RETURN
1000  FORMAT(1X,'VON LUISE ELEMENTE VON',I8,' BIS ',I10)
1001  FORMAT(1X,I5,' STRUKTUR',I10, ' EPO ',
     1 E16.8,'  INDEX ',I6,' MVK ',19I3)
1002  FORMAT(1X,'DIMENSIONIERUNGSFEHLER:  NDIM2 ZU KLEIN',/,1X,
     *       'LSUM(',I4,')/2*1:',I4,'   NDIM2:',I4)
      END
      SUBROUTINE BEGRI
      IMPLICIT double precision (A-H,O-Z)
C
C     BEGRI BESTIMMT DIE SIGMAFAKTOREN, DIE ALPHAS UND DIE VORFAKTOREN DER
C     ORTSRAUMATRIXELEMENTE, SUMMIERT UEBER ALLE RECORDS AUS LUISE EINES
C     OPERATORS AUS UND SPEICHERT DAS ERGEBNIS DES RED. ORTSRAUMMATRIXELEMENTS
C     IN DIE FELDER DMM UND DEN EXPONENTIALVORFAKTOR IN FG AB
C
      INCLUDE 'par/QUAL'
C
      PARAMETER (NOUT=12, NBAND1=10, NBAND7=4)
C
      COMMON /KOMIX/ NUM(4,NZRHOM,NZFMAX), NZC(NZFMAX),
     *          NZPAR(NZFMAX), MZPAR(NZFMAX),NZRHO(NZFMAX),
     *          LC(NZTMAX),NGRV(2,NZCMAX,NZFMAX),
     *          NSH1(NDIM5), NSH2(NDIM5), NT1(NPDC),
     *          ITV2, NTE, NZT, MFR, NSH
C
      COMMON /BEG/ IZLURE, NBAUS, JWIED, NZAVZU
C
      COMMON /BIG/ F(NDIM,NDIM),DM(NDIM,NDIM,NDIM2)
C
      COMMON /CSH/ SH(NDIM5)
C
      COMMON /CALBE/ QO(NZVMAX,NZVMAX), QOR(NZVMAX,NZVMAX),
     *               VR(NZTMAX,NZCMAX-1), WR(NZTMAX),
     *               VZ(NZVMAX,NZVMAX)
C
      COMMON /LURI/ KVK(NDIM1), JQ(NDIM1), INDPO(NDIM1),
     *             MVK(NZIQMA,NDIM1),  NDIM3,
     *             IZ, JZ, KZHV(2,2*NZAVMA), LUPAUS, IQM
C
      COMMON /FLTW/ WERT(NDIM4,NDIM4,NDIM2)
C
      COMMON /STORE/ COF(NZPARM,NZRHOM,NZFMAX),RVEC(NZTMAX,NZTMAX),
     *               SVEC(NZTMAX,NZTMAX), U(MZGMAX,MZGMAX,NPDC),
     *               VW(NZTMAX,NZCMAX),
     *               QFCL(NZVMAX,2*NZCMAX-1),
     *               RPAR(MZPARM,NZFMAX),CPAR(NZAVMA,NZPARM,NZFMAX)
C     
      COMMON /IMP/ MKC,LAMB(NDIM4*NDIM4),LAHI(NZBMAX*NZBMAX),
     *       NUHI(NZBMAX*NZBMAX),
     *       MD1, MFL, NCC, NREB, NZAV, NZV1,
     *       LPARR, NZV,MC1
C
      DIMENSION P(NZAVMA,NZVMAX),BETA(NZVMAX),S(NZCMAX*(NZAVMA+1))
      DIMENSION D(NZVMAX),OMEGA(NZVMAX)
      DIMENSION ALPHA(NZAVMA),ALPHAS(NZAVMA)
      DIMENSION ZH(NZTMAX), ZW(NZTMAX), SS(NZCMAX*(NZAVMA+1))
      DIMENSION WERTT(NDIM4*NDIM4,NDIM2), DMM(NDIM*NDIM,NDIM2),
     *          FGG(NDIM*NDIM)
      DIMENSION ZY(NZVMAX,NZAVMA+1), RHO(NZVMAX), Y(NZAVMA+1),
     *          GAM(NZAVMA), Z(NZVMAX,NZAVMA+1), H(2*(NZAVMA+1))
      EQUIVALENCE (WERT(1,1,1),WERTT(1,1)), (DM(1,1,1),DMM(1,1))
      EQUIVALENCE (F(1,1),FGG(1))
C
C     AUSDRUCK UEBERGEBENER WERTE
C
      PHI=3.1415926536
      PHI2 =  SQRT(PHI)
C     INITIALISIERUNG DES UNTERPROGRAMMS
C
      H(1)=0.
      H(2)=0.
      DO 10, N=1,NZAV
       NH1=2*(N+NZAVZU)
       NH2=NH1-1
       ALPHA(N)=0.
       ALPHAS(N)=0.
       H(NH1)=0.
       H(NH2)=0.
10    CONTINUE
C
      DO 20   M = 1,NZV
       OMEGA(M)=0.
       RHO(M)=0.
       ZY(M,1)=0.
       DO 20   N = 1,NZAV
        NH=N+NZAVZU
        P(N,M) = 0.
        ZY(M,NH)=0.
20    CONTINUE
C
      INZ1=(NZAV*(NZAV+1))/2
      DO 30 M=1,INZ1
       S(M)=0.
30    CONTINUE
C     ENDE INITIALISIERUNG
C
C
C         DO 44 N = 1,NZV
C          write(nout,1002)(QO(N,M),M=1,NZV)
C   1002  FORMAT(1X,'QO(N,M)',(/,1X,8E12.4))
C      44 CONTINUE
C         write(nout,1005) (WR(M),M=1,NZV)
C   1005  FORMAT(1X,'WR(M)  :',8E12.4)
C         DO 46, M=1, NZV
C          write(nout,1006) M, (VZ(M,N), N=1, NZV)
C   46    CONTINUE
C   1006  FORMAT(1X,'VZ(',I2,',N)  :',8E12.4)
C         DO 47, N=1, NZV
C          write(nout,1007) N, (VR(N,M), M=1, MD1)
C   47    CONTINUE
C   1007  FORMAT(1X,'VR(',I2,',M)  :',8E12.4)
C
   82 CONTINUE
C
C
C     FUER DIFFERENTIALOPERATOREN: BESTIMMUNG VON GAM UND RHO
      IF(NREB.EQ.0) GOTO 90
      DO 84   N = 1,NZAV
       GAM(N) = .0
       IF(N-MC1)   84,84,85
   85  M = N - MC1
       DO 86   K = 1,NZV
        GAM(N) = GAM(N) + VR(K,M)*WR(K)
86     CONTINUE
   84 CONTINUE
C
      DO 87   M = 1,NZV
       DO 87   K = 1,NZV
        RHO(M) = RHO(M) + WR(K)*(QOR(K,M)+QOR(M,K))
87    CONTINUE
   90 CONTINUE
C
C
C     HAUPTACHSENTRANSFORMATION IN SUBR. HAUPT; EIGENWERTE WERDEN IN QO
C     GESPEICHERT
      CALL HAUPT(MKC,NZV)
C
C
C     AUSDRUCK DER EIGENWERTE
C           DO 25 M = 1,NZV
C            write(12,1004)M,QO(M,M)
C     1004  FORMAT(' QO(',I2,') = ',E12.4)
C        25 CONTINUE
C
C     IM FOLEGENDEN WERDEN DIE GROESSEN S(I), ALPHA(I) UND BETA(I) SO
C     WEIT BESTIMMT, WIE ES OHNE KENNTNIS DER RADIALWEITEN LI MOEGLICH
C     IST
      FFF=1.
      IF(NZV1.EQ.0) GOTO 71
      DO 1   M = 1,NZV1
       BETA(M)=1./ SQRT(QO(M,M))
       FFF=PHI2*BETA(M)*FFF
1     CONTINUE
   71 CONTINUE
      BETA(NZV) = 1.
C
      DO 2 N=1,MC1
       NN=NZV-MC1+N
       DO 2 M=NN,NZV
        P(N,M)=VZ(M,NN)*BETA(M)
2     CONTINUE
      DO 4 N=1,MD1
       NN=MC1+N
       DO 4 M=1,NZV
        DO 6 K=1,M
         P(NN,M)=VR(K,N)*VZ(M,K)+P(NN,M)
6       CONTINUE
        P(NN,M)=P(NN,M)*BETA(M)
4     CONTINUE
C
      DO 8 M=1,NZV
       DO 9 K=1,M
        OMEGA(M)=OMEGA(M)+WR(K)*VZ(M,K)*BETA(M)
    9  CONTINUE
    8 CONTINUE
      OMH=OMEGA(NZV)
C
      DO 3  M = 1,NZAV
    3 ZH(M) = P(M,NZV)
      DO 5   M = 1,NZV
    5 ZW(M) = VZ(NZV,M)
      IF(NZV1.EQ.0) GOTO 14
      DO 11 M=1,NZV1
       DO 11 N=1,M
        VZ(M,N)=VZ(M,N)*BETA(M)
11    CONTINUE
C
      I=0
      DO 12 M=1,NZAV
       DO 12 N=M,NZAV
        I=I+1
        DO 12   K = 1,NZV1
         S(I)=S(I)+P(N,K)*P(M,K)
12    CONTINUE
C
      DO 212 M=1,NZAV
       DO 212 N=1,NZV1
        ALPHA(M)=ALPHA(M)+OMEGA(N)*P(M,N)
212   CONTINUE
C
   14 CONTINUE
      IF(NREB.EQ.0) GOTO 45
      IF (NZV1.EQ.0) GOTO 45
      DO 62   M = 1,NZV1
       DO 62   K = M,NZV1
        IF (NZAVZU.EQ.1) ZY(M,1)=ZY(M,1)+VZ(K,M)*OMEGA(K)
        DO 62    N = 1,NZAV
         NH=N+NZAVZU
         ZY(M,NH) = ZY(M,NH) + VZ(K,M)*P(N,K)
62    CONTINUE
  45  CONTINUE
C
C                 write(nout,1008) (OMEGA(M), M=1, NZV)
C           1008  FORMAT (1X,9HOMEGA'S :,8E12.4)
C                 write(nout,1009) (ALPHA(M), M=1, NZAV)
C           1009  FORMAT (1X,9HALPHA'S :,8E12.4)
C                 DO 48, N=1, NZV
C                  write(nout,1010) N, (P(M,N), M=1, NZAV)
C           48    CONTINUE
C           1010  FORMAT(1X,'P(M,',I2,') :',8E12.4)
C                 write(nout,1011) (S(I), I=1, (NZAV*(NZAV+1)/2))
C           1011  FORMAT(1X,'S(I) :',8E12.4)
C
C
C
55    DO 170  KPARR = 1,NCC
C     LOOP RADIALWEITEN LINKS
C
      FG = .0
      DO 150 M = 1,NZV
  150 D(M) = .0
      BETA(NZV) = 1./ SQRT(QO(NZV,NZV)+RPAR(KPARR,MFL))
      OMEGA(NZV)=OMH*BETA(NZV)
      FF = FFF *PHI2*BETA(NZV)
      FF=FF**3
C
      DO 91   M = 1,NZAV
       P(M,NZV) = ZH(M)*BETA(NZV)
91    CONTINUE
C
      DO 92    M = 1,NZV
       FG=FG+OMEGA(M)*OMEGA(M)
       VZ(NZV,M) = ZW(M)*BETA(NZV)
92    CONTINUE
      FG=FG/(-4)
C
      I = 0
      DO 93   M = 1,NZAV
       ALPHAS(M)=ALPHA(M)+OMEGA(NZV)*P(M,NZV)
       DO 93   N = M,NZAV
        I = I + 1
        SS(I) = 2.*(S(I) + P(N,NZV)*P(M,NZV))
93    CONTINUE
C
C     BESTIMMUNG DER VORFAKTOREN DES MATRIXELEMENTS
C
      H(1)=FF
      H(2)=FF
      IF (NREB.LE.0) GOTO 210
C
C     NUR DIFFERENTIALOPERATOREN
      DO 21 M=1,NZV
       IF (NZAVZU.EQ.1) Z(M,1)=ZY(M,1)+VZ(NZV,M)*OMEGA(NZV)
       DO 21 N=1,NZAV
        NH=N+NZAVZU
        Z(M,NH)=ZY(M,NH)+VZ(NZV,M)*P(N,NZV)
21    CONTINUE
C
      DO 24 N=1, (NZAV+NZAVZU)
       Y(N) = .0
       DO 24 K=1,NZV
        Y(N) = Y(N) + RHO(K)*Z(K,N)
24    CONTINUE
C
      IF (NZAVZU.NE.1) GOTO 190
      H(1)=FF*(Y(1)-1.)/2.
      H(2)=H(1)
C     AN DIESER STELLE WIRD FUER DEN KORREKTUROPERATOR DAS FALSCHE
C     VORZEICHEN DES 1. RECORDS KORRIGIERT
C
190   CONTINUE
      DO 110 N=1,NZAV
       NH=N+NZAVZU
       NH1=2*NH
       NH2=NH1-1
       H(NH1)=FF*(GAM(N)-Y(NH)*.5)
       H(NH2)=H(NH1)
110   CONTINUE
C
210   CONTINUE
C     ALLE OPERATOREN
C
C     BERECHNUNG DER SIGMAFAKTOREN UND AUFSUMMATION DER ZUSAMMENGEHOERIGEN
C     ORTSRAUMMATRIXELEMENTE IN MAT
      NQ=IZLURE
      CALL MAT(H,SS,NQ,ALPHAS)
C
C     ADREESRECHNUNG UND BESTIMMUNG VON DMM
      I0 = NDIM * LPARR + KPARR
      DO 100 M=1,NSH
       I2 = NSH2(M)
       LL1 = LAMB(I2)
       MHI = NUHI(M)
       LAHI(MHI) = LL1
       DO 100 K = 1,LL1
        C = SH(M) * WERTT(I2,K)
c        write(nout,53)MKC,SH(M) , WERTT(I2,K)
        IF (C.EQ.0.) GOTO 100
        I1 = NSH1(M) + I0
        FGG(I1) = FG
        DMM(I1,K) = SH(M) * WERTT(I2,K)
      
c        write(nout,50) M,NSH,K,LL1,SH(M),WERTT(I2,K),DMM(I1,K),FG
c        PRINT 50,M,NSH,K,LL1,SH(M),WERTT(I2,K),DMM(I1,K),FG
   50 FORMAT(' BEGRI: M, NSH, K, LL1, SH(M), WERTT, DMM, FG:',
     $ 4I4,4E16.8/)
   53 FORMAT(' BEGRI: MKC, SH(M), WERTT:',
     $ 1I4,2E16.8/)   
  100 CONTINUE
C
C     ENDE LOOP RADIALPARAMETER LI
  170 CONTINUE
      RETURN
      END
      SUBROUTINE HAUPT(MKC, NZV)
      IMPLICIT double precision (A-H,O-Z)
C    HAUPT FUEHRT HAUPTACHSENTRANSFORMATION DURCH
C    TRANSFORMATION IN VZ, HAUPTACHSENWERTE IN QO
C
      INCLUDE 'par/QUAL'
C
      COMMON /CALBE/ QO(NZVMAX,NZVMAX), QOR(NZVMAX,NZVMAX),
     *         VR(NZTMAX,NZCMAX-1), WR(NZTMAX),
     *        VZ(NZVMAX,NZVMAX)
C
C
      DO 1   K = 1,NZV
C
      VZ(K,K) = 1.
      IF(NZV-K)   1,1,2
    2 A = .5/QO(K,K)
      K1 = K + 1
      DO 7   M = K1,NZV
      VZ(M,K) = .0
      B = A*QO(M,K)
      DO 5   L = 1,K
    5 VZ(M,L) = VZ(M,L) - B*VZ(K,L)
      QO(M,M) = QO(M,M) - B*B*QO(K,K)
      IF(NZV-M)   7,7,6
    6 M1 = M + 1
      DO 4   N = M1,NZV
    4 QO(N,M) = QO(N,M) - A*QO(M,K)*QO(N,K)
    7 CONTINUE
    1 CONTINUE
      RETURN
       END
      SUBROUTINE MAT(H,SIGMA,NQ,ALPHA)
      IMPLICIT double precision (A-H,O-Z)
C     MAT BERECHNET MIT DEN SIGMA'S UND ALPHA'S AUS BEGRI
C     UND DEN MVK'S UND LI'S AUS LUISE DEN WERT DES MATR.-EL.
C
      INCLUDE 'par/QUAL'
C
      COMMON /LURI/ KVK(NDIM1), JQ(NDIM1), INDPO(NDIM1),
     *             MVK(NZIQMA,NDIM1),  NDIM3,
     *             IZ, JZ, KZHV(2,2*NZAVMA), LUPAUS, IQM
C
      COMMON /LURR/ EPO(NDIM1),WERTL(NDIM4,NDIM4)
C
      COMMON /LUP2/ NDIM4H,MUL,LINDEX(NDIM1),LINT(NZAVMA,NDIM1),
     *              LSUM(NDIM1), NDIM2H
C
      COMMON /BEG/ IZLURE,NBAUS,JWIED,NZAVZU
C
      COMMON /FLTW/ WERT(NDIM4,NDIM4,NDIM2)
      COMMON /IMP/ MKC,LAMB(NDIM4*NDIM4),LAHI(NZBMAX*NZBMAX),
     *       NUHI(NZBMAX*NZBMAX),
     *       MD1, MFL, NCC, NREB, NZAV, NZV1,
     *       LPARR, NZV,MC1
C
      DIMENSION SIGMA(NZCMAX*(NZAVMA+1)), ALPHA(NZAVMA),
     *          H(2*(NZAVMA+1)), WERTT(NDIM4*NDIM4,NDIM2)
      EQUIVALENCE (WERT(1,1,1),WERTT(1,1))
C
      IF (NDIM4H.NE.NDIM4) STOP 90
C
C      IF(LUPAUS.LT.5) GOTO 10
C      PRINT 1000,(ALPHA(I),I=1,NZAV)
C      PRINT 1001, (SIGMA(I), I=1,(NZAV*(NZAV+1)/2))
C      PRINT 1003, (H(I),I=1,(2*(NZAV+NZAVZU)))
C
10    CONTINUE
C     ANFANGSWERTZUWEISUNG
      DO 20   LR = 1,JZ
       DO 20   LL = 1,IZ
        I = (LR - 1)*NDIM4 + LL
        LAMB(I) = 0
        DO 20   K = 1, NDIM2
         WERT(LL,LR,K)=0.
20    CONTINUE
C
C     EINLESEN DER RECORDLAENGEN DES OPERATORS
C
      NL=1
      NU=NQ
      GOTO (24,24,24,22,22,24,24,24,24,24,24,24,24,21), MKC
21    STOP 91
22    CONTINUE
C     MAGN. SPINOPERATOR
      NL=NQ-JWIED
      NU=NQ-JWIED
24    DO 80, N=NL,NU
       K3ZAHL = KZHV(2,N) 
       IF(K3ZAHL.EQ.0) GOTO 80
       K1ZAHL = KZHV(1,N) + 1
       K2ZAHL = K1ZAHL + K3ZAHL - 1
       VORFAK=H(N)
       IF (MKC.EQ.8.OR.MKC.EQ.9) VORFAK=H(2*N)
       IF (MKC.EQ.12.OR.MKC.EQ.13) VORFAK=H(2*N)
C
C      BERECHNUNG DER MATRIXELELMENTE
       DO 70 K=K1ZAHL,K2ZAHL
        MM = JQ(K)
        KZOM=LSUM(K)/2+1
        IF (NZAVZU.EQ.1.AND.N.EQ.1) KZOM=(LSUM(K)+1)/2+1
c        write(12,*)'MAT:   LASUM,KZOM  ',LSUM(K),KZOM
C       BEIM 1. RECORD DES KORREKTUROPERATORS WIRD AN DIESER STELLE
C       MIT DEM WELLENVEKTOR K MULTIPLIZIERT
        SIG=VORFAK
        IF(MM.EQ.0) GOTO 40
        DO 30, M = 1,MM
         MN = MVK(M,K)
         SIG = SIG*SIGMA(MN)
30      CONTINUE
40      IF (LSUM(K) .EQ. 0) GOTO 60
        DO 50 M=1,NZAV
         IF (LINT(M,K) .EQ. 0) GOTO 50
         SIG=SIG*ALPHA(M)**LINT(M,K)
50      CONTINUE
60      M1 = KVK(K)
        LAMB(M1) = MAX(KZOM,LAMB(M1))
        WERTT(M1,KZOM) = WERTT(M1,KZOM) + EPO(K) * SIG
C        IF(LUPAUS.LT.3) GOTO 70
C        write(12,1002)K1ZAHL,K2ZAHL,K,MM,M1,KZOM,LAMB(M1),VORFAK,
C     $          SIG,EPO(K),WERTT(M1,KZOM),(LINT(M,K),M=1,NZAV)
70     CONTINUE
C
80    CONTINUE
C
      RETURN
1000  FORMAT(1X,8HALPHA'S:,8E13.5)
1001  FORMAT(1X,8HSIGMA'S:,8E13.5)
1002  FORMAT(1X,'K1,K2,K,MM,M1,KZOM,LAMB(M1),FF,SIG,EPO,WERTT,LINT:'
     *       ,/,2I3,2I4,I5,I3,I5,4E16.8,6I3)
1003  FORMAT(1X,'H-WERTE:',8E13.5,/,9X,8E13.5)
      END

      subroutine fehler(Nummer)
      IMPLICIT double precision (A-H,O-Z)
C     Statt einen Prozess mit STOP abzuschiessen, sollen alle
C     Prozesse von dem Unheil unterrichtet werden und daran 
C     eingehen. Scheint nicht richtig zu funktionieren.
      include "mpif.h"
      PARAMETER (NOUT=6)
      INTEGER iStatus( MPI_STATUS_SIZE )      
      call MPI_COMM_RANK(MPI_COMM_WORLD,myid,ierr)
      call MPI_COMM_SIZE(MPI_COMM_WORLD,nproc,ierr)

      do i=1,nproc-1 
         istop=66
         call MPI_SEND(istop,1,MPI_INTEGER,i,0,MPI_COMM_WORLD,ierr)
      enddo

      
      write (NOUT,*) 'FEHLER! Gestoppt.'
      write (NOUT,*) 'STOP ',nummer
      if (nummer.eq.261) write (NOUT,*) 'NZC?'
      if (nummer.eq.277) write (NOUT,*) 'NOL?'
      if (nummer.eq.9000) write (NOUT,*) 'MZPARM ZU KLEIN?'
      if (nummer.eq.9000) write (NOUT,*) 'NPDC?'

      istop=666
      
      call MPI_BCAST(istop,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierror)

      CLOSE(UNIT=6,STATUS='KEEP')

      call MPI_ABORT(MPI_COMM_WORLD,ifehler,ierror)
      call MPI_FINALIZE(ierr)
      STOP 666
      end