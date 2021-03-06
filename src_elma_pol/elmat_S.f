      PROGRAM ELMAT_S
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
C   MATRIXELEMENT IN WELLER DEFINITION ZWISCHEN J UND J'
C   FUER EINE MULTIPOLARITAET L
C
       INCLUDE 'par/elmat'
C
C  27.5.86   H.M.H. AUS ASFAK,MIT NEUEN FAKTOREN
C         GEAENDERT 24.6.86 H.M.H.
C         INPUT UND UEBERGABE VON VERSEM GEAENDERT 18.7.86 H.M.H.
C         AENDERUNG 11.03.87 T.S.
C        (VORFAKTOREN UND FEHLER IN NORTRA KORRIGIERT)
C         AENDERUNG 05.05.87 B.W.
C        (PHASEN AN WELLER ANGEPASST: (-)**(JWSL+JWSR+MUL) * I**MUL FUER
C         DAZU -I FUER DIE MAGNETISCHEN MATRIXELEMENTE )
C         AENDERUNG 19.05.88 M.U.
C        (DISTORSIONSKANAELE EINGEBAUT)
C      30.1.03 AUF PARAMETER UMGESTELLT H.M.H.
C      7.2.3 UEBERGABE VON S-POLE, EIN- UND AUSLAUFENDE WELLEN H.M.H.
C     17.2.3 FUER DISTORTIONSKANAELE WIRD SNORM AUS WNORM S-POLE UEBER
C            NOMMEN H.M.H
C
C
      COMPLEX *16 OPE,OPM,OPES,CPHASE, A, AA, CWF, OPWERT
      DIMENSION EBIN(NZKMAX),LWERT(NZKMAX),OPW(NDIM,NZOPER)
      DIMENSION ISTRN(NZKMAX) ,EW(NZAENM), EK(NZKMAX), REDM(NZKMAX),
     *  JWERT(3,NZKMAX),MASSE(2,NZKMAX) ,MLAD(2,NZKMAX), HQOM(NZKMAX),
     *  EKG(NZKMAX),EMOM(NZKMAX),NENTK(NZKMAX)
      DIMENSION IZQ(NZKMAX+1),IZP(NZKMAX),DF(13)
      DIMENSION PAR(NZPARM,NZKMAX),NAR(NZPARM,NZKMAX)
      DIMENSION SNORM(NZKMAX), WNORM(NZKMAX)
C     DIMENSION MREG(NZOPER),CWF(NZKBMA,NZKMAX,NZPARM)
      CHARACTER*10 METH(2)
      CHARACTER*10 IFRAK(7)
C
      COMMON /MART/ A(NZKBMA,NZKBMA,2), OPWERT(NZKBMA,NZOPER+2),
     *  MREG(NZOPER),CWF(NZKBMA,NZKMAX,NZPARM)
C
      DATA IFRAK /'N ','H ','HE','LI','BE','B ','C '/
      DATA METH/'S-DIREKT  ','S-INVERS  '/
C
      OPEN(UNIT=5,FILE='INEX',STATUS='OLD')
      OPEN(UNIT=6,FILE='OUTPUT')
      OPEN(UNIT=11,FILE='ENOUT',STATUS='OLD',FORM='UNFORMATTED')
      OPEN(UNIT=12,FILE='ELMA_SOUT',STATUS='OLD',FORM='UNFORMATTED')
      OPEN(UNIT=15,FILE='OUTEM_S',STATUS='UNKNOWN',FORM='FORMATTED')
C
      RNULL=0.
      REINS= 1.
      INPUT=5
      H2M=197.32858**2/(938.2796+939.5731)
      READ(INPUT,1000) NBAND1,NBAND2,NPRI,NPRIX,NPRICN,NZKAPH
      READ(INPUT,1000)(MREG(MKC),MKC=1,NZOPER)
C    NBAND1 VON ENELMAS, NBAND2 VON VERSEM
C    NPRI=1: RED.MATRIXELEMENTE <KANAL//OP//BOUND> WERDEN GEDRUCKT
C    NPRI=4: UEBERLAGERUNGSKOEFFIZIENTEN CWF VON VERSEM WERDEN
C            ZUSAETZLICH GEDRUCKT
C    NPRI=3: A-MATRIX VON VERSEM WIRD ZUSAETZLICH GEDRUCKT.
C    NPRIX=1;2: NOCH MEHR DRUCKEN
C    NPRICN=1: CN VON VERSEM WERDEN GEDRUCKT
C    NZKAPH: ANZAHL DER PYSIKALISCHEN KANAELE
C
1000  FORMAT(20I3)
      REWIND NBAND1
      REWIND NBAND2
      READ(NBAND1) MUL,JWSL,JWSR,NZKL,EB
      IF (NZKAPH.EQ.0) NZKAPH=NZKL
      WRITE(NOUT,*) ' ES WERDEN ',NZKAPH,' KANAELE GERECHNET'
      IF(NZKAPH.GT.NZKBMA) STOP 'NZKBMA ZU KLEIN '
      IF(NZKL.GT.NZKMAX) THEN
        WRITE(NOUT,*) 'NZKMAX= ',NZKMAX,' KLEINER NZKL ',NZKL
        STOP 'NZKL'
      ENDIF
      WRITE(NOUT,1002) MUL,(MREG(MKC),MKC=1,NZOPER)
1002  FORMAT('1ES WERDEN DIE OPERATOREN',/,1X,'L = ',I2,' :      ',
     $10I3,'  BERECHNET'//)
      READ(INPUT,1005)(SNORM(K),K=1,NZKAPH)
1005  FORMAT(6E12.4)
      WRITE(NOUT,1003)(SNORM(K),K=1,NZKAPH)
1003  FORMAT(/1X,'SNORM ',/1X,10E13.6)
      WRITE(NOUT,1004) JWSR
1004  FORMAT(/1X,'J(BIND) = ',I3,'/2')
      PI=3.141592654
      DF(1)=1.
      DO 200 I=2,12,2
200   DF(I+1)=DBLE(I+1)*DF(I-1)
C   DF(I) = I!!
      F=DBLE(MUL+1)/DBLE(MUL)/(DF(2*MUL+1)**2)
      HC=197.32858
      F=F*8.*PI*HC/137.03604
      DO 10 K=1,NZKL
      READ(NBAND1) LWERT(K),(JWERT(KH,K),KH=1,3),REDM(K),MASSE(1,K)
     $,MASSE(2,K),MLAD(1,K),MLAD(2,K)
       READ(NBAND2) MHH,LWERTH,MASSE1,MASSE2,MLAD1,MLAD2,JWERTH,
     *              JWSH,WNORM(K)
      WRITE(NOUT,*) 'kanal',K
      if(k.eq.1) then
         wh1=1./(wnorm(1)**2 *redm(1)**1.5)
         const=snorm(1)/wh1
      endif	 
      if(k.gt.nzkaph) snorm(k)=const/(wnorm(k)**2 *redm(k)**1.5)
       IF(ABS(LWERTH-LWERT(K))+ABS(JWERTH-JWERT(3,K))+
     *    ABS(MASSE(1,K)-MASSE1)+ABS(MASSE(2,K)-MASSE2)-ABS(MLAD(1,K)-
     *    MLAD1)+ABS(MLAD(2,K)-MLAD2)+ABS(JWSL-JWSH).EQ.0) GOTO 10
        WRITE(NOUT,1007) K
1007  FORMAT(' KANAL NR.',I4,' STIMMT BEI VERMAT UND ENEMS NICHT',
     *       ' UEBEREIN, VERMAT KANAL:')
	WRITE(NOUT,*) 'MHH,LWH,M1,M2,L1,L2,JH,JWH',
     * MHH,LWERTH,MASSE1,MASSE2,MLAD1,MLAD2,JWERTH,JWSH
        STOP 5
10    CONTINUE
      READ(NBAND1) IZQ(NZKL+1),IZPWM,(IZP(I),IZQ(I),I=1,NZKL)
      DO 12 K=1,NZKL
      READ(NBAND1)(NAR(M,K),M=1,IZPWM)
      READ(NBAND1)(PAR(M,K),M=1,IZPWM)
12    CONTINUE
      MM=IZQ(NZKL+1)
      IF(MM.GT.NDIM) THEN
        WRITE(NOUT,*) ' MM = ',MM,' GREATER DIMENSION ',NDIM
        STOP ' NDIM'
      ENDIF
      READ(NBAND1) ((OPW(K,MKC),K=1,MM),MKC=1,NZOPER)
      READ(NBAND2)(EBIN(K),K=1,NZKAPH)
1001  FORMAT(10E12.4)
      READ(NBAND2) NZAEN,(EW(MH),MH=1,NZAEN)
      READ(INPUT,1000) (ISTRN(K),K=1,NZKAPH),ISRB
C
C
      DO 1 MWEE=1,NZAEN
      ENERGI=EW(MWEE)+EBIN(1)
      NCN = 0
      DO 24 K=1,NZKAPH
      EMOM(K)=ENERGI-EBIN(K)
      EK(K)=(ABS(EMOM(K))*REDM(K)/H2M)**.5
      write(nout,*) 'k,ek',k,ek(k)
      HQOM(K)=EBIN(K)+EMOM(K)-EB
24    EKG(K)=HQOM(K)/HC
      I=0
      DO 90    K = 1,NZKAPH
      IF(EMOM(K).LE.0.) GOTO 5
       I = I + 1
      NENTK(K) = 1
      GO TO 90
    5 NENTK(K) =0
   90 CONTINUE
      NO = I
      WRITE(NOUT, 1078)
1078  FORMAT(1H1)
      WRITE (NOUT,1010) ENERGI
      DO 6 K=1,NZKAPH
      MLA1=MLAD(1,K)+1
      MLA2=MLAD(2,K)+1
      IF(NENTK(K).GT.0) GOTO 9
      WRITE (NOUT,1012) K,EMOM(K),HQOM(K)
     1 ,MASSE(1,K),IFRAK(MLA1),MASSE(2,K),
     2  IFRAK(MLA2),LWERT(K),JWERT(3,K) ,JWSL
      GO TO 6
9     WRITE (NOUT,1013) K,EMOM(K),HQOM(K)
     1      ,MASSE(1,K),IFRAK(MLA1),MASSE(2,K),
     2   IFRAK(MLA2),LWERT(K),JWERT(3,K) ,JWSL
    6 CONTINUE
 1010 FORMAT (////20H ENERGIE IM SYSTEM =,F10.4,4H MEV/)
 1012 FORMAT (4H DER,I3,25H TE KANAL IST GESCHLOSSEN /
     126H DIE KANALENERGIE BETRAEGT,F14.8,17H MEV, E(GAMMA) = ,F14.8
     2,4H MEV/,
     4  I2,A2,2H -,I2,A2,7H     L=,I1,5H   S=,I1,3H /2,
     5     5H   J=  ,I1,3H /2  /)
 1013 FORMAT (4H.DER,I3,19H TE KANAL IST OFFEN /
     126H DIE KANALENERGIE BETRAEGT,F14.8,17H MEV, E(GAMMA) = ,F14.8
     2,4H MEV ,2X,
     4  I2,A2,2H -,I2,A2,7H     L=,I3,5H   S=,I3,3H /2,
     5     5H   J=  ,I3,3H /2  /)
C
      DO 52 MTD=1,2
      READ(NBAND2)((A(K,L,MTD),L=1,NO),K=1,NO)
      IF(NPRI.LE.2) GOTO 52
      WRITE(NOUT,1001)((A(K,L,MTD),L=1,NO),K=1,NO)
52    CONTINUE
C
C
      DO 51 MTD=1,2
      NCN=0
      DO 30 K=1,NO
      DO 30 MKC=1,NZOPER
30    OPWERT(K,MKC)=0.
      WRITE(NOUT,1020) METH(MTD)
1020  FORMAT(///1X,A10)
      IF(NPRI.LE.3) GOTO 23
      WRITE(NOUT,1009)
1009  FORMAT(1X,'VON VERSEM :  CWF')
23    DO 26 K=1,NO
      L=0
      DO 27 LL=1,NZKL
      L=L+1
      IZPW2=IZP(LL)
      READ(NBAND2)(CWF(K,L,M),M=1,IZPW2)
      IF(NPRI.LE.3) GOTO 27
      WRITE(NOUT,1014) K,L
1014  FORMAT(1X,2I3)
      WRITE(NOUT,1015)(CWF(K,L,M),M=1,IZPW2)
1015  FORMAT(1X,10E13.6)
27    CONTINUE
26    CONTINUE
      IF(NPRIX.GT.1) WRITE(NOUT,1097)
1097  FORMAT(1X,'OPWERT, CWF*OPW, CWF, OPW')
80    DO 100 MKC=1,NZOPER
      IF(MREG(MKC).LE.0) GOTO 100
      DO 2 K=1,NO
      ILL=0
ccc   DO 3 L=1,NZKL
      DO 3 L=1,NZKAPH
      IZPW2=IZP(L)
      DO 3 I=1,IZPW2
      DREPO=LWERT(L)+NAR(I,L)
      ILL=ILL+1
      OPWERT(K,MKC)=OPWERT(K,MKC)+CWF(K,L,I)*OPW(ILL,MKC)/
     $SQRT(SNORM(L)*REDM(L)**(DREPO+3.))
C   UMRECHNUNG VON JACOBI- AUF RELATIVKOORD.
      IF(NPRIX.LT.2) GOTO 3
      IF(OPW(ILL,MKC).EQ.0.) GOTO 3
      AA=CWF(K,L,I)*OPW(ILL,MKC)
      WRITE(NOUT,1099) L,I,ILL,OPWERT(K,MKC),AA,CWF(K,L,I),OPW(ILL,MKC)
1099  FORMAT(1X,' L, I, ILL:', 3I8,4(5X,E12.6))
3     CONTINUE
      OPWERT(K,MKC)=OPWERT(K,MKC)*SQRT(4./DBLE(JWSL+1))
      IF(MKC.EQ.1.OR.MKC.EQ.2) OPWERT(K,MKC)=OPWERT(K,MKC)/EKG(K)
      IF(MKC.EQ.3.OR.MKC.EQ.4) OPWERT(K,MKC)=OPWERT(K,MKC)*EKG(K)
      IF(NPRIX.LE.0) GOTO 2
      WRITE(NOUT,1098)K,MKC,OPWERT(K,MKC)
1098  FORMAT(1X,'OPWERT(KANAL',I2,',OPERATOR',I2,')', 2E12.6)
2     CONTINUE
C
100   CONTINUE
      IF(NPRI.LE.0) GOTO 41
      WRITE(NOUT,1019)
1019  FORMAT(/1X,'KANAL    RED. MATRIXELEMENTE <KANAL//OP//BOUND>'
     $,'  1 BIS 10')
      DO 50 K=1,NO
      WRITE(NOUT,1018) K,(OPWERT(K,MKC),MKC=1,NZOPER)
50    CONTINUE
1018  FORMAT(1X,I3,/,(1X,10E12.4,/))
41    CONTINUE
      DO 40 K=1,NO
      ARG=  PI/2*(-LWERT(K)+MUL)
      CPHASE=EXP(CMPLX(RNULL,ARG)) *
     $ (-1.)**((JWSL+JWSR)/2+MUL)
C     PHASE=I**(-LWERT(K)).
C     DIESE PHASE STAMMT AUS DER STREUWELLENFKT.
C     PHASE=(-1)**[(JWSL+JWSR)/2+MUL] * I**MUL
C     DIES STAMMT AUS DEM VERGLEICH DER DIFF. WQ VON WELLER UND T.M.
      KS1=IABS(JWERT(1,K)-JWERT(2,K))
      INR=JWERT(3,K)
c      CALL ISPNR(JWERT(3,K),KS1,INR)
      IF (MREG(1)+MREG(2)+MREG(3)+MREG(4).LE.0) GOTO 42
      OPE=OPWERT(K,1)+OPWERT(K,2)+OPWERT(K,3)+OPWERT(K,4)
      OPE=OPE*CPHASE
      IP=1
      CALL PHASOD(OPE,TE,PE)
42    IF(MREG(5)+MREG(6)+MREG(7)+MREG(8).LE.0) GOTO 43
      OPM=OPWERT(K,5)+OPWERT(K,6)+OPWERT(K,7)+OPWERT(K,8)
      OPM=OPM*CPHASE*CMPLX(RNULL,-REINS)
C     FAKTOR -I BEI DEN MAGNETISCHEN OPERATOREN VON T.M.
      IP=2
      CALL PHASOD(OPM,TE,PE)
43    IF(MREG(9)+MREG(10).LE.0) GOTO 44
      OPES=OPWERT(K,9)+OPWERT(K,10)+OPWERT(K,3)+OPWERT(K,4)
      OPES=OPES*CPHASE
      IP=1
      CALL PHASOD(OPES,TES,PES)
44    CONTINUE
      XG=EKG(K)
      XX=SQRT(F*(XG**(2*MUL+1)))
      IF(MREG(1)+MREG(2)+MREG(3)+MREG(4)+MREG(5)+MREG(6)+
     *MREG(7)+MREG(8).LE.0) GOTO 140
       TE=TE*XX
      WRITE(NOUT,2000) EMOM(K),HQOM(K),ISTRN(K),ISRB,INR,IP,JWSL
     1  ,LWERT(K),MUL,TE,PE,MTD
      WRITE(15,2000) EMOM(K),HQOM(K),ISTRN(K),ISRB,INR,IP,JWSL
     1  ,LWERT(K),MUL,TE,PE,MTD
2000  FORMAT(2F12.5,7I3,F12.9,F12.4,I3)
 140  IF(MREG(9)+MREG(10).LE.0) GOTO 40
       TES=TES*XX
      WRITE(NOUT,2000) EMOM(K),HQOM(K),ISTRN(K),ISRB,INR,IP,JWSL
     1  ,LWERT(K),MUL,TES,PES,MTD+4
      WRITE(15,2000) EMOM(K),HQOM(K),ISTRN(K),ISRB,INR,IP,JWSL
     1  ,LWERT(K),MUL,TES,PES,MTD+4

40    CONTINUE
51    CONTINUE
C
C
1     CONTINUE
C
C
      WRITE(NOUT,1040)
1040  FORMAT(1X,'ENDE DER RECHNUNG VON ELMAT')
      STOP
      END
      SUBROUTINE ISPNR(J,K,I)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
       INCLUDE 'par/elmat'
      DO 10 IK=2,10,2
      IF(J.EQ.K+IK-2) GOTO20
  10  CONTINUE
      STOP
  20  I=IK/2
      RETURN
      END
      SUBROUTINE PHASOD(S,T,P)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
       INCLUDE 'par/elmat'
      COMPLEX *16 S
      DATA PI/3.141592654/
      T=ABS(S)
      C=DREAL(S)
      IF(C.NE.0.) GOTO 706
      B=0.
      GOTO 710
706   B=(90./PI)*ATAN(IMAG(S)/C)
      IF(C.GT.0.) GOTO 710
      IF(B.GT.0) GOTO 704
      B=B+90.
      GOTO 710
704   B=B-90.
710    P=B
      RETURN
      END
