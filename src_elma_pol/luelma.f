      PROGRAM LUELMA
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C     PROGRAM LUELMA PREPARES THE CALCULATION OF THE REDUCED
C     SPACE MATRIX ELEMENTS
C  
C   FUER ELEKTROMAGNETISCHE UEBERGAENGE
C
C
C                      T.M. 1985
C
C     LETZTE AENDERUNG H.M.H. 26.2.88 IQM
C     29.1.03 AUF PARAMETER UMGESTELLT UND KRECH EINGEBAUT H.M.H.
C
      INCLUDE 'par/luelma'
C
      PARAMETER (NZLMAX=3*2*(NZCMAX-1))
C     NZLMAX: MAXIMALE ANZAHL DER DREHIMPULSE IM MATR.-EL
C
      PARAMETER (PI=3.1415926536)
C
      COMMON /COMY/D(100)
C     LOGARITHMS OF FACULTIES FOR ANGULAR MOMENTUM ROUTINES
      COMMON /SDA/CW(NDIMCW), MVM(NZIQMA,NDIMCW), KAUS,
     *            N5, IQ, NKAPO
C     WIRD IN SDA1 BENOETIGT
C     N5: ANZAHL DER DREHIMPULSE IM MATR.-EL.
C     IQ: ANZAHL DER SIGMA-FAKTOREN
C     MVM: VERSCHLUESSELTE INDIZES DER SIGMAFAKTOREN
C     CW: VORFAKTOR DES PRODUKTS DER SIGMAFAKTOREN
C     KAUS: STEUERT ZUSAETZLICHEN AUSDRUCK IN SDA1
C
      COMMON/MVA/KV(2*(NZCMAX-1)),KMV(2*(NZCMAX-1),NDIMK),
     *            MMOEG,MKAUS
      COMMON /FAKTOR/          F8(2*LWMAX+1) 
C
      COMMON/ORD/EPO(NDIMOR),KZAHL, JQ(NDIMOR),
     *           INDPO(NDIMOR), MVK(NZIQMA,NDIMOR),
     *           INDEX(NDIMOR)
C     
      DIMENSION IENT(NZOPER), KRECH(NZFMAX,NZLWMA)
      DIMENSION NZLW(NZFMAX), NZC(NZFMAX), NZPO(NZFMAX)
      DIMENSION LW(2*NZCMAX-3,NZLWMA,NZFMAX),
     *          KP(NZCMAX-1,NZPOMA,NZFMAX)
C     LW(I,J,K): I-TE DREHIMPULSWERT DER J-TEN DREH-
C                IMPUSSTRUKTUR DER K-TEN ZERLEGUNG
C     KP(I,J,K): WIE LW FUER POLYNOME
C     
      DIMENSION LV(NZLMAX), MV(NZLMAX), LVMOD(NZLMAX)
      DIMENSION LKOP1(NZCMAX-2,NZLWMA,NZFMAX),
     *          LKOP2(NZCMAX-2,NZLWMA,NZFMAX)
      DIMENSION C(NDIMC), MW(NDIMC,NZLMAX)
      DIMENSION KWL(2*NZCMAX-3,5), KWR(2*NZCMAX-3,5) 
      DIMENSION F4(4), F5(4), F3(2*LWMAX+1), F7(2*LWMAX+1)
C
C     DIMENSION NZPO(NZF),F3(2*LW-1),F7(2*LW-1),F8(2*LW-1)
C     DIMENSION LW(2*NZC-3,NZLW,NZF),KP(NZC-1,NZLW,NZF),F4(KP+1),F5(KP+1)
C     DIMENSION IENT(NOP),NZLW(NZF),NZC(NZF),LKOP1(NZC-2,NZLW,NZF)
C     DIMENSION LKOP2(NZC-2,NZLW,NZF),LV(3*2*(NZC-1)+1),MV(3*2*(NZC-1)+1)
C     DIMENSION KWL(2*NZC-3,5),KWR(2*NZC-3,5),KV(2*NZC-2)
C     DIMENSION C(NDIMC),MW(NDIMC,3*2*(NZC-1)+1),LVMOD(3*2*(NZC-1)+1)
C     DIMENSION CW(NDIMCW),MVM(SUMME L,NDIMCW),KMV(2*NZC-2,NDIMK)
C     START VALUES
      DATA F4/1.,-2.,2.,-2./,F5/1.,.25,.0625,.015625/
C     F4(I+1)=(-)**I (2-DELTA(I,0)),F5(I+1)=1/4**I
C
C
      OPEN(UNIT=5,FILE='INLU',STATUS='OLD')
      OPEN(UNIT=6,FILE='OUTPUT',STATUS='UNKNOWN')
      OPEN(UNIT=9,FILE='LUOUT',STATUS='UNKNOWN',FORM='UNFORMATTED')
C
      NKAPO=0
      MMOEG=1
C     NO POLYNOMIALS FOR NKAPO=0
C
C
      D(1)=0.
      D(2)=0.
      DO 100 I=2,99
      HMH=I
100   D(I+1)=LOG(HMH)+D(I)
      INPUT = 5
      FT=SQRT(1./(4.*3.1415926536))
      F7(1)=FT
      F3(1)=1.
      F8(1)=1.
      FH=1.
      DO 102 NH=2,2*LWMAX+1
      FSH=NH-1
      FH=FSH*FH
      F3(NH)=SQRT(FH)
      F8(NH)=1./FH
      FT=FT*.5
102   F7(NH)=FT*SQRT(2.*FSH+1.)
C     F3(I)=SQRT((I-1)!),F8(I)=1/(I-1)!,F7(I)=SQRT(2I-1)/2**(I-1)/SQRT(4PI)
C
C     READ ANGULAR MOMENTUM AND POLYNOMIAL STRUCTURES AND WRITE ON TAPE
1000  FORMAT(20I3)
      READ (INPUT,1000) NBAND,LAUS,KAUS,MAUS,MKAUS
C     NBAND=9, TAPENUMBER FOR RESULTS
C     KAUS =1, ADDITIONAL OUTPUT IN SDA1,LAUS=1,MAUS=1 ADDITIONAL
C     OUTPUT IN MAIN,MKAUS=1 ADDITIONAL OUTPUT IN MVAL
      READ (INPUT,1000) (IENT(NOP),NOP=1,NZOPER)
C     IENT =1 DETERMINES OPERATORS TO BE CALCULATED
C     IENT(1..4) = EL-BAHN-,EL-SPIN-,ML-BAHN-,ML-SPIN-OPERATOR
      READ (INPUT,1000) NZF,MUL
C     NZF NUMBER OF ZERLEGUNGEN, MUL MULTIPOLARITAET
      DO 110 MH=1,NZF
110   READ (INPUT,1000) NZLW(MH),NZC(MH),NZPO(MH)
C
      REWIND NBAND
      WRITE(NBAND)NZF,MUL,(NZLW(MH),NZC(MH),MH=1,NZF),
     *      (IENT(MH),MH=1,NZOPER), (NZPO(MH),MH=1,NZF)
C
C     READ AND WRITE ACTUAL ANGULAR MOMENTA AND POLYNOMIALS
C
      WRITE(NOUT,1060)
1060  FORMAT(1H1)
      WRITE (NOUT,1001) NZF
1001  FORMAT(' ZAHL DER ZERLEGUNGEN',I3)
      IF(NZF.GT.NZFMAX) STOP 1
      DO 140 MH=1,NZF
      WRITE (NOUT,1002) MH
1002  FORMAT(//,'0 ZERLEGUNG',I5)
      N1=NZC(MH)-1
      IF(N1.GT.3) STOP 2
      N3=N1-1
C     NUMBER OF ANGULAR MOMENTA
      N2=NZLW(MH)
      WRITE (NOUT,1003) NZC(MH),NZLW(MH)
1003  FORMAT(' ANZAHL DER CLUSTER',I5,
     1 '      ANZAHL DER DREHIMPULSSTRUKTUREN',I5)
      DO 130 L=1,N2
      WRITE (NOUT,1004) L
1004  FORMAT(/,' DREHIMPULSSTRUKTUR',I5)
      READ (INPUT,1000)(LW(M,L,MH),M=1,N1),KRECH(MH,L)
C       LW(M,L,MH): M-TE DREHIMPULSWERT DER L-TEN DREHIMPULSSTRUKTUR
C                   DER MH-TEN ZERLEGUNG
C       KRECH(MH,L): ES WERDEN NUR SOLCHE MATRIXELEMENTE BERECHNET
C                    BEI DENEN DAS PRODUKT VON KRECH NICHT GROESSER NULL IST
C
      WRITE (NOUT,1005) (M,LW(M,L,MH),M=1,N1)
1005  FORMAT(I4,' TER DREHIMPULS =',I3,' RELATIVDREHIMPULS')
      IF(N3.LE.0) GOTO 130
      DO 120 M=1,N3
      MM=N1+M
      READ(INPUT,1000) LW(MM,L,MH),LKOP1(M,L,MH),LKOP2(M,L,MH)
C     LW COUPLED ANGULAR MOMENTUM FROM LKOP1 AND LKOP2
      IF(N3.GT.1) GOTO 116
      LKOP1(1,L,MH) =1
      LKOP2(1,L,MH) =2
116   WRITE (NOUT,1006) MM,LW(MM,L,MH),LKOP1(M,L,MH),LKOP2(M,L,MH)
1006  FORMAT(I4,' TER DREHIMPULS =',I3,' GEKOPPELT AUS DREHIMPULS',
     1 I3,' UND DREHIMPULS',I3)
120   CONTINUE
130   CONTINUE
      MPO=NZPO(MH)
      NZPO(MH)=MAX0(NZPO(MH),1)
      IF(MPO.GT.0)GOTO 124
      DO 122 KH=1,N1
122   KP(KH,1,MH)=0
      MPO=1
      GOTO 127
124   DO 126 KH=1,MPO
126   READ(INPUT,1000) (KP(M,KH,MH),M=1,N1)
      NKAPO=1
127   DO 128 KH=1,MPO
128   WRITE (NOUT,1007) KH,(KP(M,KH,MH),M=1,N1)
1007  FORMAT(' BEI DER',I3,' TEN POLYNOM FUNKTION SIND POLYNOME DER'
     1 ,' ORDNUNG',3I3)
      N4=N1+N3
140   WRITE (NBAND)((LW(M,L,MH),M=1,N4),L=1,N2),
     1 ((KP(M,KH,MH),M=1,N1),KH=1,MPO)
C     END OF INPUT
C     START OF CALCULATION OF MATRIX ELEMENTS
      DO 900 MFL=1,NZF
C     LOOP ZERLEGUNGEN LEFT  ------------------------------------------------------------ ECCE lupo: NZL(eft)
      N1=NZC(MFL)-1
      N8=N1+N1-1
      N10=N1+1
      N2=NZLW(MFL)
      K1=NZPO(MFL)
      DO 890 MFR=1,MFL
      WRITE(NOUT,*) ' ZWISCHEN ZERLEGUNG',MFL,' UND',MFR,' WIRD ',
     *             'GERECHNET'
C     LOOP ZERLEGUNGEN RIGHT
      N3=NZC(MFR)-1
      N7=N3+1
      N9=N3+N3-1
      N4=NZLW(MFR)
      N5=N1+N3
      NH5=(2*NKAPO+1)*N5
      NH6=NH5+1
      K2=NZPO(MFR)
      DO 880 NOP=1,4
C     LOOP OF OPERATORS
      IZREK=0
      IF(IENT(NOP).LE.0)GOTO 880
      GOTO(601,602,603,604),NOP
601   WRITE(NOUT,1031) MUL
      GOTO 605
602   WRITE(NOUT,1032) MUL
      GOTO 605
603   WRITE(NOUT,1033) MUL
      GOTO 605
604   WRITE(NOUT,1034) MUL
605   CONTINUE
1031  FORMAT(/,1X,'E',I2,'-BAHN-OPERATOR WIRD GERECHNET')
1032  FORMAT(/,1X,'E',I2,'-SPIN-OPERATOR WIRD GERECHNET')
1033  FORMAT(/,1X,'M',I2,'-BAHN-OPERATOR WIRD GERECHNET')
1034  FORMAT(/,1X,'M',I2,'-SPIN-OPERATOR WIRD GERECHNET')
      L3=1
160   KZAHL =0
C      FOR REDUCING ANGULAR MOMENTA
      GOTO (162,164,162,164),NOP
162   WRITE (NOUT,1009)L3
1009  FORMAT(/,' DER',I3,' TE DREHIMPULS WIRD ERNIEDRIGT')
164   DO 800 LL=1,N2
C      LOOP DREHIMPULSSTRUKTUREN LEFT
      NPWL=0
      DO 166 KH=1,N1
166   NPWL=NPWL+LW(KH,LL,MFL)
      NPWLI=NPWL
      NPWL=(-1)**NPWL
      DO 790 LR=1,N4
C     LOOP DREHIMPULSSTRUKTUREN RIGHT
      IF (KRECH(MFL,LL)*KRECH(MFR,LR).GT.0) GOTO 790
C     HIER WERDEN IM INPUT BESTIMMTE DREHIMPULSSTRUKTUREN UEBERSPRUNGEN.
      NPWR=0
      DO 168 KH=1,N3
168   NPWR=NPWR+LW(KH,LR,MFR)
      IF((NPWLI+NPWR+MUL)/2.GT.NZIQMA) STOP 'NZIQMA ZU KLEIN'
      NPWR=(-1)**NPWR
C     PARITY CHECK
      IF(NOP.GE.3)GOTO 501
      NPAROP=(-1)**MUL
      GOTO 502
501   NPAROP=(-1)**(MUL+1)
502   IF(NPWL*NPWR.NE.NPAROP)GOTO 790
      IDL=LW(N8,LL,MFL)-LW(N9,LR,MFR)
C     IDL TOTAL ANGULAR MOMENTUM DIFFERENCE
C     FACTOR FOR REDUCED MATRIXELEMENT
      BL=LW(N8,LL,MFL)
      F=1.
      IMUL= MUL   
      IF(NOP.EQ.4) IMUL= MUL -1
      GOTO (182,190,182,190),NOP
182    IZUS1=1
      IZUS2=-1
190   Y=CLG(2*LW(N9,LR,MFR),2*IMUL,2*LW(N8,LL,MFL),2*LW(N9,LR,MFR),
     *      2*IDL)
c      WRITE (NOUT,*) 2*LW(N9,LR,MFR),2*IMUL,2*LW(N8,LL,MFL),
c     *               2*LW(N9,LR,MFR),2*IDL
C     MAXIMAL COUPLING
      IF(Y.NE.0.) GOTO 192
      F=0.
      IF(MAUS.LE.0) GOTO 790
      WRITE (NOUT,1010) LW(N9,LR,MFR),MUL,LW(N8,LL,MFL)
1010  FORMAT(' KOPPLUNG LIEFERT NULL',3I5)
      GOTO 790
192   F=F*SQRT(2.*BL+1.)/Y
C
C     DETERMINE ANGULAR MOMENTA AND FACTOR
C
      MV(NH6)=IDL
      DO 196 MH=1,N1
      LV(MH)=LW(MH,LL,MFL)
      KH=LV(MH)
196   F=F*F7(KH+1)/FLOAT((-1)**KH)
C     FACTOR LEFT SIDE
      DO 198 MH=1,N3
      MMH=N1+MH
      LV(MMH)=LW(MH,LR,MFR)
      KH=LV(MMH)
198   F=F*F7(KH+1)/FLOAT((-1)**KH)
C     FACTOR RIGHT SIDE
      GOTO(503,504,504,503),NOP
503   LV(NH6)=MUL-1
      GOTO 505
504   LV(NH6)=MUL
505   A=FLOAT(LV(NH6))
C     LETZTER DREHIMPULS IST OPERATOR
      AMUL=FLOAT(MUL)
      F=F*F7(LV(NH6)+1)/FLOAT((-1)**LV(NH6))
      DO 780 KZL=1,K1
C     LOOP POLYNOMIALS LEFT
      DO 199 KHH=1,N1
199   KV(KHH)=KP(KHH,KZL,MFL)
      DO 770 KZR=1,K2
C     LOOP POLYNOMIALS RIGHT
      IF(NKAPO.EQ.0) GOTO 203
      DO 200 KHH=1,N3
      MHH=KHH+N1
200   KV(MHH)=KP(KHH,KZR,MFR)
      DO 201 KHH=1,N5
      IH1=KHH+N5
      IH2=IH1+N5
      LV(IH1)=KV(KHH)
201   LV(IH2)=KV(KHH)
C     HIGH L-VALUES ARE POLYNOMIALS
      CALL MVAL(N5)
203    GOTO(204,210,204,210),NOP
C     DIFFERENTIAL OPERATOR CHECK IF L-VALUE REDUCABLE
204   IF(LV(L3).LE.0) GOTO 770
210   NP=0
      IF(MAUS.GT.0) WRITE(NOUT,1011) KZL,LL,MFL,NPWL,KZR,LR,MFR,NPWR
1011  FORMAT(/,' MATRIXELEMENT ZWISCHEN ',
     1 'DEM',I3,' TEN POLYNOM DER',I3,
     2 ' TEN DREHIMPULSSTRUKTUR DER',I3,
     3 ' TEN ZERLEGUNG MIT PARITAET',I3,' UND',/,24X,
     4 'DEM',I3,' TEN POLYNOM DER',I3,
     5 ' TEN DREHIMPULSSTRUKTUR DER',I3,
     6 ' TEN ZERLEGUNG MIT PARITAET',I3)
C     SET LEFT M-VALUES TO MINUS L
      DO 222 MH=1,N1
222   MV(MH)=-LV(MH)
      DO 224 MH=N10,N5
224    MV(MH)=LV(MH)
C     DETERMINE ALL POSSIBLE M-VALUES OF LEFT SIDE WITH SUM=TOTAL L
C     AND CALCULATE FACTOR F1
230   F1=1.
      ML=0
      DO 232 MH=1,N1
      KWL(MH,1)=MV(MH)
232   ML=ML+MV(MH)
      IF(ML+LW(N8,LL,MFL).NE.0) GOTO 240
      IF(N1.LE.1) GOTO 254
C     DETERMINE IF COUPLING OF SEVERAL ANGULAR MOMENTA POSSIBLE
      DO 234 MH=N10,N8
      MMH=MH-N1
      MM1=LKOP1(MMH,LL,MFL)
      MM2=LKOP2(MMH,LL,MFL)
      KWL(MH,4)=LW(MM1,LL,MFL)
      KWL(MH,5)=LW(MM2,LL,MFL)
      KWL(MH,2)=KWL(MM1,1)
      KWL(MH,3)=KWL(MM2,1)
      KWL(MH,1)=KWL(MH,2)+KWL(MH,3)
      IF(IABS(KWL(MH,1)).GT.LW(MH,LL,MFL) ) GOTO 240
234   CONTINUE
      GOTO 250
C     INCREASE MVALUES
240   MH=N1
242   IF(MV(MH)-LV(MH)) 244,246,246
244   MV(MH)=MV(MH)+1
      GOTO 230
246   MV(MH)=-LV(MH)
      MH=MH-1
      IF(MH) 400,400,242
C     ALL COMBINATIONS FOUND
C     CALCULATE FACTOR,D(L,M) AND CLEBSCH
250   ML=MV(1)
      DO 252 MH=2,N1
      M1=LV(MH)+MV(MH)+1
      M2=LV(MH)-MV(MH)+1
      MM=N1+MH-1
      BH1=KWL(MM,4)
      BH2=KWL(MM,5)
      BH3=LW(MM,LL,MFL)
      L4=-2*KWL(MM,2)
      L5=-2*KWL(MM,3)
      F1=F1*F3(M1)*F3(M2)*
     *    CLG(2*KWL(MM,4),2*KWL(MM,5),2*LW(MM,LL,MFL),L4,L5)
252   ML=ML+MV(MH)
254   M1=LV(1)+MV(1)+1
      M2=LV(1)-MV(1)+1
      F1=((-1)**ML)*F1*F3(M1)*F3(M2)
C     END LEFT SIDE
C     DETERMINE ALL POSSIBLE M-VALUES OF RIGHT SIDE WITH SUM = TOTAL L
C     AND CALCULATE FACTOR F2
260    F2=1.
      MR=0
      DO 262 MH=N10,N5
      MMH=MH-N10+1
      KWR(MMH,1)=MV(MH)
262   MR=MR+MV(MH)
      IF(MR-LW(N9,LR,MFR).NE.0) GOTO 270
      IF(N3.LE.1) GOTO 284
      DO 264 MH=N7,N9
      MMH=MH-N3
      MM1=LKOP1(MMH,LR,MFR)
      MM2=LKOP2(MMH,LR,MFR)
      KWR(MH,4)=LW(MM1,LR,MFR)
      KWR(MH,5)=LW(MM2,LR,MFR)
      KWR(MH,2)=KWR(MM1,1)
      KWR(MH,3)=KWR(MM2,1)
      KWR(MH,1)=KWR(MH,2)+KWR(MH,3)
      IF(IABS(KWR(MH,1)).GT.LW(MH,LR,MFR) ) GOTO 270
264   CONTINUE
      GOTO 280
C     REDUCE M-VALUES
270   MH=N5
272   IF(MV(MH)+LV(MH)) 276,276,274
274   MV(MH)=MV(MH)-1
C     NEXT COMBINATION
      GOTO 260
276   MV(MH)=LV(MH)
      MH=MH-1
      IF(MH-N1) 240,240,272
C     LOOK FOR NEW COMBINATION ON LEFT SIDE
C     CALCULATE FACTOR RIGHT SIDE,D(L,M) AND CLEBSCH
280   DO 282 MH=2,N3
      MMH=MH+N1
      M1=LV(MMH)+MV(MMH)+1
      M2=LV(MMH)-MV(MMH)+1
      NN=N3+MH-1
      BH1=KWR(NN,4)
      BH2=KWR(NN,5)
      BH3=LW(NN,LR,MFR)
      L4=2*KWR(NN,2)
      L5=2*KWR(NN,3)
282   F2=F2*F3(M1)*F3(M2)*
     *    CLG(2*KWR(NN,4),2*KWR(NN,5),2*LW(NN,LR,MFR),L4,L5)
284   M1=LV(N10)+MV(N10)+1
      M2=LV(N10)-MV(N10)+1
      F2=F2*F3(M1)*F3(M2)
C     END RIGHT SIDE
C     CHECK DIMENSION OF C AND MW
      IF(NP.LT.NDIMC) GOTO 286
      WRITE (NOUT,1012) NP,NDIMC
1012  FORMAT(' DIMENSION VON C UND MW ZU KLEIN',2I5)
      STOP 6
C     LOOP M-VALUES OF POLYNOMIALS
286   DO 390 IMH=1,MMOEG
      FAPO=1.
      IF(NKAPO.EQ.0) GOTO 290
      DO 288 JH=1,N5
      JH1=JH+N5
      JH2=JH1+N5
      MV(JH1)=KMV(JH,IMH)
      KM=MV(JH1)
      KL=LV(JH1)
      FAPO=FAPO*F4(KM+1)/F8(KL+KM+1)*F5(KL+1)/F8(KL-KM+1)
288   MV(JH2)=-KMV(JH,IMH)
C     BERECHNUNG DER FAKTOREN C UND NOTIEREN VON C UND DEN M-WERTEN
290   FPO2=F2*FAPO
      GOTO (296,292,296,292),NOP
292   NP=NP+1
      DO 294 KH=1,NH6
294   MW(NP,KH)=MV(KH)
      LM=LV(NH6)-MV(NH6)+1
      LP=LV(NH6)+MV(NH6)+1
      C(NP)=F1*FPO2*F*F3(LM)*F3(LP)
      GOTO 390
296   MV(NH6)=IDL-1
      IF(IABS(MV(L3)+IZUS1).GE.LV(L3))GOTO 506
      NP=NP+1
      DO 308 KH=1,NH6
308   MW(NP,KH)=MV(KH)
      LM=LV(NH6)-MV(NH6)+1
      LP=LV(NH6)+MV(NH6)+1
      AZ=FLOAT(MW(NP,NH6))
      MW(NP,L3)=MW(NP,L3)+IZUS1
      CL=CLG(2*LV(NH6),2,2*MUL,2*MW(NP,NH6),2)
      IF(CL.NE.0.0) GOTO320
      NP=NP-1
      IF(MAUS.LE.0) GOTO 506
      WRITE(NOUT,1050)CL
1050  FORMAT(1X,'CG(LV(NH6),1,MUL) = ',F8.4)
      GOTO 506
320   C(NP)=CL*F1*FPO2*F*F3(LM)*F3(LP)*SQRT(2.)
506   MV(NH6)=IDL
      IF(IABS(MV(L3)).GE.LV(L3))GOTO 507
      NP=NP+1
      DO 310 KH=1,NH6
310   MW(NP,KH)=MV(KH)
      LM=LV(NH6)-MV(NH6)+1
      LP=LV(NH6)+MV(NH6)+1
      AZ=FLOAT(MW(NP,NH6))
      CL=CLG(2*LV(NH6),2,2*MUL,2*MW(NP,NH6),0)
      IF(CL.NE.0.0) GOTO 322
      NP=NP-1
      IF(MAUS.LE.0) GOTO 507
      WRITE(NOUT,1050)CL
      GOTO 507
322   C(NP)=-CL*F1*FPO2*F*F3(LM)*F3(LP)*2.
507   MV(NH6)=IDL+1
      IF(MV(L3)+IZUS2+LV(L3).LE.0)GOTO 390
      NP=NP+1
      DO 312 KH=1,NH6
312   MW(NP,KH)=MV(KH)
      LM=LV(NH6)-MV(NH6)+1
      LP=LV(NH6)+MV(NH6)+1
      AZ=FLOAT(MW(NP,NH6))
      MW(NP,L3)=MW(NP,L3)+IZUS2
      CL=CLG(2*LV(NH6),2,2*MUL,2*MW(NP,NH6),-2)
      IF(CL.NE.0.0) GOTO 324
      NP=NP-1
      IF(MAUS.LE.0) GOTO 390
      WRITE(NOUT,1050) CL
      GOTO 390
324   C(NP)=CL*F1*FPO2*F*F3(LM)*F3(LP)*SQRT(2.)
390   CONTINUE
      GOTO 270
C     ALL POSSIBLE M-COMBINATIONS FOUND
C     LOOK FOR POSSIBLE REPRESENTATIONS IN SDA1
400   NQ=NH6
      GOTO(408,402,408,402),NOP
402   IF(NP.LE.0) GOTO 770
      CALL SDA1(LV,MW,NP,C,MZAHL,NQ)
      GOTO 412
408   DO 410 MH=1,NH6
410   LVMOD(MH)=LV(MH)
      LVMOD(L3)=LV(L3)-1
      IF(NP.LE.0) GOTO 770
      CALL SDA1(LVMOD,MW,NP,C,MZAHL,NQ)
412   IF(MZAHL.LE.0) GOTO 770
      IF(MZAHL.LE.NDIMCW) GOTO 414
      WRITE (NOUT,1013) MZAHL,NDIMCW
1013  FORMAT(' DIMENSION VON CW UND MVM ZU KLEIN',2I5)
      STOP 11
414   IF(MAUS.LE.0) GOTO 430
      WRITE (NOUT,1014) MZAHL,IQ
1014  FORMAT(' NUMBER OF MATRIX ELEMENTS',I3,' NUMBER OF ',
     1 'SIGMA-FACTORS',I3)
      WRITE(NOUT,1015)
1015  FORMAT(' MATRIX ELEM SIGMA-FACTORS')
      DO 416 MH=1,MZAHL
416   WRITE (NOUT,1016) CW(MH),(MVM(KH,MH),KH=1,IQ)
1016  FORMAT(1X,E12.4,I5,20I3)
430   CALL PORD(LL,LR,MZAHL,KZL,KZR)
      MDIM=NDIMOR
C
C     END LOOP POLYNOMIALS RIGHT
770   CONTINUE
C     END LOOP POLYNOMIALS LEFT
780   CONTINUE
C     END LOOP DREHIMPULSSTRUKTUREN RIGHT
790   CONTINUE
C     END LOOP DREHIMPULSSTRUKTUREN LEFT
800   CONTINUE
C     DETERMIN IQM
      IQM=0
      IF(KZAHL.LE.0) GOTO 803
      DO 802  KH=1,KZAHL
      IF(IQM.GE.JQ(KH)) GOTO 802
      IQM=JQ(KH)
802   CONTINUE
803   IQM=MAX0(IQM,1)
      WRITE (NBAND) KZAHL,IQM
      WRITE (NOUT,1017) KZAHL,IQM,MDIM
      IZREK=IZREK+1
1017  FORMAT(' KZAHL=',I5,'      IQM=',I4,'      DIMENSION ',I6)
      IF(KZAHL.LE.0) GOTO 808
      IF(LAUS.GT.0) WRITE (NOUT,1018)
1018  FORMAT(//,' AUSDRUCK DER TERME WIE AUF BAND')
                     CALL WRITAP(LAUS,IQM,NBAND)
808   GOTO (810,880,810,880),NOP
C     FOR DIFFERENTIAL OPERATORS,REDUCE L-VALUE
810   L3=L3+1
      IF(L3.LE.NH6) GOTO 160
      WRITE (NOUT,879) IZREK
879   FORMAT(I10,' RECORDS GESCHRIEBEN')
C     END LOOP OPERATORS
880   CONTINUE
C     END LOOP ZERLEGUNGEN RIGHT
890   CONTINUE
C     END LOOP ZERLEGUNGEN LEFT
900   CONTINUE
      WRITE (NOUT,1981)
1981  FORMAT(' ENDE DER RECHNUNG VON LUELMA')
      END
      SUBROUTINE MVAL(N5)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C     THIS SUBROUTINE CALCULATES ALL POSSIBLE MVALUE COMBINATIONS OF
C     POLYNOMIALS
      INCLUDE 'par/luelma'
C
      COMMON/MVA/KV(2*(NZCMAX-1)),KMV(2*(NZCMAX-1),NDIMK),
     *            MMOEG,MKAUS
      DIMENSION IHV(6)
      MMOEG=1
      DO 10 I=1,N5
      IHV(I)=KV(I)
10    MMOEG=MMOEG*(KV(I)+1)
      IF(MMOEG.GT.NDIMK) GOTO 199
      DO 100 I=1,MMOEG
      DO 20 NH=1,N5
20    KMV(NH,I)=IHV(NH)
      MH=N5
25    IF(IHV(MH).LE.0)  GOTO 30
      IHV(MH)=IHV(MH)-1
      GOTO 100
30    IHV(MH)=KV(MH)
      MH=MH-1
      IF(MH) 40,40,25
40     IF(I.LT.MMOEG)GOTO 200
100   CONTINUE
      IF(MKAUS.EQ.0) GOTO 160
      DO 140  MH=1,MMOEG
140   PRINT 150,(KMV(NH,MH),NH=1,N5)
150    FORMAT(' MVALUES',10I4)
160   RETURN
198   FORMAT(' DIMENSION KMV ZU KLEIN')
199   PRINT 198
200   PRINT 201,I,MMOEG,KV,IHV
201   FORMAT(' NACH',I5,' VERSUCHEN VON',I5,' AUFGEGEBEN',
     1 ' KV,IHV= ',12I3)
      STOP 100
      END
      SUBROUTINE SDA1(L,MW,NEL,C,JJ,ND)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C     THIS SUBROUTINE CALCULATES ALL POSSIBLE REPRESENTATIONS
C     AND DETERMINES WHAT SIGMA FACTORS EXIST
C     FOR EACH REDUCED MATRIXELEMENT
      INCLUDE 'par/luelma'
C
      PARAMETER (NZLMAX=3*2*(NZCMAX-1))
C     NZLMAX: MAXIMALE ANZAHL DER DREHIMPULSE IM MATR.-EL
C
C     L CONTAINS ND L-VALUES
C     MW CONTAINS ALL M-VALUE COMBINATIONS
C     NEL IS THE NUMBER OF MATRIX ELEMENTS C
C     JJ IS DETERMINED IN SDA1 AND GIVES THE NUMBER OF SIGMA COMBINATIONS
      COMMON /SDA/CW(NDIMCW), MVM(NZIQMA,NDIMCW), KAUS,
     *            N5, IQ, NKAPO
      COMMON /FAKTOR/          F8(2*LWMAX+1) 
      DIMENSION MW(NDIMC,NZLMAX), L(NZLMAX)      
      DIMENSION IVM(NZIQMA),IVN(NZIQMA),IG(NZIQMA),C(NDIMC)
      DIMENSION JVM(NZIQMA),JVN(NZIQMA),LQ(3*(2*(NZCMAX-1)))
      DIMENSION KY(3*(2*(NZCMAX-1))),KYS(3*(2*(NZCMAX-1)))
      DIMENSION MKOM(2*(NZCMAX-1)+1,2*(NZCMAX-1)+1)
      DIMENSION IWN(NZIQMA),IWM(NZIQMA)
C
C
C     CHECK PRINTOUT
      IF(KAUS.LE.0) GOTO 20
      WRITE (NOUT,1001) NEL,ND
1001  FORMAT(' ANZAHL DER MATRIXELEMENTE',I4,'ANZAHL DER',
     1 ' DREHIMPULSE',I4)
      WRITE (NOUT,1002) (L(KH),KH=1,ND)
1002  FORMAT(' DREHIMPULSE',19I3)
      DO 10 KH=1,NEL
10    WRITE (NOUT,1003) C(KH),(MW(KH,NH),NH=1,ND)
1003  FORMAT(' MATRIXELEMENT ',E12.4,' M-WERTE',19I3)
C     PREPARATION
20    JJ=0
      NQ=N5+(ND-(2*NKAPO+1)*N5)
      ITEN=0
      IF(NQ.NE.N5) ITEN=1
C      TREATS TENSOR SEPARATELY
      IH=0
      DO 22 MH=1,NQ
      DO 22 NH=MH,NQ
      IH=IH+1
22    MKOM(MH,NH)=IH
C     MKOM GIVES THE COMBINATION OF SIGMA FACTORS
      IH=0
      DO 24  MH=1,ND
      LQ(MH)=L(MH)
24    IH=IH+L(MH)
      IF(MOD(IH,2).NE.0) GOTO 2000
C     CHECK IF SUM L EVEN
      IQ=IH/2
      IF(IH.GT.0) GOTO 28
      C1=0.
      DO 26 NH=1,NEL
26    C1=C1+C(NH)
      CW(1)=C1
      JJ=1
      MVM(1,1)=0
C     NO SIGMA FACTORS
      RETURN
28    DO 32 NH=1,NEL
      IH=0
      DO 30 MH=1,ND
30    IH=IH+MW(NH,MH)
      IF(IH.NE.0) GOTO 2001
C     CHECK IF SUM M-VALUES ZERO
32    CONTINUE
C     LOOK FOR POSSIBLE SIGMA FACTORS
C     DETERMINE INDICES JVN AND JVM
      IH=1
      MD=ND-1
      M=1
34    JVM(IH)=M
      IF(LQ(M)) 2002,54,36
36    LQ(M)=LQ(M)-1
      IF(IH-1) 2003,40,38
38    IF(M.NE.JVM(IH-1)) GOTO 40
      N=JVN(IH-1)
      GOTO 42
40    N=M+1
42    JVN(IH)=N
      IF(LQ(N)) 2002,48,44
44    LQ(N)=LQ(N)-1
      IF(IH-IQ) 46,66,2003
46    IH=IH+1
C     NEXT FACTOR SIGMA
      GOTO 34
C     LOOK FOR OTHER N
48    IF(N-ND) 50,52,2003
50    N=N+1
      GOTO 42
52    LQ(M)=LQ(M)+1
C     RESTORE L-VALUE AND LOOK FOR OTHER M
54    IF(M-MD)  56,58,2003
56    M=M+1
      GOTO 34
58    IF(IH-1)  2003,999,60
60    IH=IH-1
62    M=JVM(IH)
C     ENTRY POINT FOR SEARCH OF OTHER SIGMA COMBINATION
      N=JVN(IH)
      LQ(N)=LQ(N)+1
      GOTO 48
C     JVM AND JVN DETERMINED
C     DETERMINE INDICES IG
66    C1=0.
C     LOOP OVER ALL M-COMBINATIONS
      DO 500 KE=1,NEL
      IF(KAUS.GT.3) WRITE(NOUT,1500) KE
1500   FORMAT(' LOOP UEBER M-WERTE',I5,' TES M-ELEMENT')
      C2=C(KE)
      DO 70 NH=1,ND
      KY(NH)=L(NH)-1+MW(KE,NH)
70    KYS(NH)=L(NH)-1-MW(KE,NH)
C     KY AND KYS GIVE THE NUMBER OF ARROWS FROM N AND TO N
C     CHOOSE SIGMA FACTOR
C     PREPARATION
      IH=1
      M=JVM(IH)
      N=JVN(IH)
      GOTO 100
80    IF(IH-IQ)  82,200,2003
82    IH=IH+1
      M=JVM(IH)
      N=JVN(IH)
C     CHECK IF SAME SIGMA FACTOR
      IF(M-JVM(IH-1).NE.0) GOTO 100
      IF(N-JVN(IH-1).NE.0) GOTO 100
      IF(IG(IH-1))  140,160,100
C     TRY FOR IG=+1,ARROW FROM M TO N
100   IF(KY(N).LE.0) GOTO 140
      IF(KYS(M).LE.0) GOTO 140
      IG(IH)=1
      KY(N)=KY(N)-2
      KYS(M)=KYS(M)-2
C     NEXT SIGMA FACTOR
      GOTO 80
C     TRY FOR IG=-1,ARROW FROM N TO M
140   IF(KYS(N)) 190,162,142
142   IF(KY(M)) 190,164,144
144   IG(IH)=-1
      KY(M)=KY(M)-2
      KYS(N)=KYS(N)-2
C     NEXT SIGMA FACTOR
      GOTO 80
C     TRY FOR IG=0, NO ARROW
160   IF(KYS(N).LT.0) GOTO 190
162   IF(KY(M).LT.0) GOTO 190
164   IF(KY(N).LT.0) GOTO 190
      IF(KYS(M).LT.0) GOTO 190
      IG(IH)=0
      KY(N)=KY(N)-1
      KYS(N)=KYS(N)-1
      KY(M)=KY(M)-1
      KYS(M)=KYS(M)-1
C     NEXT SIGMA FACTOR
      GOTO 80
190   IF(IH-1) 2003,500,192
192   IH=IH-1
194   N=JVN(IH)
      M=JVM(IH)
      IF(IG(IH))  196,197,195
C     RESTORE IG=+1
195   KY(N)=KY(N)+2
      KYS(M)=KYS(M)+2
C     TRY IG=-1
      GOTO 140
C     RESTORE IG=-1
196   KY(M)=KY(M)+2
      KYS(N)=KYS(N)+2
C     TRY IG=0
      GOTO 164
C     RESTORE IG=0
197   KY(M)=KY(M)+1
      KYS(M)=KYS(M)+1
      KY(N)=KY(N)+1
      KYS(N)=KYS(N)+1
C     TRY PRIOR SIGMA FACTOR
      GOTO 190
C     END OF IG
C     UTILYSE EQUALITY OF SIGMA FACTORS
200    DO 250 IK=1,IQ
C     VM FACTOR
      KH=(JVM(IK)-1)/N5
      IVM(IK)=JVM(IK)-KH*N5*NKAPO
C      NO POLYNOMIAL = NKAPO =0,NO REDUCTION NECESSARY
      IF(ITEN*KH.EQ.3) IVM(IK)=IVM(IK)+N5
      KH=MOD(KH,3)
      IWM(IK)=KH+1
C     VN FACTORS
      KH=(JVN(IK)-1)/N5
      IVN(IK)=JVN(IK)-KH*N5*NKAPO
C      NO POLYNOMIAL = NKAPO =0,NO REDUCTION NECESSARY
      IF(ITEN*KH.EQ.3) IVN(IK)=IVN(IK)+N5
      KH=MOD(KH,3)
      IWN(IK)=KH+1
      IF(KAUS.GT.3)  WRITE (NOUT,1200) IK,JVM(IK),JVN(IK)
     1  ,IVM(IK),IVN(IK),IWM(IK),IWN(IK),IG(IK)
1200  FORMAT(' LOOP W-INDICES, IK,JVM,JVN,IVM,IVN,IWM,IWN,IG =',8I5)
C     ORDER VM .LE. VN
      IF(IVM(IK).LE.IVN(IK)) GOTO 230
      KH=IVM(IK)
      IVM(IK)=IVN(IK)
      IVN(IK)=KH
C     ORDER WM .LE. WN
230   IF(IWM(IK).LE.IWN(IK)) GOTO 250
      KH=IWM(IK)
      IWM(IK)=IWN(IK)
      IWN(IK)=KH
250   CONTINUE
C     DETERMINE CONSTANT
      C3=C2
      NF=2
      DO 320 IK=1,IQ
      IS=IK-1
      IF(IG(IK).NE.0) C3=C3*(-.5)
      IF(IS.EQ.0) GOTO 320
C     SKIP FIRST FACTOR
C     CHECK FOR EQUALITY OF ALL INDICES AND CALCULATE FACULTY
      IF(JVN(IK).NE.JVN(IS)) GOTO 310
      IF(JVM(IK).NE.JVM(IS)) GOTO 310
      IF(IG(IK).NE.IG(IS)) GOTO 310
      NF=NF+1
      GOTO 320
310   C3=C3*F8(NF)
      NF=2
320   CONTINUE
      C3=C3*F8(NF)
      C1=C1+C3
C     FACTOR DETERMINED AND SUMMED UP
      IH=IQ
C     RESTORE KY AND KYS VALUES
      IF(KAUS.GT.1) WRITE (NOUT,1320) C1,IH,(IVM(LX),IVN(LX),LX=1,IH)
1320   FORMAT(' NEUES ELEMENT C1=',G15.5,' IH=',I5,' VM VN',20I3)
      GOTO 194
500   CONTINUE
C     CHECK IF FACTOR DIFFERENT FROM ZERO
      IF(ABS(C1).LT.1.E-10) GOTO 610
      JJ=JJ+1
      IF(JJ.GT.NDIMCW) GOTO 2004
C      ORDER INDICES
      CALL SORT2(IQ,IVM,IVN)
      CW(JJ)=C1
      DO 600 IK=1,IQ
      MH=IVM(IK)
      NH=IVN(IK)
600   MVM(IK,JJ)=MKOM(MH,NH)
C     LOOK FOR DIFFERENT SIGMA FACTORS
      IF(KAUS.GT.1) WRITE(NOUT,1610)JJ,CW(JJ),(MVM(IK,JJ),IK=1,IQ)
1610  FORMAT(I5,'TES M-ELEMENT C =',G15.5,' MVM',19I3)
610   IH=IQ
      GOTO 62
999   CONTINUE
      RETURN
2000  WRITE(NOUT,1004)
1004  FORMAT(' +++ SUMME L-WERTE UNGERADE ++')
      STOP 21
2001  WRITE (NOUT,1005) NH,(MW(NH,IH),IH=1,ND)
1005  FORMAT(' ++ SUMME M-WERTE NICHT NULL, NH=',I4,' M-WERTE',19I3)
      STOP 22
2002  WRITE (NOUT,1006) M,LQ(M)
1006  FORMAT(' DER',I5,'TE L-WERT IST NEGATIV',I5)
      STOP 23
2003  WRITE (NOUT,1007) IH,M,N
1007  FORMAT(' INDICES AUSSERHALB BEREICH,IH,M,N',3I5)
      STOP 24
2004  WRITE (NOUT,1008) JJ
1008  FORMAT(' DIMENSIONIERUNG VON CW ZU KLEIN JJ=',2I5)
      STOP 25
      END
      SUBROUTINE PORD(LL,LR,MZAHL,KZL,KZR)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C     THIS ROUTINE PUTS ALL MATRIXELEMENTS INTO A SCHEME,ACCORDING TO
C     THE OCCURRING SIGMA FACTORS
      INCLUDE 'par/luelma'
C
      COMMON /SDA/CW(NDIMCW), MVM(NZIQMA,NDIMCW), KAUS,
     *            N5, IQ, NKAPO
      COMMON/ORD/EPO(NDIMOR),KZAHL, JQ(NDIMOR),
     *           INDPO(NDIMOR), MVK(NZIQMA,NDIMOR),
     *           INDEX(NDIMOR)
      IND=(((KZL-1)*10+KZR-1)*100+LL-1)*100+LR-1
      DO 200 MH=1,MZAHL
C     LOOP OVER ALL ELEMENTS DETERMINED IN SDA1
      IF(KZAHL.LE.0) GOTO 100
      DO 40 KH=1,KZAHL
      IF(IQ.NE.JQ(KH))  GOTO 40
C     CHECK IF COINCIDENCE WITH PREVIOUS ELEMENTS
      IF(IND.NE.INDPO(KH)) GOTO 40
C     CHECK IF SAME ANGULAR MOMENTUM STRUCTURE
      IF(IQ.EQ.0)  GOTO 50
      JH=0
      DO 10 IH=1,IQ
10    JH=JH+IABS(MVK(IH,KH)-MVM(IH,MH))
      IF(JH.EQ.0)  GOTO 50
40    CONTINUE
      GOTO 100
50    EPO(KH)=EPO(KH)+CW(MH)
C     MATRIX ELEMENT EXISTS ALREADY, ADD UP
      GOTO 200
100   KZAHL=KZAHL+1
C     NEW MATRIX ELEMENT
      IF(KZAHL.GT.NDIMOR)  GOTO 1001
      EPO(KZAHL)=CW(MH)
      INDPO(KZAHL)=IND
      JQ(KZAHL)=IQ
      IF(IQ.LE.0)  GOTO 190
C     STORE SIGMA INDICES FOR LATER COMPARISON
      DO 120  IH=1,IQ
120   MVK(IH,KZAHL)=MVM(IH,MH)
      GOTO 200
190   MVK(1,KZAHL)=0
200   CONTINUE
      RETURN
1001  WRITE (NOUT,1002) KZAHL,MH
1002  FORMAT(' DIMENSION IN PORD ZU KLEIN,KZAHL=',2I5,'TES ELEMENT')
      STOP 41
      END
      SUBROUTINE SORT2(IQ,JVM,JVN)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C     SORT2 ORDER THE INDICES JVM,JVN LEXICOGRAPHICALLY
      INCLUDE 'par/luelma'
C
      DIMENSION IND(NZIQMA),JVM(NZIQMA),JVN(NZIQMA)
      IF(IQ.LE.1) GOTO 200
      IQ1=IQ-1
C     ORDER ACCORDING TO JVM
      CALL SORTA(1,IQ,JVM,IND)
      CALL PUTA(1,IQ,JVN,IND)
C      ORDER JVN
      JMIN=0
      DO 30 I=1,IQ1
      IF(JVM(I+1).NE.JVM(I)) GOTO 25
      IF(JMIN.NE.0) GOTO 22
      JMIN=I
22    IF(I.NE.IQ1) GOTO 30
      JMAX=I+1
      GOTO 27
25    IF(JMIN.LE.0) GOTO 30
      JMAX=I
27    IF(JMAX.EQ.JMIN) GOTO 30
      CALL SORTA(JMIN,JMAX,JVN,IND)
      JMIN=0
30    CONTINUE
200   RETURN
      END
      SUBROUTINE SORTA(IMIN,IMAX,JQ,IND)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C     SORTA ORDERS ARRAY J AND GIVES THE ORDERING IN IND
      INCLUDE 'par/luelma'
C
      DIMENSION JQ(NZIQMA),IND(NZIQMA)
      DO 10 I=IMIN,IMAX
10    IND(I)=I
      IF(IMAX.EQ.IMIN) GOTO 200
      IQ1=IMAX-1
      DO 20 I=IMIN,IQ1
      J=I
15    IF(JQ(J+1).GE.JQ(J)) GOTO 20
      JM=JQ(J)
      ISUB=IND(J)
      JQ(J)=JQ(J+1)
      IND(J)=IND(J+1)
      JQ(J+1)=JM
      IND(J+1)=ISUB
      J=J-1
      IF(J.GE.IMIN) GOTO 15
20    CONTINUE
200   RETURN
      END
      SUBROUTINE PUTA(IMIN,IMAX,JQ,IND)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C     PUTA INTERCHANGES ELEMENTS IN JQ ACCORDING TO IND
      INCLUDE 'par/luelma'
C
      DIMENSION JQ(NZIQMA),IND(NZIQMA),KH(NZIQMA)
      DO 10 I=IMIN,IMAX
10    KH(I)=JQ(I)
      DO 20 I=IMIN,IMAX
      IH=IND(I)
20    JQ(I)=KH(IH)
      RETURN
      END
      SUBROUTINE WRITAP(LAUS,IQM,NBAND)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INCLUDE 'par/luelma'
C
      COMMON/ORD/EPO(NDIMOR),KZAHL, JQ(NDIMOR),
     *           INDPO(NDIMOR), MVK(NZIQMA,NDIMOR),
     *           INDEX(NDIMOR)
      DO 40 KH=1,KZAHL
      IF(LAUS.LE.0) GOTO 40
      II=JQ(KH)
      WRITE(NOUT,1000) II,INDPO(KH),EPO(KH),(MVK(IH,KH),IH=1,II)
1000  FORMAT(I5,' SIGMA-FAKTOREN U INDEX',I8,' EPO=',G15.5,
     1 'MVK',19I3)
40    CONTINUE
      WRITE (NBAND)(JQ(KH),INDPO(KH),EPO(KH),(MVK(IH,KH),IH=1,IQM),
     1 KH=1,KZAHL)
      DO 50 IH=1,IQM
      DO 45 KH=1,KZAHL
45    MVK(IH,KH)=0.
50    CONTINUE
      RETURN
      END
      DOUBLE PRECISION FUNCTION CLG(J1,J2,J3,M1,M2)
C
C     CLG BERECHNET DIE CLEBSCH-GORDAN-KOEFFIZIENTEN
C     (J1/2,M1/2;J2/2,M2/2|J3/2,(M1+M2)/2) NACH
C     EDMONDS 'ANGULAR MOMENTUM IN QUANTUM MECHANICS',
C     PRINCETON, 1960 GLEICHUNGEN (3.10.60), (3.7.3)
C     UND TABELLE 2 (1. GLEICHUNG)
C
C     BENUTZT COMMON /COMY/ MIT DEN LOGRITHMEN DER
C     FAKULTAETEN
C
C     M. UNKELBACH 1989
C     LETZTE AENDERUNG: 06.02.89
C
C
      INTEGER JW1, JW2, JW3, MW1, MW2, MW3, JSUM, JSUM1,
     *        JDIF1, JDIF2, JDIF3, JMSUM1, JMSUM2, JMSUM3,
     *        JMDIF1, JMDIF2, JMDIF3, JJM1, JJM2, IMAX, IMIN,
     *        I, J1, J2, J3, M1, M2
C
      DOUBLE PRECISION FAKLN, CLGH
C
      COMMON /COMY/ FAKLN(0:99)
C     FAKLN(I) = LOG(I!)
C
C
C
C
      JW1=J1
      JW2=J2
      JW3=J3
      MW1=M1
      MW2=M2
C
C     CHECK, OB CLG = 0
      CLG=0.
      IF (JW1.LT.IABS(MW1)) RETURN
      IF (JW2.LT.IABS(MW2)) RETURN
      IF (JW3.GT.JW1+JW2.OR.JW3.LT.IABS(JW1-JW2)) RETURN
      MW3=MW1+MW2
      IF (JW3.LT.IABS(MW3)) RETURN
      JMSUM1=JW1+MW1
      JMSUM2=JW2+MW2
      JMSUM3=JW3+MW3
      IF (MOD(JMSUM1,2).EQ.1) RETURN
      IF (MOD(JMSUM2,2).EQ.1) RETURN
      IF (MOD(JMSUM3,2).EQ.1) RETURN
C
C
      JSUM=(JW1+JW2+JW3)/2
      JSUM1=JSUM+1
      JDIF1=JSUM-JW1
      JDIF2=JSUM-JW2
      JDIF3=JSUM-JW3
C
      IF (IABS(MW1)+IABS(MW2).EQ.0) GOTO 100
C
C     NORMALE CLEBSCH-GORDAN-KOEFFIZIENTEN
      JMSUM1=JMSUM1/2
      JMDIF1=JMSUM1-MW1
      JMSUM2=JMSUM2/2
      JMDIF2=JMSUM2-MW2
      JMSUM3=JMSUM3/2
      JMDIF3=JMSUM3-MW3
      JJM1=JDIF1+JMDIF1
      JJM2=JDIF3-JMDIF1
      IMIN=MAX0(0,-JJM2)
      IMAX=MIN0(JMDIF1,JMDIF3)
C
      CLGH=0.
      DO 50, I=IMIN, IMAX
       CLGH=CLGH+DBLE(1-2*MOD(I,2))*
     *     EXP(FAKLN(JMSUM1+I)+FAKLN(JJM1-I)-FAKLN(I)-FAKLN(JMDIF1-I)-
     *         FAKLN(JMDIF3-I)-FAKLN(JJM2+I))
50    CONTINUE
C
      IF (IMIN.GT.IMAX) CLGH=1.
      CLGH=CLGH*EXP((FAKLN(JDIF3)+FAKLN(JMDIF1)+FAKLN(JMDIF2)+
     *             FAKLN(JMDIF3)+FAKLN(JMSUM3)-FAKLN(JSUM1)-
     *             FAKLN(JDIF1)-FAKLN(JDIF2)-FAKLN(JMSUM1)-
     *             FAKLN(JMSUM2)+FAKLN(JW3+1)-FAKLN(JW3))*.5)
      CLG=CLGH*DBLE(1-2*MOD(JMDIF1,2))
C
C     ENDE DER BERECHNUNG FUER NORMALE CLEBSCH-GORDAN-KOEFFIZIENTEN
      RETURN
C
C
C
100   CONTINUE
C     PARITAETSCLEBSCH
C
      IF (MOD(JSUM,2).EQ.1) RETURN
C
      CLGH=EXP((FAKLN(JDIF1)+FAKLN(JDIF2)+FAKLN(JDIF3)-FAKLN(JSUM1)+
     *         FAKLN(JW3+1)-FAKLN(JW3))*.5+
     *        FAKLN(JSUM/2)-FAKLN(JDIF1/2)-FAKLN(JDIF2/2)-
     *        FAKLN(JDIF3/2))
      CLG=CLGH*DBLE(1-2*MOD((JSUM+JW1-JW2)/2,2))
C
C
C     ENDE DER RECHNUNG FUER PARITAETSCLEBSCH
      RETURN
      END
      FUNCTION F6J(JD1,JD2,JD3,LD1,LD2,LD3)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C     VERSION I F6J FUNCTION CALLS S6J  FORTRAN IV
C     VEREINFACHT 27.9.95 H.M.H
      J1=JD1
      J2=JD2
      J3=JD3
      L1=LD1
      L2=LD2
      L3=LD3
C     ANGULAR MOMENTUM COUPLING TESTS FOR 6J COEFFICIENT
      F6J=0.0
      IF(J1.LT.0 .OR. J2.LT.0 .OR. J3.LT.0) RETURN
      IF(L1.LT.0 .OR. L2.LT.0 .OR. L3.LT.0) RETURN
      IF(MOD(J1+J2+J3,2).NE.0) RETURN
      IF(J3.GT.J1+J2 .OR. J3.LT.ABS(J1-J2)) RETURN
      IF(MOD(J1+L2+L3,2).NE.0) RETURN
      IF(L3.GT.J1+L2 .OR. L3.LT.ABS(J1-L2)) RETURN
      IF(MOD(L1+J2+L3,2).NE.0) RETURN
      IF(L3.GT.L1+J2 .OR. L3.LT.ABS(L1-J2)) RETURN
      IF(MOD(L1+L2+J3,2).NE.0) RETURN
      IF(J3.GT.L1+L2 .OR. J3.LT.ABS(L1-L2)) RETURN
      F6J=S6J(J1,J2,J3,L1,L2,L3)
      RETURN
      END
      FUNCTION F9J(JD1,JD2,JD3,JD4,JD5,JD6,JD7,JD8,JD9)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C     F9J VERSION I  CALLS S6J  FORTRAN IV
C     VEREINFACHT 27.9.95 H.M.H.
      DIMENSION KN(6),KX(6),NN(6)
      J1=JD1   
      J2=JD2  
      J3=JD3 
      J4=JD4
      J5=JD5
      J6=JD6
      J7=JD7
      J8=JD8
      J9=JD9
      F9J= 0.0
C     ANGULAR MOMENTUM COUPLING TESTS FOR 9J COEFFICIENT 
      IF(MOD(J1+J2+J3,2).NE.0) RETURN
      IF(J3.GT.J1+J2 .OR. J3.LT.ABS(J1-J2)) RETURN
      IF(MOD(J4+J5+J6,2).NE.0) RETURN
      IF(J6.GT.J4+J5 .OR. J6.LT.ABS(J4-J5)) RETURN
      IF(MOD(J7+J8+J9,2).NE.0) RETURN
      IF(J9.GT.J7+J8 .OR. J9.LT.ABS(J7-J8)) RETURN
      IF(MOD(J1+J4+J7,2).NE.0) RETURN
      IF(J7.GT.J1+J4 .OR. J7.LT.ABS(J1-J4)) RETURN
      IF(MOD(J2+J5+J8,2).NE.0) RETURN
      IF(J8.GT.J2+J5 .OR. J8.LT.ABS(J2-J5)) RETURN
      IF(MOD(J3+J6+J9,2).NE.0) RETURN
      IF(J9.GT.J3+J6 .OR. J9.LT.ABS(J3-J6)) RETURN
      KN(1)=MAX0(IABS(J2-J6),IABS(J1-J9),IABS(J4-J8))
      KN(2)=MAX0(IABS(J2-J7),IABS(J5-J9),IABS(J4-J3))
      KN(3)=MAX0(IABS(J6-J7),IABS(J5-J1),IABS(J8-J3))
      KN(4)=MAX0(IABS(J6-J1),IABS(J2-J9),IABS(J5-J7))
      KN(5)=MAX0(IABS(J2-J4),IABS(J3-J7),IABS(J6-J8))
      KN(6)=MAX0(IABS(J3-J5),IABS(J1-J8),IABS(J4-J9))
      KX(1)=MIN0(J2+J6,J1+J9,J4+J8)
      KX(2)=MIN0(J2+J7,J5+J9,J4+J3)
      KX(3)=MIN0(J6+J7,J5+J1,J8+J3)
      KX(4)=MIN0(J1+J6,J2+J9,J5+J7)
      KX(5)=MIN0(J2+J4,J3+J7,J6+J8)
      KX(6)=MIN0(J3+J5,J1+J8,J4+J9)
      DO 35 K=1,6
   35 NN(K)=KX(K)-KN(K)
      KSIGN=1
      I=MIN0(NN(1),NN(2),NN(3),NN(4),NN(5),NN(6))
      DO 40 K=1,6
      IF(I-NN(K))40,50,40
   40 CONTINUE
   50 KMIN=KN(K)+1
      KMAX=KX(K)+1
      GO TO(130,52,53,54,55,56),K
   52 J=J1
      J1=J5
      J5=J
      J=J3
      J3=J8
      J8=J
      J=J6
      J6=J7
      J7=J
      GO TO 130
   53 J=J2
      J2=J7
      J7=J
      J=J3
      J3=J4
      J4=J
      J=J5
      J5=J9
      J9=J
      GO TO 130
   54 J=J1
      J1=J2
      J2=J
      J=J4
      J4=J5
      J5=J
      J=J7
      J7=J8
      J8=J
      GO TO 120
   55 J=J1
      J1=J3
      J3=J
      J=J4
      J4=J6
      J6=J
      J=J7
      J7=J9
      J9=J
      GO TO 120
   56 J=J2
      J2=J3
      J3=J
      J=J5
      J5=J6
      J6=J
      J=J8
      J8=J9
      J9=J
  120 KSIGN=(1-MOD(J1+J2+J3+J4+J5+J6+J7+J8+J9,4))
C     SUMMATION OF SERIES OF EQUATION (2)  
  130 SUM=0.0                             
      SIG=(-1)**(KMIN-1)*KSIGN
      FLK=KMIN                           
      DO 200 K=KMIN,KMAX,2              
      TERM=FLK*S6J(J1,J4,J7,J8,J9,K-1)*S6J(J2,J5,J8,J4,K-1,J6)
     1*S6J(J3,J6,J9,K-1,J1,J2)
      FLK=FLK+2.0                      
  200 SUM=SUM+TERM                    
      F9J=SUM*SIG
      RETURN                       
      END                         
      FUNCTION S6J(JD1,JD2,JD3,LD1,LD2,LD3)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C     VERSION I  FORTRAN IV
      DIMENSION MA(4),MB(3),MED(12)
      COMMON /FACT/FL(322),NCALL
      DATA NCALL/1/
      J1=JD1
      J2=JD2
      J3=JD3
      L1=LD1
      L2=LD2
      L3=LD3
C     DETERMINE WHETHER TO CALCULATE FL(N) S
      IF(NCALL.EQ.0) GOTO 15
      NCALL=0
C     CALCULATE FL(N) S
      FL(1)=0.0
      FL(2)=0.0
      DO 50 N= 3,322
      FN=N-1
   50 FL(N)=FL(N-1)+LOG(FN)
   15 MED(1)=(-J1+J2+J3)/2
      MED(2)=(+J1-J2+J3)/2
      MED(3)=(+J1+J2-J3)/2
      MED(4)=(-J1+L2+L3)/2
      MED(5)=(+J1-L2+L3)/2
      MED(6)=(+J1+L2-L3)/2
      MED(7)=(-L1+J2+L3)/2
      MED(8)=(+L1-J2+L3)/2
      MED(9)=(+L1+J2-L3)/2
      MED(10)=(-L1+L2+J3)/2
      MED(11)=(+L1-L2+J3)/2
      MED(12)=(+L1+L2-J3)/2
      MA(1)=MED(1)+MED(2)+MED(3)
      MA(2)=MED(4)+MED(5)+MED(6)
      MA(3)=MED(7)+MED(8)+MED(9)
      MA(4)=MED(10)+MED(11)+MED(12)
      MB(1)=MA(1)+MED(12)
      MB(2)=MA(1)+MED(4)
      MB(3)=MA(1)+MED(8)
C     DETERMINE MAXIMUM OF (J1+J2+J3),(J1+L2+L3),(L1+J2+L3),(L1+L2+J3)
      MAX=MA(1)
      DO 30 N=2,4
      IF (MAX-MA(N)) 20,30,30
   20 MAX=MA(N)
   30 CONTINUE
C     DETERMINE MINIMUM OF (J1+J2+L1+L2), (J2+J3+L2+L3),(J3+J1+L3+L1)
      MIN=MB(1)
      DO 51 N=2,3
      IF (MIN-MB(N)) 51,51,40
   40 MIN=MB(N)
   51 CONTINUE
      MINH=MIN
      KMAX=MIN-MAX
      MINP1=MIN+1
      MINI  =MINP1-MA(1)
      MIN2=MINP1-MA(2)
      MIN3=MINP1-MA(3)
      MIN4=MINP1-MA(4)
      MIN5=MINP1+1
      MIN6=MB(1)-MIN
      MIN7=MB(2)-MIN
      MIN8=MB(3)-MIN
      UK=1.E-15
      S=1.0E-15
      IF (KMAX) 65,65,55
   55 DO 60 K=1,KMAX
      UK=-UK*DBLE(MINI-K)*DBLE(MIN2-K)*DBLE(MIN3-K)*DBLE(MIN4-K)/
     1 (DBLE(MIN5-K)*DBLE(MIN6+K)*DBLE(MIN7+K)*DBLE(MIN8+K))
C     CUT OFF SERIES AT 1.0D-25
      IF(ABS(UK)-1.E-25) 65,60,60
   60 S=S+UK
   65 S=S*1.0E+15
C     CALCULATE DELTA FUNCTIONS
      DELOG=0.0
      DO 70 N=1,12
      NUM=MED(N)
   70 DELOG=DELOG+FL(NUM+1)
      NUM1=MA(1)+2
      NUM2=MA(2)+2
      NUM3=MA(3)+2
      NUM4=MA(4)+2
      DELOG=DELOG-FL(NUM1)-FL(NUM2)-FL(NUM3)-FL(NUM4)
      DELOG=0.5*DELOG
      ULOG=FL(MIN5)-FL(MINI)-FL(MIN2)-FL(MIN3)-FL(MIN4)-FL(MIN6+1)-
     1   FL(MIN7+1)-FL(MIN8+1)
      PLOG=DELOG+ULOG
      P=EXP (PLOG)
      S6J=P*S
      IF(MOD(MINH,2).NE.0)  S6J=-S6J
      RETURN
      END
