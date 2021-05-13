      PROGRAM LUELMA (INPUT,OUTPUT,TAPE5=INPUT,TAPE6=OUTPUT,TAPE9)
C     PROGRAM LUELMA PREPARES THE CALCULATION OF THE REDUCED
C     SPACE MATRIX ELEMENTS
C
C   FUER ELEKTROMAGNETISCHE UEBERGAENGE
C
C
C                      T.M. 1985
C
C     LETZTE AENDERUNG H.M.H. 26.2.88 IQM
C
C     MAXIMAL 4 CLUSTERS CAN BE HANDLED ..NZC
C     MAXIMAL 5 ZERLEGUNGEN             ..NZF
C     MAXIMAL 10 ANGULAR MOMENTUM STRUCTURES FOR EACH ZERLEGUNG ..NZLW
C     MAXIMAL ANGULAR MOMENTUM = 10
C     MAXIMAL POLYNOMIAL =2*3
C
      COMMON /COMY/D(100)
C     LOGARITHMS OF FACULTIES FOR ANGULAR MOMENTUM ROUTINES
      COMMON /SDA/F8(21),CW(191),MVM(19,191),KAUS,N5,NDIMCW,IQ
     1 ,NKAPO
      COMMON /MVA/ KV(6),KMV(6,81),NDIMK,MMOEG,MKAUS
      COMMON /HIVA/ NOUT,NBAND,NDIMPK
      COMMON /ORD/ KZAHL, JQ(635),EPO(635),INDPO(635),MVK(19,635)
      DIMENSION F3(21),F7(21),IENT(5),NZLW(5),NZC(5),NZPO(5)
      DIMENSION LW(5,10,5),KP(3,10,5),F4(4),F5(4)
      DIMENSION LV(19),MV(19),LKOP1(2,10,5),LKOP2(2,10,5)
      DIMENSION C(123),MW(123,19),LVMOD(19),KWL(5,5),KWR(5,5)
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
      NZFMAX=5
      NDIMC=123
      NDIMCW=191
      NDIMK=81
      NDIMPK=635
      NKAPO=0
      MMOEG=1
C     NO POLYNOMIALS FOR NKAPO=0
C
C
      D(1)=0.
      D(2)=0.
      DO 100 I=2,99
100   D(I+1)=ALOG(FLOAT(I))+D(I)
      INPUT = 5
      NOUT =6
      FT=SQRT(1./(4.*3.1415926536))
      F7(1)=FT
      F3(1)=1.
      F8(1)=1.
      FH=1.
      DO 102 NH=2,21
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
      READ (INPUT,1000) (IENT(NOP),NOP=1,4)
C     IENT =1 DETERMINES OPERATORS TO BE CALCULATED
C     IENT(1..4) = EL-BAHN-,EL-SPIN-,ML-BAHN-,ML-SPIN-OPERATOR
      READ (INPUT,1000) NZF,MUL
C     NZF NUMBER OF ZERLEGUNGEN, MUL MULTIPOLARITAET
      DO 110 MH=1,NZF
110   READ (INPUT,1000) NZLW(MH),NZC(MH),NZPO(MH)
C
      REWIND NBAND
      WRITE(NBAND)NZF,MUL,(NZLW(MH),NZC(MH),MH=1,NZF),(IENT(MH),MH=1,4)
     1 ,(NZPO(MH),MH=1,NZF)
C
C     READ AND WRITE ACTUAL ANGULAR MOMENTA AND POLYNOMIALS
C
      WRITE(NOUT,1060)
1060  FORMAT(1H1)
      WRITE (NOUT,1001) NZF
1001  FORMAT(" ZAHL DER ZERLEGUNGEN",I3)
      IF(NZF.GT.NZFMAX) STOP 1
      DO 140 MH=1,NZF
      WRITE (NOUT,1002) MH
1002  FORMAT(//,"0 ZERLEGUNG",I5)
      N1=NZC(MH)-1
      IF(N1.GT.3) STOP 2
      N3=N1-1
C     NUMBER OF ANGULAR MOMENTA
      N2=NZLW(MH)
      WRITE (NOUT,1003) NZC(MH),NZLW(MH)
1003  FORMAT(" ANZAHL DER CLUSTER",I5,
     1 "      ANZAHL DER DREHIMPULSSTRUKTUREN",I5)
      DO 130 L=1,N2
      WRITE (NOUT,1004) L
1004  FORMAT(/," DREHIMPULSSTRUKTUR",I5)
      READ (INPUT,1000)(LW(M,L,MH),M=1,N1)
C    RELATIV ANGULAR MOMENTA
C
      WRITE (NOUT,1005) (M,LW(M,L,MH),M=1,N1)
1005  FORMAT(I4," TER DREHIMPULS =",I3," RELATIVDREHIMPULS")
      IF(N3.LE.0) GOTO 130
      DO 120 M=1,N3
      MM=N1+M
      READ(INPUT,1000) LW(MM,L,MH),LKOP1(M,L,MH),LKOP2(M,L,MH)
C     LW COUPLED ANGULAR MOMENTUM FROM LKOP1 AND LKOP2
      IF(N3.GT.1) GOTO 116
      LKOP1(1,L,MH) =1
      LKOP2(1,L,MH) =2
116   WRITE (NOUT,1006) MM,LW(MM,L,MH),LKOP1(M,L,MH),LKOP2(M,L,MH)
1006  FORMAT(I4," TER DREHIMPULS =",I3," GEKOPPELT AUS DREHIMPULS",
     1 I3," UND DREHIMPULS",I3)
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
1007  FORMAT(" BEI DER",I3," TEN POLYNOM FUNKTION SIND POLYNOME DER"
     1 ," ORDNUNG",3I3)
      N4=N1+N3
140   WRITE (NBAND)((LW(M,L,MH),M=1,N4),L=1,N2),
     1 ((KP(M,KH,MH),M=1,N1),KH=1,MPO)
C     END OF INPUT
C     START OF CALCULATION OF MATRIX ELEMENTS
      DO 900 MFL=1,NZF
C     LOOP ZERLEGUNGEN LEFT
      N1=NZC(MFL)-1
      N8=N1+N1-1
      N10=N1+1
      N2=NZLW(MFL)
      K1=NZPO(MFL)
      DO 890 MFR=1,MFL
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
1009  FORMAT(/," DER",I3," TE DREHIMPULS WIRD ERNIEDRIGT")
164   DO 800 LL=1,N2
C      LOOP DREHIMPULSSTRUKTUREN LEFT
      NPWL=0
      DO 166 KH=1,N1
166   NPWL=NPWL+LW(KH,LL,MFL)
      NPWL=(-1)**NPWL
      DO 790 LR=1,N4
C     LOOP DREHIMPULSSTRUKTUREN RIGHT
      NPWR=0
      DO 168 KH=1,N3
168   NPWR=NPWR+LW(KH,LR,MFR)
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
      BR=LW(N9,LR,MFR)
      OZ=IDL
      F=1.
      X=FLOAT(MUL)
      IF(NOP.EQ.4) X=X-1.
      GOTO (182,190,182,190),NOP
182    IZUS1=1
      IZUS2=-1
190   Y=YG(BR,X,BL,BR,OZ)
C     MAXIMAL COUPLING
      IF(Y.NE.0.) GOTO 192
      F=0.
      IF(MAUS.LE.0) GOTO 790
      WRITE (NOUT,1010) LW(N9,LR,MFR),MUL,LW(N8,LL,MFL)
1010  FORMAT(" KOPPLUNG LIEFERT NULL",3I5)
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
1011  FORMAT(/," MATRIXELEMENT ZWISCHEN ",
     1 "DEM",I3," TEN POLYNOM DER",I3,
     2 " TEN DREHIMPULSSTRUKTUR DER",I3,
     3 " TEN ZERLEGUNG MIT PARITAET",I3," UND",/,24X,
     4 "DEM",I3," TEN POLYNOM DER",I3,
     5 " TEN DREHIMPULSSTRUKTUR DER",I3,
     6 " TEN ZERLEGUNG MIT PARITAET",I3)
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
      AH1=-KWL(MM,2)
      AH2=-KWL(MM,3)
      F1=F1*F3(M1)*F3(M2)*YG(BH1,BH2,BH3,AH1,AH2)
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
      AH1=KWR(NN,2)
      AH2=KWR(NN,3)
282   F2=F2*F3(M1)*F3(M2)*YG(BH1,BH2,BH3,AH1,AH2)
284   M1=LV(N10)+MV(N10)+1
      M2=LV(N10)-MV(N10)+1
      F2=F2*F3(M1)*F3(M2)
C     END RIGHT SIDE
C     CHECK DIMENSION OF C AND MW
      IF(NP.LT.NDIMC) GOTO 286
      WRITE (NOUT,1012) NP,NDIMC
1012  FORMAT(" DIMENSION VON C UND MW ZU KLEIN",2I5)
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
      CL=YG(A,1.,AMUL,AZ,1.)
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
      CL=YG(A,1.,AMUL,AZ,0.)
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
      CL=YG(A,1.,AMUL,AZ,-1.)
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
1013  FORMAT(" DIMENSION VON CW UND MVM ZU KLEIN",2I5)
      STOP 11
414   IF(MAUS.LE.0) GOTO 430
      WRITE (NOUT,1014) MZAHL,IQ
1014  FORMAT(" NUMBER OF MATRIX ELEMENTS",I3," NUMBER OF ",
     1 "SIGMA-FACTORS",I3)
      WRITE(NOUT,1015)
1015  FORMAT(" MATRIX ELEM SIGMA-FACTORS")
      DO 416 MH=1,MZAHL
416   WRITE (NOUT,1016) CW(MH),(MVM(KH,MH),KH=1,IQ)
1016  FORMAT(1X,E12.4,I5,20I3)
430   CALL PORD(LL,LR,MZAHL,KZL,KZR)
      MDIM=NDIMPK
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
1017  FORMAT(" KZAHL=",I5,"      IQM=",I4,"      DIMENSION ",I6)
      IF(KZAHL.LE.0) GOTO 808
      IF(LAUS.GT.0) WRITE (NOUT,1018)
1018  FORMAT(//," AUSDRUCK DER TERME WIE AUF BAND")
                     CALL WRITAP(LAUS,IQM)
808   GOTO (810,880,810,880),NOP
C     FOR DIFFERENTIAL OPERATORS,REDUCE L-VALUE
810   L3=L3+1
      IF(L3.LE.NH6) GOTO 160
      WRITE (NOUT,879) IZREK
879   FORMAT(I10," RECORDS GESCHRIEBEN")
C     END LOOP OPERATORS
880   CONTINUE
C     END LOOP ZERLEGUNGEN RIGHT
890   CONTINUE
C     END LOOP ZERLEGUNGEN LEFT
900   CONTINUE
      WRITE (NOUT,1981)
1981  FORMAT(" ENDE DER RECHNUNG VON LUELMA")
      STOP
      END
      SUBROUTINE MVAL(N5)
C     THIS SUBROUTINE CALCULATES ALL POSSIBLE MVALUE COMBINATIONS OF
C     POLYNOMIALS
      COMMON /MVA/ KV(6),KMV(6,81),NDIMK,MMOEG,MKAUS
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
150    FORMAT(" MVALUES",10I4)
160   RETURN
198   FORMAT(" DIMENSION KMV ZU KLEIN")
199   PRINT 198
200   PRINT 201,I,MMOEG,KV,IHV
201   FORMAT(" NACH",I5," VERSUCHEN VON",I5," AUFGEGEBEN",
     1 " KV,IHV= ",12I3)
      STOP 100
      END
      SUBROUTINE SDA1(L,MW,NEL,C,JJ,ND)
C     THIS SUBROUTINE CALCULATES ALL POSSIBLE REPRESENTATIONS
C     AND DETERMINES WHAT SIGMA FACTORS EXIST
C     FOR EACH REDUCED MATRIXELEMENT
C
C     L CONTAINS ND L-VALUES
C     MW CONTAINS ALL M-VALUE COMBINATIONS
C     NEL IS THE NUMBER OF MATRIX ELEMENTS C
C     JJ IS DETERMINED IN SDA1 AND GIVES THE NUMBER OF SIGMA COMBINATIONS
      COMMON /SDA/F8(21),CW(191),MVM(19,191),KAUS,N5,NDIMCW,IQ
     1 ,NKAPO
      COMMON /HIVA/ NOUT,NBAND,NDIMPK
      DIMENSION L(19),MW(123,19),C(123)
      DIMENSION JVM(23),JVN(23),LQ(19),IVM(23),IVN(23),IG(23)
      DIMENSION KY(19),KYS(19),IWN(23),IWM(23),MKOM(7,7)
C
C     DIMENSION JVM(IQ),JVN(IQ),IVM(IQ),IVN(IQ),IG(IQ),IWM(IQ),IWN(IQ)
C     DIMENSION KY(6*NC-5),KYS(6*NC-5),LQ(6*NC-5)
C
C     CHECK PRINTOUT
      IF(KAUS.LE.0) GOTO 20
      WRITE (NOUT,1001) NEL,ND
1001  FORMAT(" ANZAHL DER MATRIXELEMENTE",I4,"ANZAHL DER",
     1 " DREHIMPULSE",I4)
      WRITE (NOUT,1002) (L(KH),KH=1,ND)
1002  FORMAT(" DREHIMPULSE",19I3)
      DO 10 KH=1,NEL
10    WRITE (NOUT,1003) C(KH),(MW(KH,NH),NH=1,ND)
1003  FORMAT(" MATRIXELEMENT ",E12.4," M-WERTE",19I3)
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
1500   FORMAT(" LOOP UEBER M-WERTE",I5," TES M-ELEMENT")
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
1200  FORMAT(" LOOP W-INDICES, IK,JVM,JVN,IVM,IVN,IWM,IWN,IG =",8I5)
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
1320   FORMAT(" NEUES ELEMENT C1=",G15.5," IH=",I5," VM VN",20I3)
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
1610  FORMAT(I5,"TES M-ELEMENT C =",G15.5," MVM",19I3)
610   IH=IQ
      GOTO 62
999   CONTINUE
      RETURN
2000  WRITE(NOUT,1004)
1004  FORMAT(" +++ SUMME L-WERTE UNGERADE ++")
      STOP 21
2001  WRITE (NOUT,1005) NH,(MW(NH,IH),IH=1,ND)
1005  FORMAT(" ++ SUMME M-WERTE NICHT NULL, NH=",I4," M-WERTE",19I3)
      STOP 22
2002  WRITE (NOUT,1006) M,LQ(M)
1006  FORMAT(" DER",I5,"TE L-WERT IST NEGATIV",I5)
      STOP 23
2003  WRITE (NOUT,1007) IH,M,N
1007  FORMAT(" INDICES AUSSERHALB BEREICH,IH,M,N",3I5)
      STOP 24
2004  WRITE (NOUT,1008) JJ
1008  FORMAT(" DIMENSIONIERUNG VON CW ZU KLEIN JJ=",2I5)
      STOP 25
      END
      SUBROUTINE PORD(LL,LR,MZAHL,KZL,KZR)
C     THIS ROUTINE PUTS ALL MATRIXELEMENTS INTO A SCHEME,ACCORDING TO
C     THE OCCURRING SIGMA FACTORS
      COMMON /SDA/F8(21),CW(191),MVM(19,191),KAUS,N5,NDIMCW,IQ
     1 ,NKAPO
      COMMON /HIVA/ NOUT,NBAND,NDIMPK
      COMMON /ORD/ KZAHL, JQ(635),EPO(635),INDPO(635),MVK(19,635)
C     DIMENSION JQ(NDIMPK),EPO(NDIMPK),INDPO(NDIMPK),MVK(IQ,NDIMPK)
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
      IF(KZAHL.GT.NDIMPK)  GOTO 1001
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
1002  FORMAT(" DIMENSION IN PORD ZU KLEIN,KZAHL=",2I5,"TES ELEMENT")
      STOP 41
      END
      SUBROUTINE SORT2(IQ,JVM,JVN)
C     SORT2 ORDER THE INDICES JVM,JVN LEXICOGRAPHICALLY
      DIMENSION IND(23),JVM(23),JVN(23)
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
C     SORTA ORDERS ARRAY J AND GIVES THE ORDERING IN IND
      DIMENSION JQ(23),IND(23)
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
C     PUTA INTERCHANGES ELEMENTS IN JQ ACCORDING TO IND
      DIMENSION JQ(23),IND(23),KH(23)
      DO 10 I=IMIN,IMAX
10    KH(I)=JQ(I)
      DO 20 I=IMIN,IMAX
      IH=IND(I)
20    JQ(I)=KH(IH)
      RETURN
      END
      SUBROUTINE WRITAP(LAUS,IQM)
      COMMON /HIVA/ NOUT,NBAND,NDIMPK
      COMMON /ORD/ KZAHL, JQ(635),EPO(635),INDPO(635),MVK(19,635)
      DO 40 KH=1,KZAHL
      IF(LAUS.LE.0) GOTO 40
      II=JQ(KH)
      WRITE(NOUT,1000) II,INDPO(KH),EPO(KH),(MVK(IH,KH),IH=1,II)
1000  FORMAT(I5," SIGMA-FAKTOREN U INDEX",I8," EPO=",G15.5,
     1 "MVK",19I3)
40    CONTINUE
      WRITE (NBAND)(JQ(KH),INDPO(KH),EPO(KH),(MVK(IH,KH),IH=1,IQM),
     1 KH=1,KZAHL)
      DO 50 IH=1,IQM
      DO 45 KH=1,KZAHL
45    MVK(IH,KH)=0.
50    CONTINUE
      RETURN
      END
