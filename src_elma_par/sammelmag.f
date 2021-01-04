      PROGRAM SAMMEL
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)

      INCLUDE 'par/QUAL'
      
C     NZOPER: ANZAHL DER OPERATOREN IN QUAF
C     NZFMAX:     "      "     "  ZERLEGUNGEN
C     NZRHOM:     "      "     "  BASISVEKTOREN PRO ZERLEGUNG
C     MZPARM:     "      "     "  RADIALPARAMETER
C
      parameter (lind=(nzfmax*(nzfmax+1))/2*nzoper)
      integer nteil
C
C
      dimension  index(lind),mop(lind),mmfr(lind)
      DIMENSION  NUM(NZRHOM,NZFMAX),kmax(nzfmax),
     *           MZPAR(NZFMAX), NZRHO(NZFMAX),kmin(nzfmax)
      DIMENSION  DM(NDIM,NDIM,NDIM2,nzfmax,nzoper),
     *           F(NDIM,NDIM,nzfmax,nzoper)
      DIMENSION  LAHI(NZBMAX*NZBMAX,nzfmax,nzoper),
     *           IND(NZRHOM,NZRHOM) 
      DIMENSION  IOP(0:NZFMAX)
      DIMENSION  LREG(NZOPER)
C
      character*80  fnumber, qname


C
      OPEN(UNIT=4,STATUS='SCRATCH',FORM='UNFORMATTED')
      OPEN(UNIT=10,FILE='QUAOUT',STATUS='OLD',FORM='UNFORMATTED')
      OPEN(UNIT=14,FILE='FINDEX',STATUS='OLD',FORM='UNFORMATTED')
      OPEN(UNIT=15,FILE='INSAM',STATUS='OLD',FORM='FORMATTED')
      OPEN(UNIT=16 ,FILE='OUTPUT',
     $        STATUS='REPLACE', FORM='FORMATTED')
      OPEN(UNIT=17 ,FILE='OUTBAND')

      INPUT=15
      NOUT=16
      noutb=17
      WRITE(NOUT, 1111)
      
1111  FORMAT('1 SAMMEL VERSION VOM 15.06.20')
C
      nband7 = 4
      nband1 = 10
      read(input,1002) jfilmax,naus
      I=0
      WRITE (NOUT,*) ' Matrizen in files OUTDM.',jfilmax
      read(input,1002) mflu,mflo
      write(nout,*) 'Es wird von Zerlegung ',mflu,' bis ',mflo,
     *    ' aufgesammelt'
      REWIND NBAND1
      READ(NBAND1) NZF,MUL,(LREG(K),K=1,NZOPER),
     *             NZBASV2,(NZRHO(K),K=1,NZF)
      write(nout,*)'NZF,MUL,lreg:',NZF,MUL,(LREG(K),K=1,NZOPER)
      write(nout,*)'NZBASV2,NZRHO:',NZBASV2,(NZRHO(K),K=1,NZF)
      DO 22  MFL = 1,NZF
      KK=NZRHO(MFL)
 1002 FORMAT(20I3)
      I=I+KK
   22 CONTINUE
      I=0
      DO 950  K = 1,NZF
      DO 4536 N=1,NZRHO(K)
      i=i+1
      NUM(N,K)=i
C reads N3, only and does not bother about mass, charge, etc. which are on this line, too
      READ(NBAND1) MZPAR(K)
 4536 CONTINUE
  950 CONTINUE
      close(unit=nband1, status='keep')
      READ(14) NGER,(INDEX(I),I=1,NGER-1)
      write(nout,*) 'nger',  NGER,(INDEX(I),I=1,NGER-1)
      nteil = 0
      iop(0)= 0
      istop = 0
      kstop = 0
      DO 140  MFL = 1,MFLO
       kmin(mfl)=jfilmax
       kmax(mfl)=1
      DO 141 MFR=1,MFL
      DO 142 MKC=1, NZOPER
         nnn = nzrho(mfl)*nzrho(mfr)
      if(lreg(mkc)*nnn .eq. 0) goto 142
      nteil=nteil+1
      if(index(nteil).eq.-1) then
         write(nout,*) ' Fuer Zerlegung links ',mfl,' und rechts ',mfr,
     *   'ist operator ',mkc,' nicht gerechnet, nteil=',nteil
          if(naus.eq.0) stop 'dm fehlt'
          istop=istop+1
      endif
      if(index(nteil).gt.jfilmax) then
         write(nout,*) ' Fuer Zerlegung links ',mfl,' und rechts ',mfr,
     *   'ist operator ',mkc,' nicht korrekt, nteil=',nteil
          if(naus.eq.0) stop 'dm falsch'
          kstop=kstop+1
      endif
      mop(nteil) = mkc
      mmfr(nteil)= mfr
      if(index(nteil).gt.0) then
      kmin(mfl)=min(index(nteil),kmin(mfl))
      kmax(mfl)=max(index(nteil),kmax(mfl))
      endif
C MKC
142   CONTINUE
C MFR
141   CONTINUE
      iop(mfl)=nteil
C MFL         
140   CONTINUE

      if(istop.gt.0) then
         write(nout,*) ' Es fehlen ',istop,' operatoren'
         stop 'dm fehlt'
      endif
      if(kstop.gt.0) then
         write(nout,*) ' Es sind ',kstop,' operatoren inkorrekt'
         stop 'dm falsch'
      endif
       write(nout,*) 'iop' ,(iop(ix) ,ix=0,nzf)
       write(nout,*) 'kmin',(kmin(ix),ix=1,nzf)
       write(nout,*) 'kmax',(kmax(ix),ix=1,nzf)
       write(nout,*) 'mop' ,(mop(ix) ,ix=1,nger)
       write(nout,*) 'mmfr',(mmfr(ix),ix=1,nger)
C
       DO 447 MFL=MFLU,MFLO
       CLOSE(UNIT=11,STATUS='KEEP')
            write(fnumber,*) MFL
            do i=1, 255
               if(fnumber(i:i).ne.' ') goto  297
            end do
  297       do j=i, 255
               if(fnumber(j:j).eq.' ') goto 298
            end do
  298       qname = 'TQUAOUT.' // fnumber(i:j)
            open(unit=11, file=qname, form='unformatted',
     $           STATUS='REPLACE')

      do 802 mkc=1,nzoper
      do 802 mfr=1,nzfmax
      do 803 i=1,nzbmax*nzbmax
803        lahi(i,mfr,mkc)=0.
      do 806 j1=1,ndim
      do 806 i1=1,ndim
      do 805 k1=1,ndim2
805     dm(i1,j1,k1,mfr,mkc)=0.
806     f(i1,j1,mfr,mkc)=0.
802   continue

       IFERTIG=0
       DO 312 JLIM=IOP(MFL-1)+1,IOP(MFL)
            IF(INDEX(JLIM).GT.0)IFERTIG=IFERTIG+1
312    CONTINUE

      write(nout,*)'DMOUT.',kmin(mfl),' bis DMOUT.',kmax(mfl),' used'
      write(nout,*)'iop(mfl-1,mfl),ifertig',iop(mfl-1),iop(mfl),ifertig

C vv
      DO 317 IDMCOUNT=KMIN(MFL),KMAX(MFL)
      NBAND=40+IDMCOUNT
      CLOSE(UNIT=NBAND,STATUS='KEEP')
      write(fnumber,*) IDMCOUNT
      do i=1, 255
         if(fnumber(i:i).ne.' ') goto  997
      end do   
  997 do j=i, 255
         if(fnumber(j:j).eq.' ') goto 998
      end do
  998 qname = 'DMOUT.' // fnumber(i:j)
      open(unit=NBAND, file=qname, form='unformatted',
     $           STATUS='OLD')
131   READ(NBAND,END=317,ERR=313 ) MTEIL,JCOUNT,INDEXR
      GOTO 314
313   CONTINUE
314   IF(INDEXR.EQ.0)  GOTO 131
      write(nout,*) MTEIL,IOP(MFL),indexr,index(mteil)
      IF(MTEIL.LT.(IOP(MFL-1)+1) .OR. MTEIL.GT.IOP(MFL)
     $                           .or.indexr.ne.index(mteil)) THEN
            READ (NBAND,END=131)
            READ (NBAND,END=131)
            READ (NBAND,END=131)
            GOTO 131
      ELSE
        read(NBAND,ERR=317) (lahi(i,mmfr(mteil),mop(mteil)),
     *                          i=1,nzbmax*nzbmax)
        read(NBAND,ERR=317)
     $   ((f(i1,j1,mmfr(mteil),mop(mteil)),i1=1,ndim),j1=1,ndim)
        read(NBAND,ERR=317)
     $   (((dm(i1,j1,k1,mmfr(mteil),mop(mteil)),
     $         i1=1,ndim),j1=1,ndim),k1=1,ndim2)

C        ifertig= ifertig - 1
c jk two statements flipped following SAMMEL-uix.f
          if(KMIN(MFL).EQ.KMAX(MFL) .AND. ifertig.eq.0) goto 318
c           if( ifertig.eq.0) goto 318
          GOTO 131
      ENDIF
317   CONTINUE
C vv
C ----- here I should be as knowledgable as at the end of the serial <76> loop
C ----- in qual.f
C
318   IRHO=NZRHO(MFL)
      IK1=MZPAR(MFL) 

      jteil=iop(mfl-1)
      DO 341 MFR=1,MFL
      JRHO=NZRHO(MFR)
      JK1=MZPAR(MFR)


      DO 342 MKC=1, NZOPER

      nnn = nzrho(mfl)*nzrho(mfr)
      if(lreg(mkc)*nnn .eq. 0) goto 342
      jteil = jteil + 1
 
C -- 'write a reasonable NULL record', i.e., an (IND(M,N) = 0) line
C -- if ...
      IF(INDEX(JTEIL).EQ.0) then
      II1=1
      A=0.
      DO 580 M=1,IRHO
         NUML=NUM(M,MFL)
         DO 581 N=1,JRHO
            NUMR=NUM(N,MFR)
            IF( NUML.LT.NUMR) GOTO 581
c            WRITE (noutb,*) ((A, NN=1, JRHO), MM=1, IRHO)
            WRITE (11) ((A, NN=1, JRHO), MM=1, IRHO)

 581     CONTINUE
 580  CONTINUE
      goto 342
      endif

      rewind nband7

C >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
      DO 480 M=1,IRHO 
      MM=M-1+MKANA
      NUML=NUM(M,MFL) 
      DO 481 N=1,JRHO
      NN=N-1+NKANA
      NUMR=NUM(N,MFR)

      IF(NUML.LT.NUMR) GOTO 481
      NUHILF = NZBMAX*(NUML-1) + NUMR
      LL1 = LAHI(NUHILF,mfr,mkc)
      
      M1=(M-1)*IK1+1
      M2=M*IK1
      N1=(N-1)*JK1+1
      N2=N*JK1
      II1 = 1
      A = 0.


      DO 730, IR1=1, IRHO
       DO 720, IS1=1, JRHO
        IND(IR1,IS1)=0
720    CONTINUE
730   CONTINUE

      DO 510 K=M1,M2
      DO 510 L=N1,N2
      DO 510 I=1,LL1
 510  A = A + ABS(DM(K,L,I,mfr,mkc))

      IF (A.NE.0.) IND(MM,NN)=1

      IF(IND(MM,NN).EQ.0) GOTO 481
      WRITE (NBAND7) NUML, NUMR, IK1, JK1, LL1,
     *     ((F(K,L,mfr,mkc),(DM(K,L,J,mfr,mkc),J=1,LL1),L=N1,N2)
     *                                        ,K=M1,M2)

      IF(NAUS.LT.1) GOTO 481
      DO 520 K=M1,M2
      DO 520 L=N1,N2
520   WRITE(NOUT,1051) F(K,L,mfr,mkc),(J-1,DM(K,L,J,mfr,mkc),J=1,LL1)
      
  481 CONTINUE
  480 CONTINUE
C >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
      WRITE (11) ((IND(MM,NN), NN=1, JRHO), MM=1, IRHO)
      WRITE (noutb,'(20I3)') ((IND(MM,NN), NN=1, JRHO), MM=1, IRHO)
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
     *     ((F(K,L,mfr,mkc),(DM(K,L,J,mfr,mkc),J=1,LL1),L=N1,N2)
     *                                        ,K=M1,M2)
      WRITE (11)    NUML, NUMR, IK1H, JK1H, LL1,
     *    ((F(K,L,mfr,mkc),(J-1,DM(K,L,J,mfr,mkc),J=1,LL1),L=N1,N2),
     *                                             K=M1,M2)
      WRITE (noutb,'(5I3)') NUML, NUMR, IK1H, JK1H, LL1
      WRITE (noutb,'(8E12.4)') ((F(K,L,mfr,mkc),L=N1,N2),K=M1,M2)
      WRITE (noutb,'(8E12.4)') 
     *  (((real(J-1),DM(K,L,J,mfr,mkc),J=1,LL1),L=N1,N2),K=M1,M2)
C
C     ENDE LOOP BASISVEKTOREN RECHTS
910   CONTINUE
C
C     ENDE LOOP BASISVEKTOREN LINKS
920   CONTINUE


C        LOOP OPERATOREN
342        continue
C loop MFR
341        continue
C loop MFL
  447 CONTINUE
      WRITE(NOUT,3011)
1051  FORMAT(' F= ',E12.5,4(' IQ= ',I3,' DM= ',E12.5)) 
3011  FORMAT(//,' ENDE DER RECHNUNG VON SAMMEL')
      CLOSE(UNIT=6)
      CLOSE(UNIT=12,STATUS='KEEP')


 666  STOP

      END