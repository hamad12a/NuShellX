CBAB PROGRAM TO CONVERT INTERACTION AND SPS FILES IN ISOSPIN FORMALISM
CBAB INTO FILES IN PROTON-NEUTRON FORMALISM
CEKW Warning, The orbits are determined by reading an isospin *.sps file.
CEKW Thus, the order of orbits in *.sps must be sequential
CEKW in orbit and proton-neutron for each orbit.
      INCLUDE 'MACFOR.INC'
      INCLUDE 'DATA.INC'
 
cekw      INCLUDE 'SHGR:PARAM.INC'
cekw ******* the parameter values in this group replace shgr:param.inc
      PARAMETER (MAX_NLWRDS=8)
      PARAMETER (MAX_NSPS=MAX_NLWRDS*MAX_BIT_NLWRDS)
      PARAMETER (MAX_NJL=36)
      PARAMETER (MAX_NOSC=6)
      PARAMETER (MAX_NPART=6000)
      PARAMETER (MAX_MSCHEME=80000)
      PARAMETER (MAX_IPROT=50)
      PARAMETER (MAX_NJT=2100)
      PARAMETER (NLWRDS_BLOCK=64)
      PARAMETER (N_WORDS=MAX_NJL)
      PARAMETER (MAX_PRED=MAX_NJT)
      PARAMETER (NO_JVAL=10)
      PARAMETER (NO_TVAL=3)
      INTEGER TWELVE
      PARAMETER (TWELVE=3+NO_TVAL*NO_JVAL)
cekw ****************** end of group
 
      PARAMETER MAX_V = 200000
      CHARACTER*9 INTFIL,SPSFIL,SHELL
      CHARACTER*36 FILEG
      CHARACTER*11 INTFILPN,SPSFILPN
      CHARACTER*11 NINTFILPN,NSPSFILPN
      CHARACTER*1 LAB,LABP,LABN
      CHARACTER*3 V_TYPE
      DIMENSION LAB(9,MAX_NJL),LABP(9,MAX_NJL),LABN(9,MAX_NJL)
      DIMENSION KNJL(MAX_NJL),KNOSC(MAX_NOSC),IFIRST(MAX_NJL)
      DIMENSION KNJLPN(MAX_NJL),KNOSCPN(MAX_NOSC)
      DIMENSION KD(MAX_NSPS),KN(MAX_NSPS),KL(MAX_NSPS)
      DIMENSION KJ(MAX_NSPS),KM(MAX_NSPS),KT(MAX_NSPS)
      DIMENSION ENL(MAX_NJL),XJJ(MAX_NJL)
      DIMENSION KV(6,MAX_V),V(MAX_V)
      CHARACTER*9 CHAR_VAR(60)
      DIMENSION INT_VAR(30),REAL_VAR(30)
 
      type 113
  113 Format(1x,'******************************************************
     1*************',/,' Warning, The orbits are determined by reading
     2an isospin *.sps file.',/,' Thus the order
     3of orbits in *.sps',/,' must be sequential in orbit and
     4proton-neutron for each orbit.',/,' ****************************
     5***************************************',/)
      TYPE 112
c112   FORMAT(1X,'ENTER CODE FOR MATRIX ELEMENTS',/,1X,'P: FOR P-P TBME
c     1 ONLY',/,1X,'N: FOR N-N TBME ONLY',/,1X,'PN: FOR P-N (T=1)
c     1TBME ONLY',/,1X,'ALL: FOR ALL ME')
112   FORMAT(1X,'Enter code for matrix elements',/,1X,'p: for p-p tbme
     1 only',/,1X,'n: for n-n tbme only',/,1X,'pn: for p-n (T=1)
     1 tbme only',/,1X,'all: for all tbme')
      ACCEPT 111,V_TYPE
      if(v_type.eq.'n') v_type = 'N'
      if(v_type.eq.'p') v_type = 'P'
      if(v_type.eq.'pn') v_type = 'PN'
      if(v_type.eq.'all') v_type = 'ALL'
      TYPE 110
110   FORMAT(1X,'sps file name (A9): '$)
      ACCEPT 111,SPSFIL
111   FORMAT(A)
      TYPE 100
100   FORMAT(1X,'int file name (A9): '$)
      ACCEPT 101,INTFIL
101   FORMAT(A)
      SPSFILPN = SPSFIL//'pn'
      INTFILPN = INTFIL//'pn'
      ISPSOUT = 1
      IINTOUT = 1

      FILEG = SPSFILPN//'.sps'
      CALL FCH2(FILEG)
      TYPE 800,fileg
c800   FORMAT(/,1X,'New SPS file name = ',A36,
c     1//,1X,'TYPE RETURN IF YOU ARE HAPPY WITH THIS'
c     1/,1X,'TYPE "NONE" IF YOU DO NOT WANT NEW FILE MADE',
c     1/,1X,'TYPE "NEWNAME" IF YOU WANT TO CHANGE THE NAME'/)
800   FORMAT(/,1X,'New SPS file name = ',A36,
     1//,1X,'Type return is you are happy with this'
     1/,1X,'Type "none" of you do not want new file name',
     1/,1X,'Type "newname" if
     1 you want to change the name to "newname"'/)
      ACCEPT 111,NSPSFILPN
      IF(NSPSFILPN.EQ.'NONE') ISPSOUT = 0
      IF(NSPSFILPN.EQ.'none') ISPSOUT = 0
      IF(NSPSFILPN.EQ.'NONE') GO TO 801
      IF(NSPSFILPN.EQ.'none') GO TO 801
      IF(NSPSFILPN.NE.' ') SPSFILPN = NSPSFILPN
801   CONTINUE

      FILEG = INTFILPN//'.int'
      CALL FCH2(FILEG)
      TYPE 810,fileg
c810   FORMAT(/,1X,'New INT file name = ',A36,
c     1//,1X,'TYPE RETURN IF YOU ARE HAPPY WITH THIS'
c     1/,1X,'TYPE "NONE" IF YOU DO NOT WANT NEW FILE MADE',
c     1/,1X,'TYPE "NEWNAME" IF YOU WANT TO CHANGE THE NAME'/)
810   FORMAT(/,1X,'New INT file name = ',A36,
     1//,1X,'Type return is you are happy with this'
     1/,1X,'Type "none" of you do not want new file name',
     1/,1X,'Type "newname" if
     1 you want to change the name to "newname"'/)
      ACCEPT 111,NINTFILPN
      IF(NINTFILPN.EQ.'NONE') IINTOUT = 0
      IF(NINTFILPN.EQ.'none') IINTOUT = 0
      IF(NINTFILPN.EQ.'NONE') GO TO 811
      IF(NINTFILPN.EQ.'none') GO TO 811
      IF(NINTFILPN.NE.' ') INTFILPN = NINTFILPN
811   CONTINUE
 
CBAB      OPEN (UNIT=10,FILE=SPSFIL//'.SPS',STATUS='OLD')
CBAB      OPEN (UNIT=11,FILE=INTFIL//'.INT',STATUS='OLD')
      CALL OPENSPS(10,SPSFIL//'.SPS',0)
      CALL OPENINT(11,AREA_SPS,AREA_INT,INTFIL//'.INT',0,IERR)
      FILEG = SPSFILPN//'.sps'
      CALL FCH2(FILEG)
      IF(ISPSOUT.EQ.1)OPEN (UNIT=20,FILE=FILEG)
      FILEG = INTFILPN//'.int'
      CALL FCH2(FILEG)
      IF(IINTOUT.EQ.1) OPEN (UNIT=21,FILE=FILEG)
 
CBAB READ OLD SPS FILE
      READ(10,*) NSPS,NJL,(KNJL(I),I=1,NJL),NOSC,(KNOSC(I),I=1,NOSC),IAC
     1,IZC
CEKW      READ(10,200) ((LAB(J,I),J=1,9),I=1,NJL)
      CALL FREE_REAL_A(10,NJL,CHAR_VAR,0,REAL_VAR,READ_ERR)
      DO 3 I=1,NJL
      SHELL=CHAR_VAR(I)
      DO 4 II = 1,7
    4 LAB(II,I) = SHELL(II:II)
    3 CONTINUE
C      READ(10,200) ((LAB(J,I),J=1,7),I=1,NJL)
C200   FORMAT(200A1)
      DO 201 I = 1,NSPS
      READ(10,*,END=202) KD(I),KN(I),KL(I),KJ(I),KM(I),KT(I)
201   CONTINUE
202   CONTINUE
 
      II = 1
      DO 220 I = 1,NJL
      XJJ(I) = KJ(II)
      XJJ(I) = XJJ(I)/2.
      II = II + KNJL(I)
220   CONTINUE
 
CBAB MAKE NEW SPS FILE
      NSPSPN = NSPS
      NJLPN = 2*NJL
CBAB      NOSCPN = 2*NOSC
      NOSCPN = 2
 
      DO 321 I = 1,NJL
CEKW      DO 322 J = 1,9
      DO 322 J = 1,7
      IF(LAB(J,I).NE.' ') IFIRST(I) = J
      IF(LAB(J,I).NE.' ') GO TO 321
322   CONTINUE
321   CONTINUE
 
      MAX_JPN = 0
      DO 301 I = 1,NJL
      KNJLPN(I)=KNJL(I)/2
      LABP(1,I) = 'P'
      LABN(1,I) = 'N'
      JJ = 1
CEKW      DO 311 J = IFIRST(I),9
      DO 311 J = IFIRST(I),7
      JJ = JJ+1
      IF(LAB(J,I).NE.' '.AND.JJ.GT.MAX_JPN) MAX_JPN=JJ
      LABP(JJ,I) = LAB(J,I)
      LABN(JJ,I) = LAB(J,I)
311   CONTINUE
301   CONTINUE
CEKW      IF(ISPSOUT.EQ.1)
      IF(ISPSOUT.EQ.1.AND.NJL.LT.10)
     1WRITE(20,300) NSPSPN,NJLPN,(KNJLPN(I),I=1,NJL),(KNJLPN(I),I=1,NJL)
     1,NOSCPN,NJL,NJL,IAC,IZC
      IF(ISPSOUT.EQ.1.AND.NJL.GT.9)
     1WRITE(20,307) NSPSPN,NJLPN,(KNJLPN(I),I=1,NJL),(KNJLPN(I),I=1,NJL)
     1,NOSCPN,NJL,NJL,IAC,IZC
  300 FORMAT(1X,I3,1X,I3,1X,<NJL>(1X,I2),1X,<NJL>(1X,I2),2X,I3,
     1<NOSCPN>(1X,I2),2X,I3,2X,I3)
  307 FORMAT(1X,I3,1X,I3,1X,<NJL>(1X,I2),1X,<NJL>(1X,I2),' - ',/,2X,I3,
     1<NOSCPN>(1X,I2),2X,I3,2X,I3)
CEKW      IF(ISPSOUT.EQ.1) WRITE(20,305) ((LABP(J,I),J=1,8),I=1,NJL),
CEKW     1((LABN(J,I),J=1,8),I=1,NJL)
      IF(ISPSOUT.EQ.1)WRITE(20,305)((LABP(J,I),J=1,7),I=1,NJL)
      IF(ISPSOUT.EQ.1)WRITE(20,306)((LABN(J,I),J=1,7),I=1,NJL)
305   FORMAT(<7*NJL>(A1),' - ')
306   FORMAT(<7*NJL>(A1))
      ID = 0
      DO 302 I = 1,NSPS
      IF(KT(I).EQ.-1) ID = ID + 1
      IF(ISPSOUT.EQ.1.AND.KT(I).EQ.-1)
     1WRITE(20,303) ID,KN(I),KL(I),KJ(I),KM(I),KT(I)
303   FORMAT(1X,I6,10I5)
302   CONTINUE
      DO 304 I = 1,NSPS
      IF(KT(I).EQ.1) ID = ID + 1
      IF(ISPSOUT.EQ.1.AND.KT(I).EQ.1)
     1WRITE(20,303) ID,KN(I),KL(I),KJ(I),KM(I),KT(I)
304   CONTINUE
 
CBAB READ OLD INT FILE
 
      READ(11,*) NV,(ENL(I),I=1,NJL)
      INV = 0
401   CONTINUE
      INV = INV + 1
      IF(INV.GT.MAX_V) TYPE 410
      IF(INV.GT.MAX_V) STOP
410   FORMAT(1X,'ERROR: INCREASE MAX_V')
      READ(11,*,END=402) (KV(I,INV),I=1,6),V(INV)
      GO TO 401
402   CONTINUE
      INV = INV-1
      IF(INV.NE.NV) TYPE 403,NV,INV
403   FORMAT(1X,'WARNING: FIRST # OF SECOND LINE OF *.INT FILE (NV)'
     1/,'DIFFERS FROM # OF MATRIX ELEMENTS READ IN (INV)',
     1/,'IN,INV = ',2I5)
 
CBAB OUTPUT TO NEW INT FILE
      IF(IINTOUT.EQ.1.AND.NJL.LE.14)THEN
      WRITE(21,405) (I,(LABP(J,I),J=1,MAX_JPN),I=1,NJL)
      WRITE(21,405) (I+NJL,(LABN(J,I),J=1,MAX_JPN),I=1,NJL)
      ENDIF
  405 FORMAT(1X'!',
     1<NJL>(I2,'=',<MAX_JPN>(A1)))
C     1<NJL>(I2,'=',<MAX_JPN>(A1)),<NJL>(I2,'=',<MAX_JPN>(A1)))
 
      INVPN = 0
      IF(V_TYPE.EQ.'ALL')THEN
      DO 600 I = 1,INV
      IF(KV(6,I).EQ.1) INVPN = INVPN+1
600   CONTINUE
      DO 602 I = 1,INV
      IF(KV(6,I).EQ.1) INVPN = INVPN+1
602   CONTINUE
      DO 603 I = 1,INV
      INVPN = INVPN+1
      IF(KV(1,I).NE.KV(2,I)) INVPN = INVPN+1
      IF(KV(3,I).NE.KV(4,I)) INVPN = INVPN+1
      IF(KV(1,I).NE.KV(2,I).AND.KV(3,I).NE.KV(4,I)) INVPN = INVPN+1
603   CONTINUE
      GOTO 700
      ELSE
      IF((V_TYPE.EQ.'P').OR.(V_TYPE.EQ.'N'))THEN
      INVPN=INV
      GOTO 690
      END IF
      IF(V_TYPE.EQ.'PN')THEN
      DO 604 I = 1,INV
CEKW      IF(KV(6,I).EQ.1)INVPN = INVPN+1
      INVPN = INVPN+1
      IF(KV(1,I).NE.KV(2,I).AND.KV(6,I).EQ.1) INVPN = INVPN+1
      IF(KV(3,I).NE.KV(4,I).AND.KV(6,I).EQ.1) INVPN = INVPN+1
      IF(KV(1,I).NE.KV(2,I).AND.KV(3,I).NE.KV(4,I).AND.KV(6,I).EQ.1)
     1INVPN = INVPN+1
  604 CONTINUE
      END IF
      END IF
 
  690 DO I=1,2*NJL
      ENL(I)=0.0
      END DO
 
 
  700 CONTINUE
      IF(IINTOUT.EQ.1.AND.NJL.LE.10)THEN
      WRITE(21,520) INVPN,(ENL(I),I=1,NJL),(ENL(I),I=1,NJL)
      ENDIF
      IF(IINTOUT.EQ.1.AND.NJL.GT.10)THEN
      WRITE(21,521) INVPN,(ENL(I),I=1,NJL),(ENL(I),I=1,NJL)
      ENDIF
  520 FORMAT(I6,5(1X,F7.3),' - ',/,6x,5(1X,F7.3),' - ',
     1/,6x,5(1X,F7.3),' - ',/,6x,5(1X,F7.3))
  521 FORMAT(I6,5(1X,F10.6),' - ',/,<2*NJL/5>(6X,5(1X,F10.6),' - ',/))
CEKW 520  FORMAT(1X,I5,20(1X,F7.3))
CEKW      TYPE 521,NJL,IINTOUT,INV,V_TYPE,KV(6,1)
      IF((V_TYPE.EQ.'ALL').OR.(V_TYPE.EQ.'P'))THEN
      DO 500 I = 1,INV
      IF(IINTOUT.EQ.1.AND.KV(6,I).EQ.1) WRITE(21,501) (KV(J,I),J=1,4),
     1(KV(J,I),J=5,6),V(I)
  501 FORMAT(1X,4I3,1X,2I3,1X,F13.5)
  509 FORMAT(1X,4I3,1X,2I3,1X,F13.5)
  500 CONTINUE
      END IF
      IF((V_TYPE.EQ.'ALL').OR.(V_TYPE.EQ.'N'))THEN
      DO 502 I = 1,INV
      IF(IINTOUT.EQ.1.AND.KV(6,I).EQ.1)THEN
      WRITE(21,501) (KV(J,I)+NJL,J=1,4),(KV(J,I),J=5,6),V(I)
      ENDIF
502   CONTINUE
      END IF
CEKW ********* ILIST = 0 IF THE RESTRICTION II(1).LE.II(3) IS DESIRED  ***
      ILIST = 1
      INVPNC = INVPN
      IF((V_TYPE.EQ.'ALL').OR.(V_TYPE.EQ.'PN'))THEN
      DO 503 I = 1,INV
      IPRINT = 0
      IF(V_TYPE.EQ.'PN')GOTO 550
      IF(KV(1,I).GT.KV(3,I).AND.ILIST.EQ.0)IPRINT = 1
      INVPNC = INVPNC - IPRINT
      IF(IINTOUT.EQ.1.AND.IPRINT.EQ.0)THEN
      WRITE(21,501) KV(1,I),KV(2,I)+NJL,KV(3,I),KV(4,I)+NJL,
     1(KV(J,I),J=5,6),V(I)
      ENDIF
      IPRINT = 0
      XJT = KV(5,I)
      XTT = KV(6,I)
      P12 = PAR(XJJ(KV(1,I))+XJJ(KV(2,I))-XJT-XTT)
      PV12 = P12*V(I)
      IF(KV(2,I).GT.KV(3,I).AND.ILIST.EQ.0)IPRINT = 1
      INVPNC = INVPNC - IPRINT
      IF(IINTOUT.EQ.1.AND.KV(1,I).NE.KV(2,I).AND.IPRINT.EQ.0)THEN
      WRITE(21,501) KV(2,I),KV(1,I)+NJL,KV(3,I),KV(4,I)+NJL,
     1(KV(J,I),J=5,6),PV12
      ENDIF
      IPRINT = 0
      P34 = PAR(XJJ(KV(3,I))+XJJ(KV(4,I))-XJT-XTT)
      PV34 = P34*V(I)
      IF(KV(1,I).GT.KV(4,I).AND.ILIST.EQ.0)IPRINT = 1
      INVPNC = INVPNC - IPRINT
      IF(IINTOUT.EQ.1.AND.KV(3,I).NE.KV(4,I).AND.IPRINT.EQ.0)THEN
      WRITE(21,501) KV(1,I),KV(2,I)+NJL,KV(4,I),KV(3,I)+NJL,
     1(KV(J,I),J=5,6),PV34
      ENDIF
      IPRINT = 0
      PV1234 = P12*P34*V(I)
      IF(KV(2,I).GT.KV(4,I).AND.ILIST.EQ.0)IPRINT = 1
      INVPNC = INVPNC - IPRINT
      IF(IINTOUT.EQ.1.AND.KV(1,I).NE.KV(2,I).AND.KV(3,I).NE.KV(4,I)
     1.AND.IPRINT.EQ.0)THEN
      WRITE(21,501) KV(2,I),KV(1,I)+NJL,KV(4,I),KV(3,I)+NJL,
     1(KV(J,I),J=5,6),PV1234
      ENDIF
      IPRINT = 0
      GOTO 503
  550 IF(KV(1,I).GT.KV(3,I).AND.ILIST.EQ.0)IPRINT = 1
cekw      IF(IINTOUT.EQ.1.AND.KV(6,I).EQ.1.AND.IPRINT.EQ.0)THEN
      INVPNC = INVPNC - IPRINT
      IF(IINTOUT.EQ.1.and.IPRINT.EQ.0)THEN
      IOPT = 1
      WRITE(21,509) KV(1,I),KV(2,I)+NJL,KV(3,I),KV(4,I)+NJL,
     1(KV(J,I),J=5,6)
      ENDIF
      IPRINT = 0
      XJT = KV(5,I)
      XTT = KV(6,I)
      P12 = PAR(XJJ(KV(1,I))+XJJ(KV(2,I))-XJT-XTT)
      PV12 = P12*V(I)
      IF(KV(2,I).GT.KV(3,I).AND.ILIST.EQ.0)IPRINT = 1
      INVPNC = INVPNC - IPRINT
      IF(IPRINT.EQ.0)THEN
cekw      IF(IINTOUT.EQ.1.AND.KV(1,I).NE.KV(2,I).AND.KV(6,I).EQ.1)THEN
      IF(IINTOUT.EQ.1.AND.KV(1,I).NE.KV(2,I))THEN
      IOPT = 2
      WRITE(21,509) KV(2,I),KV(1,I)+NJL,KV(3,I),KV(4,I)+NJL,
     1(KV(J,I),J=5,6)
      ENDIF
      ENDIF
      IPRINT = 0
      P34 = PAR(XJJ(KV(3,I))+XJJ(KV(4,I))-XJT-XTT)
      PV34 = P34*V(I)
      IF(KV(1,I).GT.KV(4,I).AND.ILIST.EQ.0)IPRINT = 1
      INVPNC = INVPNC - IPRINT
      IF(IPRINT.EQ.0)THEN
cekw      IF(IINTOUT.EQ.1.AND.KV(3,I).NE.KV(4,I).AND.KV(6,I).EQ.1)THEN
      IF(IINTOUT.EQ.1.AND.KV(3,I).NE.KV(4,I))THEN
      IOPT = 3
      WRITE(21,509) KV(1,I),KV(2,I)+NJL,KV(4,I),KV(3,I)+NJL,
     1(KV(J,I),J=5,6)
      ENDIF
      ENDIF
      IPRINT = 0
      PV1234 = P12*P34*V(I)
      IF(KV(2,I).GT.KV(4,I).AND.ILIST.EQ.0)IPRINT = 1
      INVPNC = INVPNC - IPRINT
      IF(IINTOUT.EQ.1.AND.KV(1,I).NE.KV(2,I).AND.KV(3,I).NE.KV(4,I).
     1AND.IPRINT.EQ.0)THEN
cekw     1AND.KV(6,I).EQ.1.AND.IPRINT.EQ.0)THEN
      IOPT = 4
      WRITE(21,509) KV(2,I),KV(1,I)+NJL,KV(4,I),KV(3,I)+NJL,
     1(KV(J,I),J=5,6)
      ENDIF
  503 CONTINUE
      IF(ILIST.EQ.0)TYPE 998
      TYPE 999,INVPN,INVPNC
  998 FORMAT(1X,'ILIST = 0 which means only "pn" tbme with
     1I1.le.i3 are printed')
  999 FORMAT(1X,'original # tbme = ',I5,',',2x,'final # tbme = ',I5)
      END IF
      STOP
      END

