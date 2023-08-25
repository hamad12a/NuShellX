module TLMB
IMPLICIT NONE

REAL*8,PRIVATE BIN(0:99,0:99),TIN(0:99,0:99,0:99)


contains

! W D M Rae Garsington 2007. Interface to TMSH

  function moshinsky(cN,lam,sn1,sl1,sn2,sl2,m1,m2)
     implicit none
     integer,intent(in):: cN,lam,sn1,sl1,sn2,sl2,m1,m2
     real*8:: moshinsky,d
     logical,save::first=.true.
     integer:: ee,e1,e2
     
     if (first) then
        first=.false.
        call init()
     end if   
     if (2*cN+lam/=2*(sn1+sn2)+sl1+sl2) then
        moshinsky=real(0,8)
        return
     end if
     if (cN<0.or.lam<0.or.sn1<0.or.sn2<0.or.sl1<0.or.sl2<0) then
        moshinsky=real(0,8)
        return
     end if

          e1=2*sn1+sl1
          e2=2*sn2+sl2
          ee=2*cN+lam
     
     if (m1>=m2) then
          d=real(m1,8)/real(m2,8)
          moshinsky=TMB(ee,lam,0,0,e1,sl1,e2,sl2,lam,d)
     else     
          d=real(m2,8)/real(m1,8)
          moshinsky=TMB(ee,lam,0,0,e2,sl2,e1,sl1,lam,d)
     end if
     end function moshinsky
     


!** Version 1.0: May 2001 
!   Fortran 95 Version W D M rae, Garsington, 2007
!   http://knollhouse.org
!** E-mail: Gintautas_Kamuntavicius@fc.vdu.lt
!** Reference: nucl-th/0105009
!** WWW: http://www.nuclear.physics.vdu.lt 
!   Global Variables to Module BIN TIN
!** BIN-array of binomial coefficients
!** TIN-array of trinomial coefficients
!** TMB-real*8 function of HOBs:
!** Input: EE-centre of mass energy, LL-centre of mass angular moment
!**        ER-relative energy, LR-relative angular moment!
!**        E1-energy of the first particle, L1-angular moment of the first particle
!**        E2-energy of the second particle, L2-angular moment of the second particle
!**        LM-total angular moment   
!** Output: TMB-HARMONIC OSCILLATOR BRACKET

    SUBROUTINE INIT()
    IMPLICIT NONE
    REAL tim1,tim2,timd
	REAL*8 T1,skirtm
    INTEGER LAST
!c	LAST=0,1,2,3,4...
    LAST=3
!	CALL CPU_TIME(tim1)
    CALL BINOM
	CALL TRINOM
    END SUBROUTINE INIT

      SUBROUTINE BINOM
!  	  THE ARRAY OF BINOMIAL COEFFICIENTS
!     BIN(I,J)= = I!/J!/(I-J)! 
      IMPLICIT NONE
      INTEGER I,K
    	DO I=0,99
	        BIN(I,0)=1.D0
	        BIN(I,I)=1.D0
	        DO K=1,I/2
                BIN(I,K)=DNINT(BIN(I,K-1)/DFLOAT(K)*DFLOAT(I-K+1))
	            BIN(I,I-K)=BIN(I,K)
	        END DO
        END DO
	    RETURN
      END SUBROUTINE BINOM

      INTEGER FUNCTION TRI(I,J,K)
!     TRIADIC CONDITION FOR MOMENTS I/2,J/2,K/2:
!     I+J>=K, I+K>=J, J+K>=I,
!     I/2+J/2+K/2 = INTEGER.
!     TRI=1, WHEN TRIADIC CONDITION IS FULFILLED, TRI=0 OTHERWISE	 
      IMPLICIT NONE
      INTEGER I,J,K,L
        TRI=0
    	L=I+J+K
        IF(L/2*2.NE.L) RETURN
	    L=L/2
        IF((L-I)*(L-J)*(L-K).LT.0) RETURN
        TRI=1
      RETURN
      END FUNCTION TRI

      REAL*8 FUNCTION C6J(I,J,K,L,M,N)
!     6J - COEFFICIENT
!     ( I/2  J/2  K/2 )
!     ( L/2  M/2  N/2 ) 
!     [JB 65] (22.1.4)
	  IMPLICIT NONE
      INTEGER I,J,K,L,M,N,I1,I2,I3,I4,I5,IZ,JZ,KZ
 	  REAL*8 T,DZ
	  C6J=0.D0
	  IF(TRI(I,J,K)*TRI(I,M,N)*TRI(J,L,N)*TRI(K,L,M).EQ.0) RETURN
	  I1=(I+J+K)/2
      I2=(I+M+N)/2
      I3=(J+L+N)/2
      I4=(K+L+M)/2
      I5=(I+K+L+N)/2
	  T=DSQRT(DFLOAT((I1+1)*(I4+1))/DFLOAT((I2+1)*(I3+1))*&
         BIN(I2,I)*BIN(I,I2-N)*BIN(I4,L)/&
         BIN(I3,N)*BIN(L,I4-K)*BIN(I1,K)*BIN(K,I1-J)/&
         BIN(N,I3-L))/DFLOAT(I2-N+1)
	  JZ=MAX0(0,(I+L-J-M)/2)
	  DZ=1.D0
	  IF((JZ+I5)/2*2.NE.(JZ+I5)) DZ=-1.D0
      KZ=MIN0(I2-M,I3-J,I5-K)
	  DO IZ=JZ,KZ
	      C6J=C6J+DZ*T*BIN(I2-M,IZ)/BIN(I2,I5-K-IZ)*&
              BIN(I2-I,I3-J-IZ)/BIN(I4-I3+J+IZ+1,I2-N+1)
          DZ=-DZ
      END DO
	  RETURN
      END FUNCTION C6J

      REAL*8 FUNCTION C9J(J1,J2,J3,L1,L2,L3,K1,K2,K3)
!     9J COEFICIENT
!     (J1/2 J2/2 J3/2)
!     (L1/2 L2/2 L3/2)
!  	  (K1/2 K2/2 K3/2)  	 
!     [JB 65] (24.33)
      IMPLICIT NONE
      INTEGER J1,J2,J3,L1,L2,L3,K1,K2,K3,I,J,K,L
      C9J=0.D0
      L=TRI(J1,J2,J3)*TRI(L1,L2,L3)*TRI(K1,K2,K3)*&
          TRI(J1,L1,K1)*TRI(J2,L2,K2)*TRI(J3,L3,K3)
      IF(L.EQ.0) RETURN
      J=MAX0(IABS(J1-K3),IABS(J2-L3),IABS(L1-K2))
      K=MIN0(J1+K3,J2+L3,L1+K2)
      DO I=J,K,2
	      C9J=C9J+DFLOAT(I+1)*C6J(J1,J2,J3,L3,K3,I)*&
              C6J(L1,L2,L3,J2,I,K2)*C6J(K1,K2,K3,I,J1,L1)
      END DO
	  IF(J/2*2.NE.J) C9J=-C9J
 	  RETURN
      END FUNCTION C9J
      
	  REAL*8 FUNCTION KL0(I,J,K)
!  	  KLEBS-GORDAN COEFFICIENT WITH ZERO PROJECTIONS OF MOMENTA
!	  (I, J, K)
!	  (0, 0, 0)  
!	   I,J,K - MOMENTA = INTEGER NUMBERS
!	   [JB,65] (15.10)
      IMPLICIT NONE
	  REAL*8 T
      INTEGER I,J,K,L,M
      KL0=0.D0
      IF(TRI(I,J,K).EQ.0) RETURN
 	  L=(I+J+K)/2
	  M=L-K
	  T=1.D0
	  IF(M/2*2.NE.M) T=-1.D0
      KL0=T*BIN(K,L-J)*BIN(L,K)/&
           DSQRT(BIN(2*K,2*(L-J))*BIN(2*L+1,2*K+1))
      RETURN
      END FUNCTION KL0

   	  SUBROUTINE TRINOM
!	  THE ARRAY OF TRINOMIAL COEFFICIENTS
!  	  TIN(I,J,K)=I!!/J!!/K!!
	  IMPLICIT NONE
	  INTEGER I,J,K,M,N
	  TIN(0,0,0)=1.D0
	  TIN(1,1,1)=1.D0
	  DO I=2,99
	      M=I-I/2*2
	      TIN(I,I,M)=1.D0
	      TIN(I,M,I)=1.D0
	      N=M+2
	      DO J=I,N,-2
	          DO K=N,J,2
	              TIN(I,J,K)=TIN(I,J,K-2)/DFLOAT(K)
	              TIN(I,K,J)=TIN(I,J,K)
	          END DO
	          TIN(I,J-2,M)=TIN(I,J,M)*DFLOAT(J)
	          TIN(I,M,J-2)=TIN(I,J-2,M)
	      END DO
	  END DO
	  RETURN
	  END SUBROUTINE TRINOM
	
	  REAL*8 FUNCTION G(E1,L1,EA,LA,EB,LB)
	  IMPLICIT NONE
	  INTEGER E1,L1,EA,LA,EB,LB
	  G=KL0(LA,LB,L1)*DSQRT(DFLOAT((2*LA+1)*(2*LB+1))*&
        TIN(E1-L1,EA-LA,EB-LB)*TIN(E1+L1+1,EA+LA+1,EB+LB+1))
	  RETURN
	  END FUNCTION G

	  REAL*8 FUNCTION TMB(EE,LL,ER,LR,E1,L1,E2,L2,LM,D)
!     	   TALMI-MOSHINSKY BRACKET
!	       (EE,LL;ER,LR:LM/E1,L1;E2,L2:LM)D
	  IMPLICIT NONE
	  REAL*8 S,D,T
	  INTEGER EE,LL,ER,LR,E1,L1,E2,L2,LM
	  INTEGER L,M,ED,LD,EB,LB,EC,LC,EA,LA
 	  TMB=0.D0
	  IF(EE+ER.NE.E1+E2) RETURN
	  IF(TRI(2*LL,2*LR,2*LM)*TRI(2*L1,2*L2,2*LM).EQ.0) RETURN
	  T=DSQRT((D**(E1-ER))/((1.D0+D)**(E1+E2)))
	  M=MIN0(ER,E2)
	  S=1.D0
	  DO 1 ED=0,M
	      EB=ER-ED
	      EC=E2-ED
	      EA=E1-ER+ED
	      DO 2 LD=ED,0,-2 
	          DO 3 LB=EB,0,-2
	               IF(TRI(LD,LB,LR).EQ.0) GO TO 3
	               DO 4 LC=EC,0,-2
	                   IF(TRI(LD,LC,L2).EQ.0) GO TO 4
	                   DO 5 LA=EA,0,-2
	                        IF((TRI(LA,LB,L1).EQ.0).OR.(TRI(LA,LL,LC).EQ.0)) GO TO 5
	                        TMB=TMB+S*T*&
                                C9J(2*LA,2*LB,2*L1,2*LC,2*LD,2*L2,2*LL,2*LR,2*LM)*&
                                G(E1,L1,EA,LA,EB,LB)*G(E2,L2,EC,LC,ED,LD)*&
                                G(EE,LL,EA,LA,EC,LC)*G(ER,LR,EB,LB,ED,LD)
5 	                   CONTINUE
4 	               CONTINUE
3 	          CONTINUE
2 	      CONTINUE
	      S=S*(-D)
1 	  CONTINUE
	  RETURN
	  END FUNCTION TMB
	  
	  SUBROUTINE ORTNTMB(LAST,skirtm)
      IMPLICIT NONE 
      REAL*8 D,skirt,skirtm,atsak,check,tarp1,tarp2
      INTEGER LAST,LMB,EE,EE_,LL,LL_,e,e_,l,l_,l1,l2,e1,e2
      INTEGER*4 skaitl
      skaitl=0
      D=1.D0
      DO LMB=0,LAST
           DO EE=0,LAST
               LL=EE
               DO WHILE(LL.GE.0)
                    DO EE_=0,LAST
                       LL_=EE_
                       DO WHILE(LL_.GE.0)
                           DO e=0,LAST
                                l=e
                                DO WHILE(l.GE.0)
                                    IF(IABS(LL-l).GT.LMB) GOTO 42
                                    IF(LL+l.LT.LMB) GOTO 42
                                    DO e_=0,LAST
                                        l_=e_
                                        DO WHILE(l_.GE.0)
                                            IF(IABS(LL_-l_).GT.LMB) GOTO 43
                                            IF((LL_+l_).LT.LMB) GOTO 43
                                            atsak=0.D0
                                            DO e1=0,EE+e
                                                l1=e1
                                                DO WHILE(l1.GE.0)
                                                    DO e2=0,EE+e
                                                    IF((EE+e).NE.(e1+e2)) CYCLE
                                                    IF((EE_+e_).NE.(e1+e2)) CYCLE
                                                    l2=e2
                                                    DO WHILE(l2.GE.0)
                                                        IF(IABS(l1-l2).GT.LMB) GOTO 41
                                                        IF(l1+l2.LT.LMB) GOTO 41
                                                        tarp1=TMB(EE,LL,e,l,e1,l1,e2,l2,LMB,D)
                                                        tarp2=TMB(e1,l1,e2,l2,EE_,LL_,e_,l_,LMB,D)
                                                        atsak=atsak+tarp1*tarp2
                                                        skaitl=skaitl+2
41                                                      l2=l2-2
                                                    END DO
                                                END DO
                                            l1=l1-2
                                         END DO
                                     END DO
                                     IF((EE.EQ.EE_).AND.(LL.EQ.LL_).AND.&
                                            (e.EQ.e_).AND.(l.EQ.l_)) THEN
                                           check=1.0D0
                                     ELSE
                                           check=0.0D0
                                     END IF
                                     skirt=DABS(atsak-check)
                                     IF(skirtm<skirt) skirtm=skirt
43                                   l_=l_-2
                                  END DO
                               END DO
42                             l=l-2
                           END DO
                       END DO
                       LL_=LL_-2
                   END DO
               END DO
               LL=LL-2
             END DO
        END DO
      END DO
      write(6,4) skaitl
4     FORMAT(' TMB viso =',I10)
      RETURN
      END SUBROUTINE ORTNTMB
      
end module TLMB