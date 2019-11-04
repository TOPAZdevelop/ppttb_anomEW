      SUBROUTINE SUUB_TTB(P1,ANS)
C  
C FUNCTION GENERATED BY MADGRAPH
C RETURNS AMPLITUDE SQUARED SUMMED/AVG OVER COLORS
C AND HELICITIES
C FOR THE POINT IN PHASE SPACE P(0:3,NEXTERNAL)
C  
C FOR PROCESS : u u~ -> t t~  
C  
C Crossing   1 is u u~ -> t t~  
      IMPLICIT NONE
C  
C CONSTANTS
C  
      INTEGER    NEXTERNAL,   NCOMB,     NCROSS         
      PARAMETER (NEXTERNAL=4, NCOMB= 16, NCROSS=  1)
      INTEGER    THEL
      PARAMETER (THEL=NCOMB*NCROSS)
C  
C ARGUMENTS 
C  
      REAL*8 P1(0:3,NEXTERNAL),ANS(NCROSS)
C  
C LOCAL VARIABLES 
C  
      INTEGER NHEL(NEXTERNAL,NCOMB),NTRY
      REAL*8 T, P(0:3,NEXTERNAL)
      REAL*8 UUB_TTB
      INTEGER IHEL,IDEN(NCROSS),IC(NEXTERNAL,NCROSS)
      INTEGER IPROC,JC(NEXTERNAL)
      LOGICAL GOODHEL(NCOMB,NCROSS)
      DATA GOODHEL/THEL*.FALSE./
      DATA NTRY/0/
      DATA (NHEL(IHEL,  1),IHEL=1,4) / -1, -1, -1, -1/
      DATA (NHEL(IHEL,  2),IHEL=1,4) / -1, -1, -1,  1/
      DATA (NHEL(IHEL,  3),IHEL=1,4) / -1, -1,  1, -1/
      DATA (NHEL(IHEL,  4),IHEL=1,4) / -1, -1,  1,  1/
      DATA (NHEL(IHEL,  5),IHEL=1,4) / -1,  1, -1, -1/
      DATA (NHEL(IHEL,  6),IHEL=1,4) / -1,  1, -1,  1/
      DATA (NHEL(IHEL,  7),IHEL=1,4) / -1,  1,  1, -1/
      DATA (NHEL(IHEL,  8),IHEL=1,4) / -1,  1,  1,  1/
      DATA (NHEL(IHEL,  9),IHEL=1,4) /  1, -1, -1, -1/
      DATA (NHEL(IHEL, 10),IHEL=1,4) /  1, -1, -1,  1/
      DATA (NHEL(IHEL, 11),IHEL=1,4) /  1, -1,  1, -1/
      DATA (NHEL(IHEL, 12),IHEL=1,4) /  1, -1,  1,  1/
      DATA (NHEL(IHEL, 13),IHEL=1,4) /  1,  1, -1, -1/
      DATA (NHEL(IHEL, 14),IHEL=1,4) /  1,  1, -1,  1/
      DATA (NHEL(IHEL, 15),IHEL=1,4) /  1,  1,  1, -1/
      DATA (NHEL(IHEL, 16),IHEL=1,4) /  1,  1,  1,  1/
      DATA (  IC(IHEL,  1),IHEL=1,4) /  1,  2,  3,  4/
      DATA (IDEN(IHEL),IHEL=  1,  1) /  36/
C ----------
C BEGIN CODE
C ----------
      NTRY=NTRY+1
      DO IPROC=1,NCROSS
      CALL SWITCHMOM(P1,P,IC(1,IPROC),JC,NEXTERNAL)
      DO IHEL=1,NEXTERNAL
         JC(IHEL) = +1
      ENDDO
      ANS(IPROC) = 0D0
      DO IHEL=1,NCOMB
          IF (GOODHEL(IHEL,IPROC) .OR. NTRY .LT. 2) THEN
             T=UUB_TTB(P ,NHEL(1,IHEL),JC(1))            
             ANS(IPROC)=ANS(IPROC)+T
              IF (T .GT. 0D0 .AND. .NOT. GOODHEL(IHEL,IPROC)) THEN
                  GOODHEL(IHEL,IPROC)=.TRUE.
C             WRITE(*,*) IHEL,T
              ENDIF
          ENDIF
      ENDDO
      ANS(IPROC)=ANS(IPROC)/DBLE(IDEN(IPROC))
      ENDDO
      END
       
       
      REAL*8 FUNCTION UUB_TTB(P,NHEL,IC)
C  
C FUNCTION GENERATED BY MADGRAPH
C RETURNS AMPLITUDE SQUARED SUMMED/AVG OVER COLORS
C FOR THE POINT WITH EXTERNAL LINES W(0:6,NEXTERNAL)
C  
C FOR PROCESS : u u~ -> t t~  
C  
      IMPLICIT NONE
C  
C CONSTANTS
C  
      INTEGER    NGRAPHS,    NEIGEN,    NEXTERNAL       
      PARAMETER (NGRAPHS=   2,NEIGEN=  1,NEXTERNAL=4)   
      INTEGER    NWAVEFUNCS     , NCOLOR
      PARAMETER (NWAVEFUNCS=   6, NCOLOR=   1) 
      REAL*8     ZERO
      PARAMETER (ZERO=0D0)
C  
C ARGUMENTS 
C  
      REAL*8 P(0:3,NEXTERNAL)
      INTEGER NHEL(NEXTERNAL), IC(NEXTERNAL)
C  
C LOCAL VARIABLES 
C  
      INTEGER I,J
      COMPLEX*16 ZTEMP
      REAL*8 DENOM(NCOLOR), CF(NCOLOR,NCOLOR)
      COMPLEX*16 AMP(NGRAPHS), JAMP(NCOLOR)
      COMPLEX*16 W(6,NWAVEFUNCS)
C  
C GLOBAL VARIABLES
C  
      INCLUDE "coupl.inc"
C  
C COLOR DATA
C  
      DATA Denom(1  )/            1/                                       
      DATA (CF(i,1  ),i=1  ,1  ) /     9/                                  
C               T[2,1]T[3,4]                                               
C ----------
C BEGIN CODE
C ----------
      CALL IXXXXX(P(0,1   ),ZERO ,NHEL(1   ),+1*IC(1   ),W(1,1   ))        
      CALL OXXXXX(P(0,2   ),ZERO ,NHEL(2   ),-1*IC(2   ),W(1,2   ))        
      CALL OXXXXX(P(0,3   ),TMASS ,NHEL(3   ),+1*IC(3   ),W(1,3   ))       
      CALL IXXXXX(P(0,4   ),TMASS ,NHEL(4   ),-1*IC(4   ),W(1,4   ))       
      CALL JIOXXX(W(1,1   ),W(1,2   ),GAU ,ZERO    ,ZERO    ,W(1,5   ))    
      CALL IOVXXX(W(1,4   ),W(1,3   ),W(1,5   ),GAU ,AMP(1   ))            
      CALL JIOXXX(W(1,1   ),W(1,2   ),GZU ,ZMASS   ,ZWIDTH  ,W(1,6   ))    
      CALL IOVXXX(W(1,4   ),W(1,3   ),W(1,6   ),GZU ,AMP(2   ))            
      JAMP(   1) = -AMP(   1) -AMP(   2) 
      UUB_TTB = 0.D0 
      DO I = 1, NCOLOR
          ZTEMP = (0.D0,0.D0)
          DO J = 1, NCOLOR
              ZTEMP = ZTEMP + CF(J,I)*JAMP(J)
          ENDDO
          UUB_TTB =UUB_TTB+ZTEMP*DCONJG(JAMP(I))/DENOM(I)   
      ENDDO
C      CALL GAUGECHECK(JAMP,ZTEMP,EIGEN_VEC,EIGEN_VAL,NCOLOR,NEIGEN) 
      END
