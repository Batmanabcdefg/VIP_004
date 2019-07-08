C     PROGRAM MESHGEN
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      DIMENSION XROW(500), YROW(500), X1(5000), X2(5000), LM1(4)
     1 ,LM2(4)
C
C---- GENERATES MESH FOR Rectangular domain
C
      open(unit=5,file='meshgen.inp',form='formatted',access=
     1 'sequential')
      open(unit=6,file='mesh.dat',form='formatted',
     1 access='sequential')
c---- input # nodes along x (NROW) and y directions (NCOL)
      READ(5,*) NROW,NCOL,xmax,ymax
c---- set up # elems in a column (y-dirn) NELCOL & in a row (x-dirn) NELROW
      NELCOL = NCOL - 1                                                         !nelcol is 40, ncol is 41, y direction
      NELROW = NROW - 1                                                         !nelrow is 20, nrow is 21, x direction
      NUMNP = NROW*NCOL
      NUMEL = NELROW*NELCOL
      delx = xmax/dfloat(nelrow)                                                 !x length of each element
      dely = ymax/dfloat(nelcol)                                                 !y length of each element
      WRITE(6,2000) NUMNP,NUMEL
c     READ(5,1001) (XROW(I), I=1,NROW) 
c     READ(5,1001) (YROW(I), I=1,NCOL) 
      II = 0
c---- number nodes by first going x dirn (holding y fixed) & then y
      DO 100 J = 1,NCOL                                                         !y direction
        DO 200 I = 1,NROW                                                       !x direction
          II = II+1
          X1(II) = dfloat(I-1)*delx + 0.d0                                        !x coordinate of successive nodes in x direction
          X2(II) = dfloat(J-1)*dely + 0.d0                                        !y coordinate of successive nodes in y direction
200     CONTINUE
100   CONTINUE
      LINC = 0
      DO 300 I =1,NUMNP
        WRITE(6,2001) I,LINC,X1(I),X2(I)                                        !write the node number, x coordinate and y coordinate
300   CONTINUE
      LM1(1) = 1                                                                !LM1(1) is 1, LM1(2) is 2, LM1(3) is 22, LM1(4) is 23...
      LM1(2) = LM1(1) + 1
      LM1(4) = LM1(1) + NROW 
      LM1(3) = LM1(4) + 1
      DO 400 I=1,4
C400   LM2(I) = LM1(I) +(NYROW-2)*NXROW
400   LM2(I) = LM1(I) +NELROW-1                                                 !LM2(1) is 20, LM2(2) is 21, LM2(3) is 41, LM2(4) is 42...
      MAT = 1
      LX1 = 1
      LX2 = 0
      IEL1 = 1
      IEL2 = IEL1+NELROW-1                                                      !IEL2 is 20...
      DO 500 I = 1,NELCOL
        WRITE(6,2002) IEL1,MAT,(LM1(K),K=1,4),LX1
        WRITE(6,2002) IEL2,MAT,(LM2(K),K=1,4),LX2
        IF(I.EQ.NELCOL) GO TO 500
        IEL1 = IEL1+NELROW                                                      !IEL1 is 1,21,41,...
        IEL2 = IEL2+NELROW                                                      !IEL2 is 20,40,60,...
        DO 600 K=1,4
          LM1(K) = LM1(K)+NROW
          LM2(K) = LM2(K)+NROW
600     CONTINUE
500   CONTINUE
C
      STOP
C
1000  FORMAT(2I5)
1001  FORMAT(7E10.0)
1002  FORMAT(E10.0)
2000  FORMAT(2I5)
2001  FORMAT(2I5,2E10.3)
2002  FORMAT(7I5)
2003  FORMAT(3I5)
2004  FORMAT(E10.3)
C
      END 
C

