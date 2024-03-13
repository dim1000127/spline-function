program spline_function
   use Environment
   
   implicit none
   character(*), parameter    :: output_file = "output.txt"
   integer                    :: Out = 0, i = 0, x_min = 1, x_max = 4, length = 0, NOFUN = 0
   real(R_)                   :: h = 0.375, lower_limit = 0, upper_limit = 20
   real(R_)                   :: RELERR = 0.1, ABSERR = 0.1, ERREST, FLAG, RESULT, X_f
   real (R_), allocatable     :: F(:), X(:), B(:), C(:), D(:), X2(:)
   real (R_), allocatable     :: SevalF(:), LagrangeF(:), QUNC8_2(:)

   length = (x_max - x_min) / h + 1

   allocate(F(length), X(length), B(length), C(length), D(length))
   allocate(X2(length), SevalF(length), LagrangeF(length), QUNC8_2(length))

   X(1) = x_min;
   DO i = 2, length
      X(i) = X(i-1) + h 
   END DO

   DO i = 1, length
     X2(i) = 1.1875 + 0.375 * (i-1)
   END DO

   open (file=output_file, encoding=E_, newunit=Out)
      write (Out, *) "1 практическое задание"
      write (Out, *)
      write (Out, *) "Расчет интеграла с помощью QUANC8"
   close (Out)

   DO i = 1, length
      X_f = X(i)
      call QUANC8(FUNC, lower_limit, upper_limit, ABSERR, RELERR, RESULT, ERREST, NOFUN, FLAG)
      F(i) = RESULT
      open (file=output_file, encoding=E_, newunit=Out, position="append")
         write (Out, '(a, T4, "= ", f9.6)') "X", X(i)
         write (Out, '(a, T4, "= ", f9.6)') "F", RESULT
         write (Out, *)
      close (Out)
   END DO

   call SPLINE(length, X, F, B, C, D)

   open (file=output_file, encoding=E_, newunit=Out, position="append")
      write (Out, *)
      write (Out, *) "SEVAL"
   close (Out)

   DO i = 1, length
      SevalF(i) = SEVAL(length, X2(i), X, F, B, C, D)
      open (file=output_file, encoding=E_, newunit=Out, position="append")
         write (Out, '(a, T4, "= ", f9.6)') "X", X2(i)
         write (Out, '(a, T4, "= ", f9.6)') "Seval", SevalF(i)
         write (Out, *)
      close (Out)
   END DO

   open (file=output_file, encoding=E_, newunit=Out, position="append")
      write (Out, *)
      write (Out, *) "LAGRANGE"
   close (Out)

   DO i = 1, length
      LagrangeF(i) = LAGRANGE(X2(i), X, F)
      open (file=output_file, encoding=E_, newunit=Out, position="append")
         write (Out, '(a, T4, "= ", f9.6)') "X", X2(i)
         write (Out, '(a, T4, "= ", f9.6)') "Lagrange", LagrangeF(i)
         write (Out, *)
      close (Out)
   END DO

   open (file=output_file, encoding=E_, newunit=Out, position="append")
      write (Out, *)
      write (Out, *) "QUNC8 со вторым набором точек"
   close (Out)

   DO i = 1, length
      X_f = X2(i)
      call QUANC8(FUNC, lower_limit, upper_limit, ABSERR, RELERR, RESULT, ERREST, NOFUN, FLAG)
      QUNC8_2(i) = RESULT
      open (file=output_file, encoding=E_, newunit=Out, position="append")
         write (Out, '(a, T4, "= ", f9.6)') "X", X2(i)
         write (Out, '(a, T4, "= ", f9.6)') "F", RESULT
         write (Out, *)
      close (Out)
   END DO

contains
   !Вычисление подинтегральной функции
   REAL FUNCTION FUNC(z)
      real(R_), intent(in)  :: z
      
      FUNC = 1 / (exp(z)*(z+X_f))
      
   END FUNCTION FUNC

   REAL FUNCTION LAGRANGE(x0, x, y)
      real(R_), intent(in)   :: x0, x(:), y(:)
      real (R_)              :: temp
      integer                :: i, j
      
      LAGRANGE = 0
      DO i = 1, length
         temp = 1
         DO j = 1, length
            if (i /= j) then
               temp = temp * (x0 - x(j)) / (x(i) - x(j))            
            end if
         END DO 
         LAGRANGE = LAGRANGE + y(i) * temp
      END DO
   END FUNCTION LAGRANGE

   SUBROUTINE QUANC8(FUN,A,B,ABSERR,RELERR,RESULT, ERREST,NOFUN,FLAG)
       REAL FUN,A,B,ABSERR,RELERR,RESULT,ERREST,FLAG
       INTEGER NOFUN

       REAL W0,W1,W2,W3,W4,AREA,X0,F0,STONE,STEP,COR11,TEMP
       REAL QPREV,QNOW,QDIFF,QLEFT,ESTERR,TOLERR
       REAL QRIGHT(31),F(16),X(16),FSAVE(8,30),XSAVE(8,30)

       INTEGER LEVMIN,LEVMAX,LEVOUT,NOMAX,NOFIN,LEV,NIM,I,J
       LEVMIN=1
       LEVMAX=30
       LEVOUT=6
       NOMAX=5000
       NOFIN=NOMAX-8*(LEVMAX-LEVOUT+2**(LEVOUT+1))

       W0=3956.0/14175.0
       W1=23552.0/14175.0
       W2=-3712.0/14175.0
       W3=41984.0/14175.0
       W4=-18160.0/14175.0

       FLAG=0.0
       RESULT=0.0
       COR11=0.0
       ERREST=0.0
       AREA=0.0
       NOFUN=0
       IF(A.EQ.B)RETURN

       LEV=0
       NIM=1
       X0=A
       X(16)=B
       QPREV=0.0
       F0=FUN(X0)
       STONE=(B-A)/16.0
       X(8)=(X0+X(16))/2.0
       X(4)=(X0+X(8))/2.0
       X(12)=(X(8)+X(16))/2.0
       X(2)=(X0+X(4))/2.0
       X(6)=(X(4)+X(8))/2.0
       X(10)=(X(8)+X(12))/2.0
       X(14)=(X(12)+X(16))/2.0
       DO 25 J=2,16,2
       F(J)=FUN(X(J))
    25 CONTINUE
       NOFUN=9

    30 X(1)=(X0+X(2))/2.0
       F(1)=FUN(X(1))
       DO 35 J=3,15,2
       X(J)=(X(J-1)+X(J+1))/2.0
       F(J)=FUN(X(J))
    35 CONTINUE
       NOFUN=NOFUN+8
       STEP=(X(16)-X0)/16.0
       QLEFT=(W0*(F0+F(8))+W1*(F(1)+F(7))+W2*(F(2)+F(6))+W3*(F(3)+F(5))+W4*F(4))*STEP
       QRIGHT(LEV+1)=(W0*(F(8)+F(16))+W1*(F(9)+F(15))+W2*(F(10)+F(14))+W3*(F(11)+F(13))+W4*F(12))*STEP
       QNOW=QLEFT+QRIGHT(LEV+1)
       QDIFF=QNOW-QPREV
       AREA=AREA+QDIFF

       ESTERR=ABS(QDIFF)/1023.0
       TOLERR=AMAX1(ABSERR,RELERR*ABS(AREA))*(STEP/STONE)
       IF(LEV.LT.LEVMIN)GO TO 50
       IF(LEV.GE.LEVMAX)GO TO 62
       IF(NOFUN.GT.NOFIN)GO TO 60
       IF(ESTERR.LE.TOLERR)GO TO 70

    50 NIM=2*NIM
       LEV=LEV+1

       DO 52 I=1,8
       FSAVE(I,LEV)=F(I+8)
       XSAVE(I,LEV)=X(I+8)
    52 CONTINUE

       QPREV=QLEFT
       DO 55 I=1,8
       J=-I
       F(2*J+18)=F(J+9)
       X(2*J+18)=X(J+9)
    55 CONTINUE
       GO TO 30

    60 NOFIN=2*NOFIN
       LEVMAX=LEVOUT
       FLAG=FLAG+(B-X0)/(B-A)
       GO TO 70

    62 FLAG=FLAG+1.0

    70 RESULT=RESULT+QNOW
       ERREST=ERREST+ESTERR
       COR11=COR11+QDIFF/1023.0

    72 IF(NIM.EQ.2*(NIM/2))GO TO 75
       NIM=NIM/2
       LEV=LEV-1
       GO TO 72
    75 NIM=NIM+1
       IF(LEV.LE.0)GO TO 80

       QPREV=QRIGHT(LEV)
       X0=X(16)
       F0=F(16)
       DO 78 I=1,8
       F(2*I)=FSAVE(I,LEV)
       X(2*I)=XSAVE(I,LEV)
    78 CONTINUE
       GO TO 30

    80 RESULT=RESULT+COR11

       IF(ERREST.EQ.0.0)RETURN
    82 TEMP=ABS(RESULT)+ERREST
       IF(TEMP.NE.ABS(RESULT))RETURN
       ERREST=2.0*ERREST
       GO TO 82
   END SUBROUTINE QUANC8

   REAL FUNCTION SEVAL(N,U,X,Y,B,C,D)
       INTEGER N
       REAL U,X(N),Y(N),B(N),C(N),D(N)

       INTEGER I,J,K
       REAL DX
       DATA I/1/
       IF(I.GE.N) I=1
       IF(U.LT.X(I)) GO TO 10
       IF(U.LE.X(I+1)) GO TO 30

  10   I=1
       J=N+1
  20   K=(I+J)/2
       IF(U.LT.X(K))J=K
       IF(U.GE.X(K))I=K
       IF(J.GT.I+1)GO TO 20
 
  30   DX=U-X(I)
       SEVAL=Y(I)+DX*(B(I)+DX*(C(I)+DX*D(I)))
       RETURN
   END FUNCTION SEVAL
 
   SUBROUTINE SPLINE(N,X,Y,B,C,D)
        INTEGER N
        REAL X(N),Y(N),B(N),C(N),D(N)
        INTEGER NM1,IB,I
        REAL T

        NM1=N-1
        IF(N.LT.2) RETURN
        IF(N.LT.3) GO TO 50

        D(1)=X(2)-X(1)
        C(2)=(Y(2)-Y(1))/D(1)
        DO 10 I=2,NM1
           D(I)=X(I+1)-X(I)
           B(I)=2.*(D(I-1)+D(I))
           C(I+1)=(Y(I+1)-Y(I))/D(I)
           C(I)=C(I+1)-C(I)
 10     CONTINUE

        B(1)=-D(1)
        B(N)=-D(N-1)
        C(1)=0.
        C(N)=0.
        IF(N.EQ.3) GO TO 15
        C(1)=C(3)/(X(4)-X(2))-C(2)/(X(3)-X(1))
        C(N)=C(N-1)/(X(N)-X(N-2))-C(N-2)/(X(N-1)-X(N-3))
        C(1)=C(1)*D(1)**2/(X(4)-X(1))
        C(N)=-C(N)*D(N-1)**2/(X(N)-X(N-3))

 15     DO 20 I=2,N
           T=D(I-1)/B(I-1)
           B(I)=B(I)-T*D(I-1)
           C(I)=C(I)-T*C(I-1)
 20     CONTINUE

        C(N)=C(N)/B(N)
        DO 30 IB=1,NM1
           I=N-IB
           C(I)=(C(I)-D(I)*C(I+1))/B(I)
 30     CONTINUE

        B(N)=(Y(N)-Y(NM1))/D(NM1)+D(NM1)*(C(NM1)+2.*C(N))
        DO 40 I=1,NM1
           B(I)=(Y(I+1)-Y(I))/D(I)-D(I)*(C(I+1)+2.*C(I))
           D(I)=(C(I+1)-C(I))/D(I)
           C(I)=3.*C(I)
 40     CONTINUE
         C(N)=3.*C(N)
         D(N)=D(N-1)
         RETURN

 50     B(1)=(Y(2)-Y(1))/(X(2)-X(1))
        C(1)=0.
        D(1)=0.
        B(2)=B(1)
        C(2)=0.
        D(2)=0.
        RETURN
      END SUBROUTINE SPLINE
end program spline_function
