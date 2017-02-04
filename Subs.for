C     ###############################################################################################
C
C                                                     AMPLITUDES
C
C     ###############################################################################################
    SUBROUTINE UAMP(
    *     ampName, time, ampValueOld,  dt,  nSvars, svars,
    *     lFlagsInfo,
    *     nSensor, sensorValues, sensorNames, jSensorLookUpTable,
    *     AmpValueNew,
    *     lFlagsDefine,
    *     AmpDerivative, AmpSecDerivative, AmpIncIntegral,
    *     AmpDoubleIntegral)
C
      INCLUDE 'ABA_PARAM.INC'

C     time indices
      parameter (iStepTime        = 1,
     *           iTotalTime       = 2,
     *           nTime            = 2)
C     flags passed in for information
      parameter (iInitialization   = 1,
     *           iRegularInc       = 2,
     *           iCuts             = 3,
     *           ikStep            = 4,
     *           nFlagsInfo        = 4)
C     optional flags to be defined
      parameter (iComputeDeriv       = 1,
     *           iComputeSecDeriv    = 2,
     *           iComputeInteg       = 3,
     *           iComputeDoubleInteg = 4,
     *           iStopAnalysis       = 5,
     *           iConcludeStep       = 6,
     *           nFlagsDefine        = 6)
      dimension time(nTime), lFlagsInfo(nFlagsInfo),
     *          lFlagsDefine(nFlagsDefine)
      dimension jSensorLookUpTable(*)
      dimension sensorValues(nSensor), svars(nSvars)
      character*80 sensorNames(nSensor)
      character*80 ampName

      REAL*8 TOTAL_EPS_X1, TOTAL_EPS_X2, TOTAL_GAMMA
      REAL*8 L, PI, RATIO, X1, X2
      REAL*8 GAMMA, EPS_X1, EPS_X2, VAL
      INTEGER I, MODE

C     1:Axial;   2:Shear;    3:Coupled
      MODE = 1
            
      PI= 2*DACOS(0.0D0)
      L = 10.4D0
      
      IF (MODE .EQ. 1) THEN ! Axial
        TOTAL_EPS_X1 = 4.0D-2
        TOTAL_EPS_X2 = 4.0D-2
        TOTAL_GAMMA  = 0 ! Degree 
      ELSEIF (MODE .EQ. 2) THEN ! Shear
        TOTAL_EPS_X1 = 0.0D0
        TOTAL_EPS_X2 = 0.0D0
        TOTAL_GAMMA  = 40 ! Degree 
      ELSE 
        TOTAL_EPS_X1 = 4.0D-2
        TOTAL_EPS_X2 = 4.0D-2
        TOTAL_GAMMA  = 40 ! Degree 
      ENDIF  
     
      EPS_X1 = time(iStepTime) * TOTAL_EPS_X1 + 1.0D-6
      EPS_X2 = time(iStepTime) * TOTAL_EPS_X2 + 1.0D-6
      GAMMA  = time(iStepTime) * TOTAL_GAMMA * PI / 180 ! To radians
         
      IF (ampName .EQ. 'RP1-X1') THEN
        CALL DEFORMATION_X1(L/2, L/2, VAL,
     2      EPS_X1, EPS_X2, GAMMA)
        AmpValueNew = VAL
        
        WRITE(*,*) "TIME=", time(iStepTime), "  DTIME=", dt
              
      ELSEIF (ampName .EQ. 'RP1-X2') THEN
        CALL DEFORMATION_X2(L/2, L/2, VAL,
     2      EPS_X1, EPS_X2, GAMMA)
        AmpValueNew = VAL
        
      ELSEIF (ampName .EQ. 'RP2-X1') THEN
        CALL DEFORMATION_X1(L/2, -L/2, VAL,
     2     EPS_X1, EPS_X2, GAMMA)
        AmpValueNew = VAL
     
      ELSEIF (ampName .EQ. 'RP2-X2') THEN
        CALL DEFORMATION_X2(L/2, -L/2, VAL,
     2      EPS_X1, EPS_X2, GAMMA)
        AmpValueNew = VAL
     
      ELSEIF (ampName .EQ. 'RP3-X1') THEN
        CALL DEFORMATION_X1(-L/2, L/2, VAL,
     2      EPS_X1, EPS_X2, GAMMA)
        AmpValueNew = VAL
     
      ELSEIF (ampName .EQ. 'RP3-X2') THEN
        CALL DEFORMATION_X2(-L/2, L/2, VAL,
     2      EPS_X1, EPS_X2, GAMMA)
        AmpValueNew = VAL
        
      ELSEIF (ampName .EQ. 'RP4-X1') THEN
        CALL DEFORMATION_X1(-L/2, -L/2, VAL,
     2      EPS_X1, EPS_X2, GAMMA)
        AmpValueNew = VAL
     
      ELSEIF (ampName .EQ. 'RP4-X2') THEN
        CALL DEFORMATION_X2(-L/2, -L/2, VAL,
     2      EPS_X1, EPS_X2, GAMMA)
        AmpValueNew = VAL

      ELSEIF (ampName .EQ. 'XN-NODE-X1') THEN
        CALL DEFORMATION_X1(-L/2, 0, VAL,
     2      EPS_X1, EPS_X2, GAMMA)
        AmpValueNew = VAL  
                 
      ELSEIF (ampName .EQ. 'XN-NODE-X2') THEN
        CALL DEFORMATION_X2(-L/2, 0, VAL,
     2      EPS_X1, EPS_X2, GAMMA)
        AmpValueNew = VAL   

      ELSEIF (ampName .EQ. 'YN-NODE-X1') THEN
        CALL DEFORMATION_X1(0, -L/2, VAL,
     2      EPS_X1, EPS_X2, GAMMA)
        AmpValueNew = VAL  
        
      ELSEIF (ampName .EQ. 'YN-NODE-X2') THEN
        CALL DEFORMATION_X2(0, -L/2, VAL,
     2      EPS_X1, EPS_X2, GAMMA)
        AmpValueNew = VAL           

      ELSE
        WRITE (*,*) "NO AMP!"
        WRITE (*,*) ampName
     
      END IF  
      
      RETURN
      END
C     
C
C     The defined functions for the deformation of the field
C      
C
      SUBROUTINE DEFORMATION_X1(X1_INIT,X2_INIT,X1_NEW,E1,E2,G)
C      X1_NEW is the AmpValueNew
      REAL*8 X1_INIT, X2_INIT, X1_NEW, E1, E2, G
      
      X1_NEW = X1_INIT*(1+E1)*DCOS(G/2) + X2_INIT*(1+E2)*DSIN(G/2)
     2  - X1_INIT
      
      END SUBROUTINE DEFORMATION_X1


    SUBROUTINE DEFORMATION_X2(X1_INIT,X2_INIT,X2_NEW,E1,E2,G)
C      X2_NEW is the AmpValueNew
      REAL*8 X1_INIT, X2_INIT, X2_NEW, E1, E2, G
      
      X2_NEW = X1_INIT*(1+E1)*DSIN(G/2) + X2_INIT*(1+E2)*DCOS(G/2)
     2  - X2_INIT      
      
    END SUBROUTINE DEFORMATION_X2

C
C     ###############################################################################################
C
C                                                     MPC (MULTI-POINT CONSTRAINTS)
C
C     ###############################################################################################
C
    SUBROUTINE MPC(UE,A,JDOF,MDOF,N,JTYPE,X,U,UINIT,MAXDOF,
    * LMPC,KSTEP,KINC,TIME,NT,NF,TEMP,FIELD,LTRAN,TRAN)
C
    INCLUDE 'ABA_PARAM.INC'
C
    DIMENSION A(N),JDOF(N),X(6,N),U(MAXDOF,N),UINIT(MAXDOF,N),
    * TIME(2),TEMP(NT,N),FIELD(NF,NT,N),LTRAN(N),TRAN(3,3,N)

    REAL*8 L, RATIO
    INTEGER I

    L = 10.4

C     ##############################################
C
C     IF USED FOR X FACE SETS
C     U3, NODE-3, xNzP;        U4, NODE-1, xPzP
C     U5, NODE-4, xNzN;        U6, NODE-2, xPzN
C
C     ##############################################
C
C     IF USED FOR Z FACE SETS
C     U3, NODE-1, xPzP;        U4, NODE-2, xPzN
C     U5, NODE-3, xNzP,        U6, NODE-4, xNzN
    IF (JTYPE .EQ. 1) THEN
C     THE CONSTRAINT FOR DOF(1)
    DO I = 1,N
    JDOF(I) = 1
    END DO
    ELSEIF (JTYPE .EQ. 2) THEN
C     THE CONSTRAINT FOR DOF(2)
    DO I = 1,N
    JDOF(I) = 2
    END DO
    ELSE
    WRITE(*,*) "ERROR!"
    END IF

    RATIO = SQRT((X(1,1)-X(1,3))**2 + (X(2,1)-X(2,3))**2)/L

    A(1) =  1
    A(2) = -1
    A(3) = -1 * (1-RATIO)
    A(4) = -1 * RATIO
    A(5) =  1 * (1-RATIO)
    A(6) =  1 * RATIO

    RETURN
    END

C
C     ###############################################################################################
C
C                                                     MATERIAL
C
C     ###############################################################################################
C
    SUBROUTINE UMAT(STRESS,STATEV,DDSDDE,SSE,SPD,SCD,
     1 RPL,DDSDDT,DRPLDE,DRPLDT,
     2 STRAN,DSTRAN,TIME,DTIME,TEMP,DTEMP,PREDEF,DPRED,CMNAME,
     3 NDI,NSHR,NTENS,NSTATV,PROPS,NPROPS,COORDS,DROT_2,PNEWDT,
     4 CELENT,DFGRD0,DFGRD1,NOEL,NPT,LAYER,KSPT,KSTEP,KINC)
C
	INCLUDE 'ABA_PARAM.INC'	
C
      CHARACTER*80 CMNAME
      DIMENSION STRESS(NTENS),STATEV(NSTATV),
     1 DDSDDE(NTENS,NTENS),DDSDDT(NTENS),DRPLDE(NTENS),
     2 STRAN(NTENS),DSTRAN(NTENS),TIME(2),PREDEF(1),DPRED(1),
     3 PROPS(NPROPS),COORDS(3),DROT_2(3,3),DFGRD0(3,3),DFGRD1(3,3)

C
	REAL*8 e0_1(3,1), e0_2(3,1), e0_3(3,1), eb1(3,1), eb2(3,1), eb3(3,1)
	REAL*8 f1(3,1), f2(3,1), f3(3,1), f1v(3,1), f2v(3,1), f2b(1,1)
	REAL*8 MAG_f1v, MAG_f2v, DUM(1,1)
	REAL*8 THETA(3,3), R(3,3), F(3,3)
	REAL*8 USQUARED(3,3), U(3,3), U_Inv(3,3), STATEV_VEC(6,1)
	REAL*8 STRAIN_INC(6,1), NEW_STRESS(6,1), STRAIN_VEC(6,1)
	REAL*8 STRAIN_VEC_F(6,1), NEW_STRESS_F(6,1), STRAIN_INC_F(6,1)
	REAL*8 FABRIC_DDSDDE(6,6), M_MAT(6,6)
	REAL*8 G, EComp, E, E0_lat, E0_22, E0_33, p
	INTEGER ERR

      E     = PROPS(2)
      E0_22 = PROPS(3)
	E0_33 = PROPS(4)
	p	    = PROPS(5)
	
	E0_lat   = 5.0D6
	EComp = E/20
	G  = 25.0D6 * PROPS(6) ! PROPS(7) is The G factor in Isight

	F(1,1) = DBLE(DFGRD1(1,1))
	F(1,2) = DBLE(DFGRD1(1,2))
	F(1,3) = DBLE(DFGRD1(1,3))
	
	F(2,1) = DBLE(DFGRD1(2,1))
	F(2,2) = DBLE(DFGRD1(2,2))
	F(2,3) = DBLE(DFGRD1(2,3))
	
	F(3,1) = DBLE(DFGRD1(3,1))
	F(3,2) = DBLE(DFGRD1(3,2))
	F(3,3) = DBLE(DFGRD1(3,3))

	USQUARED = MATMUL(TRANSPOSE(F),F)

	CALL MATSQRT(USQUARED,U)

	CALL FINDInv(U, U_Inv, 3, ERR)

	R = MATMUL(F, U_Inv)

	F = MATMUL(R,MATMUL(F,TRANSPOSE(R)))

	e0_1 = RESHAPE((/1, 0, 0/), (/3, 1/))
	e0_2 = RESHAPE((/0, 1, 0/), (/3, 1/))
	e0_3 = RESHAPE((/0, 0, 1/), (/3, 1/))
	
	eb1= MATMUL(R, e0_1)
	eb2= MATMUL(R, e0_2)
	eb3= MATMUL(R, e0_3)
	
	f1v     = MATMUL(F ,e0_1)
	MAG_f1v = DSQRT(f1v(1,1)**2 + f1v(2,1)**2 + f1v(3,1)**2)
	f1      = f1v / MAG_f1v

	f2b     = MATMUL(TRANSPOSE(MATMUL(F ,e0_2)),f1)
	f2v     = MATMUL(F, e0_2) - f2b(1,1) * f1
	MAG_f2v = DSQRT(f2v(1,1)**2 + f2v(2,1)**2 + f2v(3,1)**2)
	f2      = f2v / MAG_f2v

	CALL CROSSPRODUCT(f3, f1, f2)

	DUM        = MATMUL(TRANSPOSE(f1),eb1)
	THETA(1,1) = DUM(1,1)
	DUM        = MATMUL(TRANSPOSE(f1),eb2)
	THETA(2,1) = DUM(1,1)
	DUM        = MATMUL(TRANSPOSE(f1),eb3)
	THETA(3,1) = DUM(1,1)

	DUM        = MATMUL(TRANSPOSE(f2),eb1)
	THETA(1,2) = DUM(1,1)
	DUM        = MATMUL(TRANSPOSE(f2),eb2)
	THETA(2,2) = DUM(1,1)
	DUM        = MATMUL(TRANSPOSE(f2),eb3)
	THETA(3,2) = DUM(1,1)

	DUM        = MATMUL(TRANSPOSE(f3),eb1)
	THETA(1,3) = DUM(1,1)
	DUM        = MATMUL(TRANSPOSE(f3),eb2)
	THETA(2,3) = DUM(1,1)
	DUM        = MATMUL(TRANSPOSE(f3),eb3)
	THETA(3,3) = DUM(1,1)

	DO I = 1,6
	DO J = 1,6
		FABRIC_DDSDDE(I,J) = 0.0D0
		M_MAT(I,J) = 0.0D0
	END DO
	END DO
	
	DO I=1,NSTATV
	  STATEV(I) = 0.0D0
	END DO

	DO I = 1,3
		NEW_STRESS(I,1) = DBLE(STRESS(I))
		STRAIN_INC(I,1) = DBLE(DSTRAN(I))
		STRAIN_VEC(I,1)	= DBLE(STRAN(I))
		M_MAT(I,I) = 1.0D0
	J = I + 3		
		NEW_STRESS(J,1) = DBLE(STRESS(J))
		STRAIN_INC(J,1) = DBLE(DSTRAN(J))/2
		STRAIN_VEC(J,1)	= DBLE(STRAN(J))/2
		M_MAT(J,J) = 2.0D0
	END DO

	CALL ROTATE(STRAIN_VEC, STRAIN_VEC_F, THETA)
	
	STR  = DABS(STRAIN_VEC_F(1,1))
      
	IF (STRAIN_VEC_F(1,1) .GE. 0) THEN
		FABRIC_DDSDDE(1,1) = E
	ELSE
		FABRIC_DDSDDE(1,1) = EComp
	END IF
	
	G = G * (1 + STR * 1.0D3)
	
	FABRIC_DDSDDE(2,2) = E0_22 * STR * 
	1  STRAIN_VEC_F(2,1)**2 * DABS(STRAIN_VEC_F(3,1)) + E0_lat 
	FABRIC_DDSDDE(3,3) = E0_33 * STR * 
	1  STRAIN_VEC_F(2,1)**2 * DABS(STRAIN_VEC_F(3,1)) + E0_lat 

	FABRIC_DDSDDE(4,4) = 2*G
	FABRIC_DDSDDE(5,5) = 2*G
	FABRIC_DDSDDE(6,6) = 2*G

	IF (PROPS(1) .EQ. 1.0) THEN ! TSM
		CALL ROTATEC(FABRIC_DDSDDE,DDSDDE,THETA)
		NEW_STRESS = NEW_STRESS + MATMUL(DDSDDE, 
	1		MATMUL(M_MAT,STRAIN_INC))
		CALL ROTATE(NEW_STRESS, STATEV_VEC, THETA)

	ELSE IF (PROPS(1) .EQ. 2.0) THEN ! TSS
		CALL ROTATEC(FABRIC_DDSDDE,DDSDDE,THETA)
		CALL ROTATE(NEW_STRESS, NEW_STRESS_F, THETA)
		CALL ROTATE(STRAIN_INC, STRAIN_INC_F, THETA)

		NEW_STRESS_F = NEW_STRESS_F + MATMUL(FABRIC_DDSDDE,
	1		MATMUL(M_MAT,STRAIN_INC_F))
		STATEV_VEC = NEW_STRESS_F
		CALL ROTATE(NEW_STRESS_F, NEW_STRESS, TRANSPOSE(THETA))

	ELSE IF (PROPS(1) .EQ. 3.0) THEN ! REGULAR STRESS UPDATE
		DDSDDE = FABRIC_DDSDDE
		NEW_STRESS = NEW_STRESS + MATMUL(FABRIC_DDSDDE,
	1		MATMUL(M_MAT,STRAIN_INC))
		STATEV_VEC = NEW_STRESS
	END IF

	DO I = 1,6
		STATEV(I) = STATEV_VEC(I,1)
		STRESS(I) = NEW_STRESS(I,1)
	END DO
      
	STATEV(7) = G
	STATEV(8) = FABRIC_DDSDDE(2,2) 
	STATEV(9) = FABRIC_DDSDDE(3,3)
	
	RETURN
      END

C	###########################################################################
C	THE SUBROUTINES DEVELOPED OR MODIFIED BY MOJTABA
C
C	THE SUBROUTINES WHICH ARE USED WITHOUT ANY CHANGES ARE ADDED 
C	IN THE INCLUDED FILE "EISPACK_Functions.for"

	SUBROUTINE FINDInv(matrix, inverse, n, errorflag)
	!Subroutine to find the inverse of a square matrix
	!Author : Louisda16th a.k.a Ashwith J. Rego
	!Reference : Algorithm has been well explained in:
	!http://math.uww.edu/~mcfarlat/inverse.htm           
	!http://www.tutor.ms.unimelb.edu.au/matrix/matrix_inverse.html

	INTEGER, INTENT(IN) :: n
	INTEGER, INTENT(OUT) :: errorflag  !Return -1 for error, 0 for normal
	REAL*8, INTENT(IN), DIMENSION(n,n) :: matrix  !Input matrix
	REAL*8, INTENT(OUT), DIMENSION(n,n) :: inverse !Inverted matrix
	
	LOGICAL :: FLAG = .TRUE.
	INTEGER :: i, j, k, l
	REAL*8  m
	REAL*8, DIMENSION(n,2*n) :: augmatrix !augmented matrix
	
	!Augment input matrix with an identity matrix
	DO i = 1, n
		DO j = 1, 2*n
			IF (j <= n ) THEN
				augmatrix(i,j) = matrix(i,j)
			ELSE IF ((i+n) == j) THEN
				augmatrix(i,j) = 1
			Else
				augmatrix(i,j) = 0
			ENDIF
		END DO
	END DO
	
	!Reduce augmented matrix to upper traingular form
	DO k =1, n-1
		IF (augmatrix(k,k) == 0) THEN
			FLAG = .FALSE.
			DO i = k+1, n
				IF (augmatrix(i,k) /= 0) THEN
					DO j = 1,2*n
						augmatrix(k,j) = augmatrix(k,j)+augmatrix(i,j)
					END DO
					FLAG = .TRUE.
					EXIT
				ENDIF
				IF (FLAG .EQV. .FALSE.) THEN
					PRINT*, "Matrix is non - invertible"
					inverse = 0
					errorflag = -1
					return
				ENDIF
			END DO
		ENDIF
		DO j = k+1, n			
			m = augmatrix(j,k)/augmatrix(k,k)
			DO i = k, 2*n
				augmatrix(j,i) = augmatrix(j,i) - m*augmatrix(k,i)
			END DO
		END DO
	END DO
	
	!Test for invertibility
	DO i = 1, n
		IF (augmatrix(i,i) == 0) THEN
			PRINT*, "Matrix is non - invertible"
			inverse = 0
			errorflag = -1
			return
		ENDIF
	END DO
	
	!Make diagonal elements as 1
	DO i = 1 , n
		m = augmatrix(i,i)
		DO j = i , (2 * n)				
			   augmatrix(i,j) = (augmatrix(i,j) / m)
		END DO
	END DO
	
	!Reduced right side half of augmented matrix to identity matrix
	DO k = n-1, 1, -1
		DO i =1, k
		m = augmatrix(i,k+1)
			DO j = k, (2*n)
				augmatrix(i,j) = augmatrix(i,j) -augmatrix(k+1,j) * m
			END DO
		END DO
	END DO				
	
	!store answer
	DO i =1, n
		DO j = 1, n
			inverse(i,j) = augmatrix(i,j+n)
		END DO
	END DO
	errorflag = 0
	END SUBROUTINE FINDinv

C	NORMALIZES A VECOTR (EIGENVECTOR)
		
	SUBROUTINE NORMALIZE(INP, OUT)
	REAL*8 INP(3,1), OUT(3,1)
	REAL*8 MAG

	MAG = DSQRT(INP(1,1)**2 + INP(2,1)**2 + INP(3,1)**2)
	OUT = INP / MAG

	END SUBROUTINE NORMALIZE

	SUBROUTINE CROSSPRODUCT(ANS, A, B)
		REAL*8 ANS(3,1), A(3,1), B(3,1)
	  
		ANS(1,1) =  A(2,1)*B(3,1) - A(3,1)*B(2,1)
		ANS(2,1) =-(A(1,1)*B(3,1) - A(3,1)*B(1,1))
		ANS(3,1) =  A(1,1)*B(2,1) - A(2,1)*B(1,1)
	                 
	END SUBROUTINE CROSSPRODUCT

C	SQUARE ROOT OF A MATRIX
C	/////////////////////////////////////////////////////////
	SUBROUTINE MATSQRT(USQUARED, U)
	REAL*8 USQUARED(3,3), TEMP_USQUARED(6), U(3,3), DIAGONALIZED(3,3)
	REAL*8 EIGVAL(3), EIGVECS(3,3), TEMP(3)
	REAL*8 EIGVEC1(3,1), EIGVEC2(3,1), EIGVEC3(3,1), Q(3,3)
	INTEGER ERR, I

	TEMP_USQUARED(1) = USQUARED(1,1)
	TEMP_USQUARED(2) = USQUARED(2,2)
	TEMP_USQUARED(3) = USQUARED(3,3)
	TEMP_USQUARED(4) = USQUARED(1,2)
	TEMP_USQUARED(5) = USQUARED(1,3)
	TEMP_USQUARED(6) = USQUARED(2,3)

	CALL SPRIND(TEMP_USQUARED, EIGVAL, EIGVECS, 1, 3, 3)

	DO I=1,3
		EIGVEC1(I,1) = EIGVECS(1,I)
		EIGVEC2(I,1) = EIGVECS(2,I)
		EIGVEC3(I,1) = EIGVECS(3,I)
	END DO
	
	CALL NORMALIZE(EIGVEC1, EIGVEC1)
	CALL NORMALIZE(EIGVEC2, EIGVEC2)
	CALL NORMALIZE(EIGVEC3, EIGVEC3)

		
	Q(1,1) = EIGVEC1(1,1)
	Q(2,1) = EIGVEC1(2,1)
	Q(3,1) = EIGVEC1(3,1)

	Q(1,2) = EIGVEC2(1,1)
	Q(2,2) = EIGVEC2(2,1)
	Q(3,2) = EIGVEC2(3,1)

	Q(1,3) = EIGVEC3(1,1)
	Q(2,3) = EIGVEC3(2,1)
	Q(3,3) = EIGVEC3(3,1)

	DIAGONALIZED = MATMUL(MATMUL(TRANSPOSE(Q),USQUARED),Q)
	
	DIAGONALIZED = DSQRT(DABS(DIAGONALIZED))
	
	U = MATMUL(MATMUL(Q,DIAGONALIZED),TRANSPOSE(Q))
               
	END SUBROUTINE MATSQRT

C	SUBROUTINE FOR ROTATING DDSDDE MATRIX
C	/////////////////////////////////////////
	SUBROUTINE ROTATEC(C,C_PRIME, TH)
	REAL*8 C(6,6), C_PRIME(6,6), T(6,6), TP(6,6), TH(3,3)
	REAL*8 M(6,6), M_INV(6,6)
	REAL*8 F11, F12, F13
	REAL*8 F21, F22, F23
	REAL*8 F31, F32, F33
	INTEGER I, J

	F11 = TH(1,1)
	F12 = TH(2,1)
	F13 = TH(3,1)
	F21 = TH(1,2)
	F22 = TH(2,2)
	F23 = TH(3,2)
	F31 = TH(1,3)
	F32 = TH(2,3)
	F33 = TH(3,3)

	DO I = 1,6
	DO J = 1,6
		C_PRIME(I,J) = 0.0D0
		M(I,J)		 = 0.0D0
		M_INV(I,J)	 = 0.0D0
	END DO
	END DO

C	DEFINING T MATRIX

	T(1,1) = F11**2
	T(1,2) = F21**2
	T(1,3) = F31**2
	T(2,1) = F12**2
	T(2,2) = F22**2
	T(2,3) = F32**2
	T(3,1) = F13**2
	T(3,2) = F23**2
	T(3,3) = F33**2

	T(1,4) = 2 * F11 * F21
	T(1,5) = 2 * F11 * F31
	T(1,6) = 2 * F31 * F21
		
	T(2,4) = 2 * F12 * F22
	T(2,5) = 2 * F12 * F32
	T(2,6) = 2 * F32 * F22

	T(3,4) = 2 * F13 * F23
	T(3,5) = 2 * F13 * F33
	T(3,6) = 2 * F33 * F23	
	
	T(4,1) = F11 * F12
	T(4,2) = F21 * F22
	T(4,3) = F31 * F32	

	T(5,1) = F11 * F13
	T(5,2) = F21 * F23
	T(5,3) = F31 * F33

	T(6,1) = F12 * F13
	T(6,2) = F22 * F23
	T(6,3) = F32 * F33
	
	T(4,4) = F12 * F21 + F11 * F22
	T(4,5) = F12 * F31 + F11 * F32
	T(4,6) = F22 * F31 + F21 * F32
	
	T(5,4) = F13 * F21 + F11 * F23
	T(5,5) = F13 * F31 + F11 * F33
	T(5,6) = F23 * F31 + F21 * F33	

	T(6,4) = F13 * F22 + F12 * F23
	T(6,5) = F13 * F32 + F12 * F33
	T(6,6) = F23 * F32 + F22 * F33	

C	DEFINING THE T_PRIME MATRIX
	
	TP(1,1) = F11**2
	TP(1,2) = F12**2
	TP(1,3) = F13**2
	TP(2,1) = F21**2
	TP(2,2) = F22**2
	TP(2,3) = F23**2
	TP(3,1) = F31**2
	TP(3,2) = F32**2
	TP(3,3) = F33**2

	TP(1,4) = 2 * F11 * F12
	TP(1,5) = 2 * F11 * F13
	TP(1,6) = 2 * F12 * F13

	TP(2,4) = 2 * F21 * F22
	TP(2,5) = 2 * F21 * F23
	TP(2,6) = 2 * F22 * F23

	TP(3,4) = 2 * F31 * F32
	TP(3,5) = 2 * F31 * F33
	TP(3,6) = 2 * F32 * F33	

	TP(4,1) = F11 * F21
	TP(4,2) = F12 * F22
	TP(4,3) = F13 * F23

	TP(5,1) = F11 * F31
	TP(5,2) = F12 * F32
	TP(5,3) = F13 * F33

	TP(6,1) = F21 * F31
	TP(6,2) = F22 * F32
	TP(6,3) = F23 * F33

	TP(4,4) = F12 * F21 + F11 * F22
	TP(4,5) = F13 * F21 + F11 * F23
	TP(4,6) = F13 * F22 + F12 * F23

	TP(5,4) = F12 * F31 + F11 * F32
	TP(5,5) = F13 * F31 + F11 * F33
	TP(5,6) = F13 * F32 + F12 * F33

	TP(6,4) = F22 * F31 + F21 * F32
	TP(6,5) = F23 * F31 + F21 * F33
	TP(6,6) = F23 * F32 + F22 * F33

	DO I = 1,3
		M(I,I)     = 1.0D0
		M_INV(I,I) = 1.0D0
	J = I + 3
		M(J,J)     = 2.0D0
		M_INV(J,J) = 0.5D0
	END DO	
	
	C_PRIME = MATMUL(T,MATMUL(C,MATMUL(M,MATMUL(TP,M_INV))))

	END SUBROUTINE ROTATEC


	SUBROUTINE ROTATE(S,S_PRIME, TH)
	REAL*8 S(6,1), S_PRIME(6,1), TP(6,6), TH(3,3)
	REAL*8 F11, F12, F13
	REAL*8 F21, F22, F23
	REAL*8 F31, F32, F33
	INTEGER I, J

	F11 = TH(1,1)
	F12 = TH(2,1)
	F13 = TH(3,1)
	F21 = TH(1,2)
	F22 = TH(2,2)
	F23 = TH(3,2)
	F31 = TH(1,3)
	F32 = TH(2,3)
	F33 = TH(3,3)

	TP(1,1) = F11**2
	TP(1,2) = F12**2
	TP(1,3) = F13**2
	TP(2,1) = F21**2
	TP(2,2) = F22**2
	TP(2,3) = F23**2
	TP(3,1) = F31**2
	TP(3,2) = F32**2
	TP(3,3) = F33**2

	TP(1,4) = 2 * F11 * F12
	TP(1,5) = 2 * F11 * F13
	TP(1,6) = 2 * F12 * F13

	TP(2,4) = 2 * F21 * F22
	TP(2,5) = 2 * F21 * F23
	TP(2,6) = 2 * F22 * F23

	TP(3,4) = 2 * F31 * F32
	TP(3,5) = 2 * F31 * F33
	TP(3,6) = 2 * F32 * F33	

	TP(4,1) = F11 * F21
	TP(4,2) = F12 * F22
	TP(4,3) = F13 * F23

	TP(5,1) = F11 * F31
	TP(5,2) = F12 * F32
	TP(5,3) = F13 * F33

	TP(6,1) = F21 * F31
	TP(6,2) = F22 * F32
	TP(6,3) = F23 * F33

	TP(4,4) = F12 * F21 + F11 * F22
	TP(4,5) = F13 * F21 + F11 * F23
	TP(4,6) = F13 * F22 + F12 * F23

	TP(5,4) = F12 * F31 + F11 * F32
	TP(5,5) = F13 * F31 + F11 * F33
	TP(5,6) = F13 * F32 + F12 * F33

	TP(6,4) = F22 * F31 + F21 * F32
	TP(6,5) = F23 * F31 + F21 * F33
	TP(6,6) = F23 * F32 + F22 * F33

	S_PRIME = MATMUL(TP,S)

	END SUBROUTINE ROTATE     
	
