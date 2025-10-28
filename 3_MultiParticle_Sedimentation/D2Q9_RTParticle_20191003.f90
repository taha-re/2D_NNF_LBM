! By :: Taha Rezaee (rezaee.taha@gmail.com)
MODULE PARAM
	DOUBLE PRECISION,PARAMETER :: PI=3.14159265359
	LOGICAL,PARAMETER :: NON_NEWTONIAN= .FALSE.
	LOGICAL,PARAMETER :: READ_INPUT= .FALSE.
	LOGICAL,PARAMETER :: WRITE_OUTPUT= .FALSE.

    INTEGER,PARAMETER :: NX=512 !Number of Cells in X-direction (nodes number 0 to nx)
    INTEGER,PARAMETER :: NY=512 !Number of Cells in Y-direction (nodes number 0 to ny)
    INTEGER,PARAMETER :: Q=9    !Number of discrete microscopic velocity
    DOUBLE PRECISION,PARAMETER :: W(0:Q-1)=(/4.0/9.0,1.0/9.0,1.0/9.0,1.0/9.0,1.0/9.0,&
                                &1.0/36.0,1.0/36.0,1.0/36.0,1.0/36.0/) !Weighting function for Feq
    INTEGER,PARAMETER :: CX(0:Q-1)=(/0,1,0,-1,0,1,-1,-1,1/) !X-component of Discrete Velocity
    INTEGER,PARAMETER :: CY(0:Q-1)=(/0,0,1,0,-1,1,1,-1,-1/) !Y-component of Discrete Velocity
	
	!DOUBLE PRECISION,PARAMETER :: RE_S=0.78
	!DOUBLE PRECISION,PARAMETER :: BN_S=0.033
	DOUBLE PRECISION,PARAMETER :: RE_T=1.0
	DOUBLE PRECISION,PARAMETER :: BN_T=1.0
	DOUBLE PRECISION,PARAMETER :: Y_g=0.0	!Lattice Units (Yeild Number or Bingham Number
	DOUBLE PRECISION,PARAMETER :: RHO_S=1.01
	DOUBLE PRECISION,PARAMETER :: RHO_F=1.0 !Lattice Units
	DOUBLE PRECISION,PARAMETER :: NU=0.4915/3.0	!Lattice Units (Kinematic Viscosity of La
	DOUBLE PRECISION,PARAMETER :: NU_PHYSICAL=0.01 !cm^2/sec
	DOUBLE PRECISION,PARAMETER :: DIAM_PHYSICAL=0.0625 !cm
	DOUBLE PRECISION,PARAMETER :: L_C=16.0	!Lattice Units Characteristic Length of the problem
	DOUBLE PRECISION,PARAMETER :: Omg_PHYSICAL=0.0 !rad/s
	!DOUBLE PRECISION,PARAMETER :: U_C=1.0 !RE_T*NU/L_C   ! Characteristic Velocity of Problem
	!DOUBLE PRECISION,PARAMETER :: A_G=-RE_S*NU*NU/((L_C**3)*(RHO_S-RHO_F)*RHO_F)

	DOUBLE PRECISION,PARAMETER :: C_L=DIAM_PHYSICAL/L_C !cm
	DOUBLE PRECISION,PARAMETER :: C_T=(C_L**2)*NU/NU_PHYSICAL !sec
	DOUBLE PRECISION,PARAMETER :: C_U=C_L/C_T   !cm/sec
	DOUBLE PRECISION,PARAMETER :: C_G=C_L*0.01/(C_T**2) !m^2/sec
	
	DOUBLE PRECISION,PARAMETER :: A_C=1.5*L_C !Amplitude in LU
	DOUBLE PRECISION,PARAMETER :: OMG_OSC=OMG_PHYSICAL*C_T !Angular Frequency of Oscillation=Max Vel LB/(2*pi*amplitude)
	!DOUBLE PRECISION,PARAMETER :: PERIOD_OSC=2.0*PI/OMG_OSC
	DOUBLE PRECISION,PARAMETER :: U_C=A_C*OMG_OSC
	INTEGER,PARAMETER :: TIME_MAX=NINT(24.0/C_T) !NINT(5.01*PERIOD_OSC)
	INTEGER,PARAMETER :: TIME_SAMPLE=NINT(1.0/C_T) !NINT(0.25*PERIOD_OSC)
	! ! 1 SHOULD BE CHANGED FOR EACH CASE
	!DOUBLE PRECISION,PARAMETER :: T_C=1.0/0.142784
	!DOUBLE PRECISION,PARAMETER :: U_C=1.0/1.400714
	!DOUBLE PRECISION,PARAMETER :: G_C=2.0*PI/T_C
	!DOUBLE PRECISION,PARAMETER :: F_C=Nu*U_C
	DOUBLE PRECISION,PARAMETER :: A_G=-9.81/C_G !Lattice Units
	DOUBLE PRECISION,PARAMETER :: Tau_y=Y_G*(RHO_S-RHO_F)*(-A_G)*L_C !BN_T*NU*U_C/L_C ! !yield stress BN_T*NU*U_T/L_C

	!***Fluid Reological Properties
	DOUBLE PRECISION,PARAMETER :: Mu_PL=NU !Power-Law Viscosity
	DOUBLE PRECISION,PARAMETER :: N_PL=0.9  !Power-Law power
	DOUBLE PRECISION,PARAMETER :: Mu_0=NU !viscosity
	DOUBLE PRECISION,PARAMETER :: Mu_R=1000.0*Mu_0 !unyeilded viscosity FOR BI-VISCOSITY MODEL
	DOUBLE PRECISION,PARAMETER :: M_PAP=1.0E6 !FOR PAPANASTASIOU MODEL 
	DOUBLE PRECISION,PARAMETER :: Gamma_C=Tau_y/(Mu_R-Mu_0)
	DOUBLE PRECISION           :: OMEGA
	DOUBLE PRECISION           :: OMEGA_MRT(0:8)
	DOUBLE PRECISION           :: M1SM(0:8,0:8)
    DOUBLE PRECISION,PARAMETER :: RHO0=1.0 !Initial Density
	DOUBLE PRECISION,PARAMETER :: Grad_P=0 ! Flow Occurs if Grad_P > Tau_y*Ny/2
    DOUBLE PRECISION,PARAMETER :: dx=1.0
    DOUBLE PRECISION,PARAMETER :: dy=1.0
    DOUBLE PRECISION,PARAMETER :: dt=1.0
	DOUBLE PRECISION,PARAMETER :: G_X=0.0 !3.3333333e-6 !Gravity Acceleration (gx=g)
	DOUBLE PRECISION,PARAMETER :: G_Y=0.0 !Gravity Acceleration (gx=g)
	
	INTEGER :: N_PARTICLE=504
	TYPE PARTICLE
	DOUBLE PRECISION :: CENTER(1:2)
	DOUBLE PRECISION :: ANGLE
	DOUBLE PRECISION :: VELOCITY(1:3)
	DOUBLE PRECISION :: DIAM
	DOUBLE PRECISION :: FORCE(1:3)
	DOUBLE PRECISION :: F_COL(1:2)
	DOUBLE PRECISION :: RHO_P !density
	DOUBLE PRECISION :: M_P   !Mass
	DOUBLE PRECISION :: I_P	  !moment of inertia
	DOUBLE PRECISION :: PHI(0:NX,0:NY)
	DOUBLE PRECISION :: ZITA !Interface Thickness
	END TYPE

END MODULE PARAM 
!=====================================================================
!MAIN PROGRAM FOR COMPUTATIONS OF LATTICE BOLTZMANN METHOD
!=====================================================================
PROGRAM MAIN
    USE PARAM
    IMPLICIT NONE
    DOUBLE PRECISION F(0:Q-1,0:NX,0:NY),FEQ(0:Q-1,0:NX,0:NY),F_POST(0:Q-1,0:NX,0:NY),G(1:2,0:NX,0:NY),STRESS(1:3,0:NX,0:NY)
    DOUBLE PRECISION RHO(0:NX,0:NY),UX(0:NX,0:NY),UY(0:NX,0:NY),U0(0:NX,0:NY),V0(0:NX,0:NY),Omega_BN(0:NX,0:NY)
	LOGICAL :: IS_SOLID_NODE(0:NX,0:NY),IS_WALL_NODE(0:NX,0:NY),IS_PARTICLE_NODE(0:NX,0:NY)
    DOUBLE PRECISION :: AvgDen,SSUM2,ER,START,FINISH,MA,INIT_WALL_CENTER(1:2),PHI_T(0:NX,0:NY),DUMMY,E_R
    DOUBLE PRECISION :: RHO_0_RHO_P,INV_MP,INV_IP
	INTEGER :: KK,I,J,IP,IW,Domain_Image(0:NX,0:NY)
	DOUBLE PRECISION :: F_H_IN(1:3,1:504),VP_OLD(1:3,1:504)
	TYPE(PARTICLE) :: DISK(1:504) !???????
	!DOUBLE PRECISION :: CLUSTER_CENTER(1:2),CLUSTER_U_AVG(1:2),CLUSTER_U_FLC(1:2),CLUSTER_C,CLUSTER_D(1:2)
	!********************************************************************************
    !***COMPUTING RELAXTION FREQUENCY BASED ON DOUBLE PRECISION FLOW REYNOLDS NUMBER
	PRINT*,'In the name of God';
    OMEGA=1.0/(3.*NU+0.5) !RELAXATION FREQUENCY, omega=1/tau
	Omega_BN=Omega
    MA=U_C/SQRT(1.0/3.0)  !Mach Number Ma=Umax/Cs
	OPEN(12,FILE='00_Simulation_Parameters.txt')
	WRITE(12,*) '**********************';WRITE(12,*) 'Lattice Parameter'
	WRITE(12,*) 'TAU_LB=',1.0/OMEGA;WRITE(12,*) 'U_LB_max=',U_C;WRITE(12,*) 'Tau_Y_LB=',TAU_Y;WRITE(12,*) 'AngFreq_LB=',OMG_OSC;
	WRITE(12,*) 'MAX_ITR =',TIME_MAX;WRITE(12,*) '**********************';WRITE(12,*) 'Physical Parameter'
	!WRITE(12,*) 'FREQ (Hz)=',FREQ_OSC/C_T;WRITE(12,*) 'D (cm)=',DIAM_PHYSICAL;WRITE(12,*) 'Rho_S (gr/cm^3)=',RHO_S;
	!WRITE(12,*) 'Ac (cm)=',A_C*C_L;WRITE(12,*) 'NU_PHYSICAL (cm^2/s)=',NU_PHYSICAL;
	WRITE(12,*) '**********************';WRITE(12,*) 'Transformation Coefficients'
	WRITE(12,*) 'C_L (cm)=',C_L;WRITE(12,*) 'C_T (sec)=',C_T;WRITE(12,*) 'C_U (cm/s)=',C_U;WRITE(12,*) 'C_g (m/s^2)=',C_g;
	PRINT*,'Check Simulation Parameters...';
	PRINT *,'If Simulation Parameters are OK continue';PAUSE

    IF(1.0/OMEGA < 0.5) THEN !CHECKING STABILITY CRITERIA FOR RELAXTION FREQUENCY
        PRINT *,'RELAXATION TIME EXCEEDED 0.5 -> STABILITY PROBLEM'
        STOP
    END IF

	!CALL RESULTSPROFILE(UX,UY,RHO);stop
	!********************************************************************************
    OPEN(7,FILE='01_ErrorHistory.PLT')
    OPEN(8,FILE='01_AveragedDensity.PLT')
	OPEN(9,FILE='00_Simulation.TXT')
	WRITE(7,*)'VARIABLES=Iteration,Error'
    WRITE(8,*)'VARIABLES=Iteration,Density_Averaged'
	OPEN(11,FILE='01_ClusterTimehistory.plt')
	WRITE(11,*)'VARIABLES=T,UX_AVG,UY_AVG,UX_FLC,UY_FLC,PHI,DX,DY'
	!OPEN(12,FILE='01_ParticleTimehistory.plt')
	!WRITE(12,*)'VARIABLES=T,Y1,Y2,X1,X2,Ux1,Ux2,Uy1,Uy2,W1,W2,Cd1,Cd2,Dr,DX,DY'
	OPEN(10,FILE='02_Domain_Image.plt')
	WRITE(10,*)'VARIABLES =X, Y, Image'
	OPEN(1001,FILE='01_centerDomain.plt')
	WRITE(1001,*)'VARIABLES =kk,ux,uy'
    CALL CPU_TIME(START)

    !***INITIALIZING VARIABLES**************************
    CALL INITIALIZE(RHO,UX,UY,F,FEQ,G) !Initial fluid variables
    U0=UX
    V0=UY
	IS_SOLID_NODE=.FALSE.
	Domain_Image=0
	!***'Particles"**************************************************************
	OPEN(200,FILE='0_Particle_Initial_Position.txt');
	DO IP=1,N_PARTICLE
	READ(200,*) DISK(IP)%CENTER !,DISK(IP)%DIAM
	!DISK(IP)%CENTER(1)=DISK(IP)%CENTER(1)+WALL%EDGE(4)
	!DISK(IP)%CENTER(2)=DISK(IP)%CENTER(2)+WALL%EDGE(1)
	END DO
	
	DISK%ZITA=1.0
	DISK%DIAM=L_C
	DISK%RHO_P=RHO_S
	DISK%M_P=DISK%RHO_P*PI*(0.5*DISK%DIAM)*(0.5*DISK%DIAM)
	DISK%I_P=0.5*DISK%M_P*(0.5*DISK%DIAM)*(0.5*DISK%DIAM)
	!DO IP=1,N_PARTICLE
	!WRITE(*,*) DISK(IP)%CENTER,DISK(IP)%DIAM,DISK(IP)%RHO_P,DISK(IP)%M_P,DISK(IP)%I_P
	!END DO
	!STOP %%%%%%%%%%%%%%%DKT_VALIDATION%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	!DISK(1)%CENTER=(/0.999/C_L,7.2/C_L/) !initial particle center position
	!DISK(2)%CENTER=(/1.0/C_L,6.8/C_L/)
	!DISK%ANGLE=0D0
	!DISK(1)%VELOCITY=(/0.0,0.0,0.0/) !initial velocity of particle
	!DISK(2)%VELOCITY=(/0.0,0.0,0.0/)
	!DISK%DIAM=L_C
	!DISK(1)%FORCE=(/0.0,0.0,0.0/) !initial forces on particle
	!DISK(2)%FORCE=(/0.0,0.0,0.0/)
	!DISK%ZITA=1.0
	!DISK%RHO_P=RHO_S
	!DISK%M_P=DISK%RHO_P*PI*(0.5*DISK%DIAM)*(0.5*DISK%DIAM)
	!DISK%I_P=0.5*DISK%M_P*(0.5*DISK%DIAM)*(0.5*DISK%DIAM)
	!DO IP=1,N_PARTICLE
	!WRITE(*,*) DISK(IP)%CENTER,DISK(IP)%DIAM,DISK(IP)%RHO_P,DISK(IP)%M_P,DISK(IP)%I_P
	!END DO
	!STOP
	!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	IF(.NOT. READ_INPUT)THEN
	CALL ComputeRelaxation_MRT(OMEGA)
	!call WriteResults(UX,UY,RHO,0) !Writing Initialized Fields
	!PRINT *,'Initial Field is Written'
	!PAUSE
		
	CALL ComputeConcenterationFunction(DISK,DOMAIN_IMAGE,IS_SOLID_NODE,IS_PARTICLE_NODE,IS_WALL_NODE)
	!PHI_t=DISK(1)%PHI+DISK(2)%PHI
	!CALL WRITEIMAGE(WALL%PHI,0)
	!DOMAIN_IMAGE=5;
	!DO J=0,NY
	!DO I=0,NX
	!IF(IS_WALL_NODE(I,J))THEN
	!DOMAIN_IMAGE(I,J)=7
	!ELSEIF(IS_PARTICLE_NODE(I,J))THEN
	!DOMAIN_IMAGE(I,J)=10
	!ELSEIF(IS_SOLID_NODE(I,J))THEN
	!DOMAIN_IMAGE(I,J)=-1
	!END IF
	!END DO
	!END DO
	CALL WRITE_DOMAINIMAGE(DOMAIN_IMAGE,0,10)
	!CALL WRITE_DOMAINIMAGE(IS_WALL_NODE,1,10)
	!CALL WRITE_DOMAINIMAGE(IS_SOLID_NODE,2,10)
	STOP
	IF(0)THEN
	ER=1.0;KK=0;
	!PRINT*,PERIOD_OSC;PAUSE
	DO WHILE(KK<1 )!5.0*PERIOD_OSC)
		KK=KK+1
		PRINT*,KK
	!***PERFORMING COLLISION ON ALL NODES
        !CALL COLLISION(UX,UY,RHO,F,FEQ,G,Omega_BN,IS_SOLID_NODE,DOMAIN_IMAGE)
		CALL Collision_MRT(UX,UY,RHO,F,FEQ,G,Omega_BN,DOMAIN_IMAGE,STRESS)

		!***PERFORMING PROPAGATION ON ALL NODES
        F_POST=F !Post-Collision DFs
        CALL STREAMING(F,F_POST,IS_SOLID_NODE)
        
		!***IMPOSING BOUNDARY CONDITION ON BOUNDARY NODES
        CALL BOUNDARYCONDITION(F,F_POST,RHO,UX,UY,KK)
		
		!***Predicting Velocity
		CALL MACROSCOPICVARIABLES(F,RHO,UX,UY)
		!***Computing Source Feild of SPM if we have Paarticles

		CALL ComputeForceField_SPM(G,DISK,UX,UY)

		!***Computing Hydrodynamics Force on Particles
		CALL ComputeForce_HYDRO(DISK,RHO,UX,UY,G)! Hydro force should be computed with predicted velocity not corrected
		
		!***Correcting Velocity (for Probems with Source-Term)
		CALL VelocityCorrector(UX,UY,G)

		WRITE(1001,*) kk,ux(nx/2,ny/2),uy(nx,ny/2)

		!***Averaged Density
        AvgDen=Sum(Rho)/Dble( (NX+1)*(NY+1) )
		WRITE(8,*) KK,AvgDen

	END DO
	CALL WriteResults(UX,UY,RHO,DOMAIN_IMAGE,STRESS,kk)
	!***WRITING RESULTS
	OPEN(51,FILE='0_F_output.txt');
	OPEN(52,FILE='0_UX_output.txt');
	OPEN(53,FILE='0_UY_output.txt');
	OPEN(54,FILE='0_RHO_output.txt');
	OPEN(55,FILE='0_Gspm_output.txt');
	OPEN(56,FILE='0_OmegaBN_output.txt');
	OPEN(57,FILE='0_Particle_output.txt');
	OPEN(59,FILE='0_OmegaMRT_output.txt');
	OPEN(60,FILE='0_Wall_output.txt');
	OPEN(61,FILE='0_DomainImage_output.txt');

	DO I=0,NX
	DO J=0,NY
	WRITE(51,*) F(:,I,J)
	WRITE(52,*) UX(I,J)
	WRITE(53,*) UY(I,J)
	WRITE(54,*) RHO(I,J)
	WRITE(55,*) G(:,I,J)
	WRITE(56,*) OMEGA_BN(I,J)
	WRITE(61,*) DOMAIN_IMAGE(I,J)
	END DO
	END DO
	DO IP=1,N_PARTICLE
	WRITE(57,*) DISK(IP)%CENTER,DISK(IP)%VELOCITY,DISK(IP)%FORCE
	END DO
	WRITE(59,*) OMEGA_MRT
	END IF
	!***************************************************************************************************
    !***MAIN LOOP***************************************
	 !***READING INPUT
	IF(READ_INPUT)THEN
	OPEN(51,FILE='0_F_output.txt');
	OPEN(52,FILE='0_UX_output.txt');
	OPEN(53,FILE='0_UY_output.txt');
	OPEN(54,FILE='0_RHO_output.txt');
	OPEN(55,FILE='0_Gspm_output.txt');
	OPEN(56,FILE='0_OmegaBN_output.txt');
	OPEN(57,FILE='0_Particle_output.txt');
	OPEN(59,FILE='0_OmegaMRT_output.txt');
	OPEN(61,FILE='0_DomainImage_output.txt');

	DO I=0,NX
	DO J=0,NY
	READ(51,*) F(:,I,J)
	READ(52,*) UX(I,J)
	READ(53,*) UY(I,J)
	READ(54,*) RHO(I,J)
	READ(55,*) G(:,I,J)
	READ(56,*) OMEGA_BN(I,J)
	READ(61,*) DOMAIN_IMAGE(I,J)
	END DO
	END DO
	DO IP=1,N_PARTICLE
	READ(57,*) DISK(IP)%CENTER,DISK(IP)%VELOCITY,DISK(IP)%FORCE
	END DO
	READ(59,*) OMEGA_MRT
	PRINT*,'Reading Initial Data Was Succsessful...';PAUSE
	END IF
	END IF

	KK=0;ER=1.0 !KK=439040;ER=1.32400148569192     
   
    !DO WHILE(ER>1.0E-6) !INT(2.0/C_T)INT(1.0/C_T)INT(10.0/FREQ_OSC)
    DO WHILE(KK<TIME_MAX) !DO WHILE(kk<TIME_MAX)
	PRINT *,KK
		KK=KK+1; !this loop computes variables for time kk+1 (velocities,...)

		!***Computing Particle-Particle and Particle-Wall Collision Force
		CALL ComputeCollisionForce(DISK)

		!***Particle Dynamics Calculation****************
		DO IP=1,N_PARTICLE
		!IP=1
		!***Computing Internal Hydrodynamics force
		F_H_IN(1,IP)=(RHO0/DISK(IP)%RHO_P)*DISK(IP)%M_P*(DISK(IP)%VELOCITY(1)-VP_OLD(1,IP))
		F_H_IN(2,IP)=(RHO0/DISK(IP)%RHO_P)*DISK(IP)%M_P*(DISK(IP)%VELOCITY(2)-VP_OLD(2,IP))
		F_H_IN(3,IP)=(RHO0/DISK(IP)%RHO_P)*DISK(IP)%I_P*(DISK(IP)%VELOCITY(3)-VP_OLD(3,IP))
		!F_H_IN(1,IP)=0d0
		!F_H_IN(2,IP)=0d0
		!F_H_IN(3,IP)=0d0
		VP_OLD(:,IP)=DISK(IP)%VELOCITY
		DISK(IP)%VELOCITY(1)=DISK(IP)%VELOCITY(1)+DISK(IP)%FORCE(1)/DISK(IP)%M_P+F_H_IN(1,IP)/DISK(IP)%M_P+DISK(IP)%F_COL(1)/DISK(IP)%M_P
		DISK(IP)%VELOCITY(2)=DISK(IP)%VELOCITY(2)+DISK(IP)%FORCE(2)/DISK(IP)%M_P+F_H_IN(2,IP)/DISK(IP)%M_P+(1.0-RHO0/DISK(IP)%RHO_P)*A_G+DISK(IP)%F_COL(2)/DISK(IP)%M_P
		DISK(IP)%VELOCITY(3)=DISK(IP)%VELOCITY(3)+DISK(IP)%FORCE(3)/DISK(IP)%I_P+F_H_IN(3,IP)/DISK(IP)%I_P
		DISK(IP)%CENTER(1)=DISK(IP)%CENTER(1)+0.5*( DISK(IP)%VELOCITY(1)+VP_OLD(1,IP) )
		DISK(IP)%CENTER(2)=DISK(IP)%CENTER(2)+0.5*( DISK(IP)%VELOCITY(2)+VP_OLD(2,IP) )
		!DISK(IP)%ANGLE=DISK(IP)%ANGLE+0.5*( DISK(IP)%VELOCITY(3)+VP_OLD(3,IP) )
		END DO
		!WRITE(12,*) KK*C_T,DISK(1)%CENTER(2)*C_L,DISK(2)%CENTER(2)*C_L,(DISK(1)%CENTER(1))*C_L,DISK(2)%CENTER(1)*C_L,&
		!			&DISK(1)%VELOCITY(1)*C_U,DISK(2)%VELOCITY(1)*C_U,DISK(1)%VELOCITY(2)*C_U,DISK(2)%VELOCITY(2)*C_U,&
		!			&DISK(1)%VELOCITY(3),DISK(2)%VELOCITY(3),&
		!			&(DISK(1)%FORCE(1)+F_H_IN(1,1)),(DISK(1)%FORCE(2)+F_H_IN(2,1)),(DISK(1)%FORCE(3)+F_H_IN(3,1)),&
		!			&A_C*OMG_OSC*SIN(OMG_OSC*DBLE(KK))*C_U,C_U
		!***Cluster Calculation**********************************************************
		!CLUSTER_U_AVG=0D0
		!DO IP=1,N_PARTICLE
		!	CLUSTER_U_AVG(1)=DISK(IP)%VELOCITY(1)+CLUSTER_U_AVG(1)
		!	CLUSTER_U_AVG(2)=DISK(IP)%VELOCITY(2)+CLUSTER_U_AVG(2)
		!END DO
		!CLUSTER_U_AVG=CLUSTER_U_AVG/DBLE(N_PARTICLE)
		!CLUSTER_U_FLC=0D0
		!DO IP=1,N_PARTICLE
		!	CLUSTER_U_FLC(1)=(DISK(IP)%VELOCITY(1)-CLUSTER_U_AVG(1))**2+CLUSTER_U_FLC(1)
		!	CLUSTER_U_FLC(2)=(DISK(IP)%VELOCITY(2)-CLUSTER_U_AVG(2))**2+CLUSTER_U_FLC(2)
		!END DO
		!CLUSTER_U_FLC=SQRT(CLUSTER_U_FLC/DBLE(N_PARTICLE))
		
		!DUMMY=(MAXVAL(DISK%CENTER(1))-MINVAL(DISK%CENTER(1))+DISK(1)%DIAM)*&
		!	 &(MAXVAL(DISK%CENTER(2))-MINVAL(DISK%CENTER(2))+DISK(1)%DIAM)
		
		!CLUSTER_C=( DBLE(N_PARTICLE)*PI*(DISK(1)%DIAM**2)/4.0 )/DUMMY

		!CLUSTER_CENTER(1)=SUM(DISK%CENTER(1))/DBLE(N_PARTICLE)
		!CLUSTER_CENTER(2)=SUM(DISK%CENTER(2))/DBLE(N_PARTICLE)
		!CLUSTER_D=0D0
		!DO IP=1,N_PARTICLE
		!	CLUSTER_D(1)=(DISK(IP)%CENTER(1)-CLUSTER_CENTER(1))**2+CLUSTER_D(1)
		!	CLUSTER_D(2)=(DISK(IP)%CENTER(2)-CLUSTER_CENTER(2))**2+CLUSTER_D(2)
		!END DO
		!CLUSTER_D=SQRT(CLUSTER_D/DBLE(N_PARTICLE))

		!WRITE(11,*) KK*C_T,CLUSTER_U_AVG(1)*C_U,CLUSTER_U_AVG(2)*C_U,CLUSTER_U_FLC(1)*C_U,CLUSTER_U_FLC(2)*C_U,&
		!			&CLUSTER_C,CLUSTER_D(1)*C_L,CLUSTER_D(2)*C_L

		!PRINT *,KK,CLUSTER_CENTER(1),CLUSTER_CENTER(2)
					
		!***HydroDynamics Calculation****************
		CALL ComputeConcenterationFunction(DISK,DOMAIN_IMAGE,IS_SOLID_NODE,IS_PARTICLE_NODE,IS_WALL_NODE)
		
		!***PERFORMING COLLISION ON ALL NODES
        !CALL COLLISION(UX,UY,RHO,F,FEQ,G,Omega_BN,IS_SOLID_NODE,DOMAIN_IMAGE)
		CALL Collision_MRT(UX,UY,RHO,F,FEQ,G,Omega_BN,DOMAIN_IMAGE,STRESS)

		!***PERFORMING PROPAGATION ON ALL NODES
        F_POST=F !Post-Collision DFs
        CALL STREAMING(F,F_POST,IS_SOLID_NODE)
        
		!***IMPOSING BOUNDARY CONDITION ON BOUNDARY NODES
        CALL BOUNDARYCONDITION(F,F_POST,RHO,UX,UY,KK)
		
		!***Predicting Velocity
		CALL MACROSCOPICVARIABLES(F,RHO,UX,UY)
		
		!***Computing Source Feild of SPM if we have Particles
		CALL ComputeForceField_SPM(G,DISK,UX,UY)

		!***Computing Hydrodynamics Force on Particles
		CALL ComputeForce_HYDRO(DISK,RHO,UX,UY,G)! Hydro force should be computed with predicted velocity not corrected
		
		!***Correcting Velocity (for Probems with Source-Term)
		CALL VelocityCorrector(UX,UY,G)

		!***Averaged Density
        AvgDen=Sum(Rho)/Dble( (NX+1)*(NY+1) )
		WRITE(8,*) KK,AvgDen

		IF(MOD(KK,1000)==0)THEN	
            CALL COMPUTEERROR(ER,UX,UY,U0,V0)
            PRINT *,'ITERATION=',KK,',','ERROR=',ER
            WRITE(7,*) KK,ER
        END IF

		!IF( KK>=NINT(4.0*PERIOD_OSC) .AND. MOD(KK,TIME_SAMPLE)==0)THEN
		IF(MOD(KK,TIME_SAMPLE)==0)THEN	
			IF(NON_NEWTONIAN)THEN
			DO J=0,NY
			DO I=0,NX
			IF(IS_PARTICLE_NODE(I,J))THEN
			DOMAIN_IMAGE(I,J)=10
			END IF
			END DO
			END DO
			END IF
			CALL WRITE_DOMAINIMAGE(DOMAIN_IMAGE,kk,10)		
			CALL WriteResults(UX,UY,RHO,DOMAIN_IMAGE,STRESS,kk)
		END IF

    END DO
	!************************************************END OF MAIN LOOP*****************
    
	PRINT *,'CONVERGENCE REACHED AT ITERATION=',KK
    !***END OF THE MAIN LOOP

    !***WRITING RESULTS
	IF(WRITE_OUTPUT)THEN
	OPEN(51,FILE='0_F_output.txt');
	OPEN(52,FILE='0_UX_output.txt');
	OPEN(53,FILE='0_UY_output.txt');
	OPEN(54,FILE='0_RHO_output.txt');
	OPEN(55,FILE='0_Gspm_output.txt');
	OPEN(56,FILE='0_OmegaBN_output.txt');
	OPEN(57,FILE='0_Particle_output.txt');
	OPEN(59,FILE='0_OmegaMRT_output.txt');
	OPEN(60,FILE='0_Wall_output.txt');
	OPEN(61,FILE='0_DomainImage_output.txt');

	DO I=0,NX
	DO J=0,NY
	WRITE(51,*) F(:,I,J)
	WRITE(52,*) UX(I,J)
	WRITE(53,*) UY(I,J)
	WRITE(54,*) RHO(I,J)
	WRITE(55,*) G(:,I,J)
	WRITE(56,*) OMEGA_BN(I,J)
	WRITE(61,*) DOMAIN_IMAGE(I,J)
	END DO
	END DO
	DO IP=1,N_PARTICLE
	WRITE(57,*) DISK(IP)%CENTER,DISK(IP)%VELOCITY,DISK(IP)%FORCE
	END DO
	WRITE(59,*) OMEGA_MRT
	!WRITE(60,*) WALL%CENTER,WALL%VELOCITY,WALL%EDGE
	END IF
	IF(NON_NEWTONIAN)THEN
		DO J=0,NY
		DO I=0,NX
		IF(IS_PARTICLE_NODE(I,J))THEN
		DOMAIN_IMAGE(I,J)=10
		END IF
		END DO
		END DO
	END IF
    call WriteResults(UX,UY,RHO,DOMAIN_IMAGE,STRESS,KK);!WRITE(1001,*) KK,WALL%VELOCITY(2)
	!CALL WRITE_DOMAINIMAGE(DOMAIN_IMAGE,KK,10)
	!CALL RESULTSPROFILE(UX,UY,RHO)
    CALL CPU_TIME(FINISH)
    PRINT *,'RUNTIME=',FINISH-START
	WRITE(9,*) 'LastIteration=',KK,'RunTime=',FINISH-START,'LastError=',ER
    stop
end
! end of the main program
!######################################################################################
!####################End Of Main Program#############################################
!######################################################################################
!===================================================================
! Computing solid fraction
!===================================================================
SUBROUTINE ComputeConcenterationFunction(DISK,DOMAIN_IMAGE,IS_SOLID_NODE,IS_PARTICLE_NODE,IS_WALL_NODE)
	USE PARAM,ONLY:NX,NY,DX,PARTICLE,PI,N_PARTICLE
	IMPLICIT NONE
	DOUBLEPRECISION :: X_PC,Y_PC,R_P,D_I,EDGE_TOP,EDGE_BOTTOM,EDGE_RIGHT,EDGE_LEFT,D_X,D_Y,T_P,XPR,YPR,D_P
	INTEGER :: I,J,II,JJ,IP,IW,DOMAIN_IMAGE(0:NX,0:NY)
	TYPE(PARTICLE) :: DISK(1:N_PARTICLE)
	LOGICAL :: IS_SOLID_NODE(0:NX,0:NY),IS_WALL_NODE(0:NX,0:NY),IS_PARTICLE_NODE(0:NX,0:NY)

	DOMAIN_IMAGE=5;IS_SOLID_NODE=.FALSE.;IS_WALL_NODE=.FALSE.;IS_PARTICLE_NODE=.FALSE.
	
	!***Particles
	DO IP=1,N_PARTICLE !LOOP ON PARTICLES
	X_PC=DISK(IP)%CENTER(1) !Particle Center Coordinates
	Y_PC=DISK(IP)%CENTER(2)
	R_P=0.5*DISK(IP)%DIAM !Particle Radius
	T_P=DISK(IP)%ANGLE
	XPR=X_PC+0.75*R_P*COS(T_P)
	YPR=Y_PC+0.75*R_P*SIN(T_P)
	
	DO J=0,NY
		DO I=0,NX
			D_I=R_P-SQRT( (FLOAT(I)-X_PC)*(FLOAT(I)-X_PC)+(FLOAT(J)-Y_PC)*(FLOAT(J)-Y_PC) )
			
			IF( D_I<=-0.5*DISK(IP)%ZITA )THEN
				DISK(IP)%PHI(I,J)=0.0
			ELSEIF( ABS(D_I)<0.5*DISK(IP)%ZITA )THEN
				DISK(IP)%PHI(I,J)=0.5*( SIN( PI*D_I/DISK(IP)%ZITA )+1.0 )
			ELSEIF( D_I>=0.5*DISK(IP)%ZITA )THEN
				DISK(IP)%PHI(I,J)=1.0
			END IF
			IF( DISK(IP)%PHI(I,J)>=0.5 ) THEN
				DOMAIN_IMAGE(I,J)=10
				IS_PARTICLE_NODE(I,J)=.TRUE.
				!D_P=25.D0-( (DBLE(I)-XPR)**2+(DBLE(J)-YPR)**2 )
				!IF(D_P>=0)THEN
				!	DOMAIN_IMAGE(I,J)=20
				!END IF

			END IF
		ENDDO
	ENDDO

	END DO !LOOP ON PARTICLES
END SUBROUTINE ComputeConcenterationFunction
!===================================================================
!Compute Force Field for SPM Method
!===================================================================
SUBROUTINE ComputeForceField_SPM(G,DISK,UX,UY)
	USE PARAM,ONLY:NX,NY,DT,PARTICLE,N_PARTICLE
	IMPLICIT NONE
	DOUBLEPRECISION :: G(1:2,0:NX,0:NY),UX(0:NX,0:NY),UY(0:NX,0:NY)
	DOUBLEPRECISION :: X_PC,Y_PC,UX_PC,UY_PC,OMG_P,UX_P,UY_P,Phi_t(0:NX,0:NY)
	INTEGER :: I,J,IP,IW
	TYPE(PARTICLE) :: DISK(1:N_PARTICLE)
	
	DO J=0,NY
		DO I=0,NX
			G(:,I,J)=0.0
			DO IP=1,N_PARTICLE
				IF(DISK(IP)%PHI(I,J)>0.0)THEN
					UX_P=DISK(IP)%VELOCITY(1)-DISK(IP)%VELOCITY(3)*( FLOAT(J)-DISK(IP)%CENTER(2) )
					UY_P=DISK(IP)%VELOCITY(2)+DISK(IP)%VELOCITY(3)*( FLOAT(I)-DISK(IP)%CENTER(1) )
					G(1,I,J)=DISK(IP)%PHI(I,J)*( UX_P-UX(I,J) )
					G(2,I,J)=DISK(IP)%PHI(I,J)*( UY_P-UY(I,J) )
				END IF
			END DO
			!IF(WALL%PHI(I,J)>0.0)THEN
			!	UX_P=WALL%VELOCITY(1)-WALL%VELOCITY(3)*( FLOAT(J)-WALL%VELOCITY(2) )
			!	UY_P=WALL%VELOCITY(2)+WALL%VELOCITY(3)*( FLOAT(I)-WALL%VELOCITY(1) )
			!	G(1,I,J)=WALL%PHI(I,J)*( UX_P-UX(I,J) )
			!	G(2,I,J)=WALL%PHI(I,J)*( UY_P-UY(I,J) )
			!END IF
		END DO
	END DO
	!DO J=0,NY
	!DO I=0,NX
	!	Phi_t(I,J)=0.0
	!	G(:,I,J)=0.0
	!	DO IP=1,N_PARTILCE
	!	Phi_t(I,J)=Phi_t(I,J)+DISK(IP)%PHI(I,J)
	!	END DO
	!END DO
	!END DO
	
	!DO J=0,NY
	!	DO I=0,NX
			!IF(PHI_T(I,J)>0.0)THEN
				!UX_P=0.0;UY_P=0.0;
				!DO IP=1,2
				!UX_P=UX_P+DISK(IP)%PHI(I,J)*( DISK(IP)%VELOCITY(1)-DISK(IP)%VELOCITY(3)*( FLOAT(J)-DISK(IP)%CENTER(2) ) ) !current point veloctiy
				!UY_P=UY_P+DISK(IP)%PHI(I,J)*( DISK(IP)%VELOCITY(2)+DISK(IP)%VELOCITY(3)*( FLOAT(I)-DISK(IP)%CENTER(1) ) )
				!END DO
				!G(1,I,J)=( UX_P-PHI_T(I,J)*UX(I,J) )
				!G(2,I,J)=( UY_P-PHI_T(I,J)*UY(I,J) )
			
			!ELSE
			!	G(:,I,J)=0.0
			!ENDIF
	!	ENDDO
	!ENDDO
END SUBROUTINE ComputeForceField_SPM
!===================================================================
! Computing Hydrodynamic Force on Particle
!===================================================================
SUBROUTINE ComputeForce_HYDRO(DISK,RHO,UX,UY,G)
	USE PARAM
	IMPLICIT NONE
	DOUBLEPRECISION,INTENT(IN)::G(1:2,0:NX,0:NY),RHO(0:NX,0:NY),UX(0:NX,0:NY),UY(0:NX,0:NY)
	TYPE(PARTICLE) :: DISK(1:N_PARTICLE)
	DOUBLEPRECISION :: F_P(1:3),UX_P,UY_P
	INTEGER :: I,J,IP

	DO IP=1,N_PARTICLE
	F_P=0.0
	DO J=0,NY
	DO I=0,NX
		IF(DISK(IP)%PHI(I,J)>0.0)THEN
			UX_P=DISK(IP)%VELOCITY(1)-DISK(IP)%VELOCITY(3)*( FLOAT(J)-DISK(IP)%CENTER(2) )
			UY_P=DISK(IP)%VELOCITY(2)+DISK(IP)%VELOCITY(3)*( FLOAT(I)-DISK(IP)%CENTER(1) )
			F_P(1)=F_P(1)+RHO(I,J)*DISK(IP)%PHI(I,J)*(UX(I,J)-UX_P)
			F_P(2)=F_P(2)+RHO(I,J)*DISK(IP)%PHI(I,J)*(UY(I,J)-UY_P)
			F_P(3)=F_P(3)+RHO(I,J)*DISK(IP)%PHI(I,J)*( ( FLOAT(I)-DISK(IP)%CENTER(1) )*(UY(I,J)-UY_P)-( FLOAT(J)-DISK(IP)%CENTER(2) )*(UX(I,J)-UX_P) )
		ENDIF
	ENDDO
	ENDDO
	DISK(IP)%FORCE=F_P
	END DO

END SUBROUTINE ComputeForce_HYDRO

!===================================================================
!INITIALIZING MACROSCOPIC QUANTINTIES
!===================================================================
SUBROUTINE INITIALIZE(RHO,UX,UY,F,FEQ,G)
    USE PARAM
    IMPLICIT NONE
    DOUBLE PRECISION,INTENT(INOUT) :: rho(0:nx,0:ny),ux(0:nx,0:ny),uy(0:nx,0:ny),F(0:Q-1,0:NX,0:NY),G(1:2,0:NX,0:NY),FEQ(0:Q-1,0:NX,0:NY)
    INTEGER :: i,j,K
    DOUBLE PRECISION :: CFEQ(0:8)
    do J=0,NY
        do I=0,NX
            rho(i,j)=rho0
            ux(i,j)=0.0
            uy(i,j)=0.0
			CALL COMPUTEFEQNEW(UX(I,J),UY(I,J),RHO(I,J),FEQ(:,I,J))
			F(:,I,J)=FEQ(:,I,J)
        end do
    end do
	G=0.0

END SUBROUTINE INITIALIZE
!===================================================================
!
!===================================================================
SUBROUTINE ComputeCollisionForce(DISK)  
	USE PARAM
	IMPLICIT NONE
	TYPE(PARTICLE) :: DISK(1:N_PARTICLE)
	DOUBLEPRECISION :: ZITA,R_P,R_W(1:2),R_PJ,D_PW,CIJ,D_PP,F_PP(1:2),F_PW(1:2),EPS_P,EPS_W,E_P,E_W
	INTEGER :: IP,JP,IW
	
	DO IP=1,N_PARTICLE
	DISK(IP)%F_COL=0.0
	END DO
	!DISK(2)%F_COL=0.0
	!*******************L-J MODEL**************
	!IF(0)THEN
	!ZITA=1.0;R_P=0.5*DISK(1)%DIAM;

	!IP=1;JP=2;
	!		D_PP=SQRT( ( DISK(IP)%CENTER(1)-DISK(JP)%CENTER(1) )**2+( DISK(IP)%CENTER(2)-DISK(JP)%CENTER(2) )**2 )
	!		IF( D_PP>(0.5*DISK(IP)%DIAM+0.5*DISK(JP)%DIAM+ZITA) )THEN
	!			F_PP=0.0
	!		ELSEIF( D_PP<=(0.5*DISK(IP)%DIAM+0.5*DISK(JP)%DIAM+ZITA) )THEN
	!			F_PP(1)=2.4*(R_P**2)*(2.0*(2.0*R_P/D_PP)**14-(2.0*R_P/D_PP)**8)*(DISK(IP)%CENTER(1)-DISK(JP)%CENTER(1))/(4.0*R_P*R_P)
	!			F_PP(2)=2.4*(R_P**2)*(2.0*(2.0*R_P/D_PP)**14-(2.0*R_P/D_PP)**8)*(DISK(IP)%CENTER(2)-DISK(JP)%CENTER(2))/(4.0*R_P*R_P)
	!		ENDIF
	!		DISK(IP)%F_COL=+F_PP
	!		DISK(JP)%F_COL=-F_PP 
	!DO IP=1,N_PARTICLE
	!	DO IW=1,4 !Index of each Wall
	!		IF(IW==1)THEN
	!		R_W(1)=DISK(IP)%CENTER(1);R_W(2)=0.0 ! DISK(IP)%CENTER(2) !Glowinsky
	!		ELSEIF(IW==2)THEN
	!		R_W(1)=DBLE(NX);R_W(2)=DISK(IP)%CENTER(2)  !2.0*DBLE(NX)-DISK(IP)%CENTER(1)
	!		ELSEIF(IW==3)THEN
	!		R_W(1)=DISK(IP)%CENTER(1);R_W(2)=2.0*DBLE(NY)-DISK(IP)%CENTER(2)
	!		ELSEIF(IW==4)THEN
	!		R_W(1)=-DISK(IP)%CENTER(1);R_W(2)=DISK(IP)%CENTER(2)
	!		END IF
	!		D_PW=SQRT( (DISK(IP)%CENTER(1)-R_W(1))**2+(DISK(IP)%CENTER(2)-R_W(2))**2 ) !Particle-Wall Distance
			
	!		IF( D_PW>(DISK(IP)%DIAM+ZITA) )THEN
	!			F_PW=0.0
	!		ELSEIF( D_PW<=(DISK(IP)%DIAM+ZITA) )THEN
	!			F_PW(1)=2.4*(R_P**2)*(2.0*(R_P/D_PW)**14-(R_P/D_PW)**8)*(DISK(IP)%CENTER(1)-R_W(1))/(R_P*R_P)
	!			F_PW(2)=2.4*(R_P**2)*(2.0*(R_P/D_PW)**14-(R_P/D_PW)**8)*(DISK(IP)%CENTER(2)-R_W(2))/(R_P*R_P)
	!		ENDIF
	!		DISK(IP)%F_COL=DISK(IP)%F_COL+F_PW
	!	END DO
	!END DO
	!END IF
	!*******************W-T MODEL**************
	!IF(0)THEN !WAN AND TUREK MODEL
	!ZITA=1.0;CIJ=1.0;EPS_P=1.0;EPS_W=2.0*EPS_P;
	!DO IP=1,N_PARTICLE
		!DO JPP=IP+1,N_PARTICLE
		!	IF(IP==N_PARTICLE)THEN
		!		JP=JPP-N_PARTICLE+1
		!	END IF
	!IF(N_PARTICLE>1)THEN
	!IP=1;JP=2;
	!		D_PP=SQRT( ( DISK(IP)%CENTER(1)-DISK(JP)%CENTER(1) )**2+( DISK(IP)%CENTER(2)-DISK(JP)%CENTER(2) )**2 )	
	!		IF( D_PP>(0.5*DISK(IP)%DIAM+0.5*DISK(JP)%DIAM+ZITA) )THEN
	!			F_PP=0.0
	!		ELSEIF( D_PP<(0.5*DISK(IP)%DIAM+0.5*DISK(JP)%DIAM) )THEN
	!			F_PP(1)=EPS_P*(DISK(IP)%CENTER(1)-DISK(JP)%CENTER(1))*(0.5*DISK(IP)%DIAM+0.5*DISK(JP)%DIAM-D_PP)
	!			F_PP(2)=EPS_P*(DISK(IP)%CENTER(2)-DISK(JP)%CENTER(2))*(0.5*DISK(IP)%DIAM+0.5*DISK(JP)%DIAM-D_PP)
	!		ELSE
	!			F_PP(1)=EPS_P*(DISK(IP)%CENTER(1)-DISK(JP)%CENTER(1))*(0.5*DISK(IP)%DIAM+0.5*DISK(JP)%DIAM+ZITA-D_PP)**2
	!			F_PP(2)=EPS_P*(DISK(IP)%CENTER(2)-DISK(JP)%CENTER(2))*(0.5*DISK(IP)%DIAM+0.5*DISK(JP)%DIAM+ZITA-D_PP)**2
	!		ENDIF
	!		DISK(IP)%F_COL=+F_PP
	!		DISK(JP)%F_COL=-F_PP
	!	END DO
	!END DO
	!END IF
	!DO IP=1,N_PARTICLE
	!	DO IW=1,4 !Index of each Wall
	!		IF(IW==1)THEN
			!R_W(1)=DISK(IP)%CENTER(1);R_W(2)=WALL%HEIGHT-0.5*DISK(IP)%DIAM !Glowinsky
			!R_W(1)=DISK(IP)%CENTER(1);R_W(2)=2.0*WALL%EDGE(IW)-DISK(IP)%CENTER(2)
	!		R_W(1)=DISK(IP)%CENTER(1);R_W(2)=-DISK(IP)%CENTER(2)
	!		ELSEIF(IW==2)THEN
			!R_W(1)=DBLE(NX)+0.5*DISK(IP)%DIAM;R_W(2)=DISK(IP)%CENTER(2)
			!R_W(1)=2.0*WALL%EDGE(IW)-DISK(IP)%CENTER(1);R_W(2)=DISK(IP)%CENTER(2)
	!		R_W(1)=DBLE(NX)+DISK(IP)%CENTER(1);R_W(2)=DISK(IP)%CENTER(2)
	!		ELSEIF(IW==3)THEN
			!R_W(1)=DISK(IP)%CENTER(1);R_W(2)=DBLE(NY)+0.5*DISK(IP)%DIAM
			!R_W(1)=DISK(IP)%CENTER(1);R_W(2)=2.0*WALL%EDGE(IW)-DISK(IP)%CENTER(2)
	!		ELSEIF(IW==4)THEN
			!R_W(1)=-0.5*DISK(IP)%DIAM;R_W(2)=DISK(IP)%CENTER(2)
			!R_W(1)=2.0*WALL%EDGE(IW)-DISK(IP)%CENTER(1);R_W(2)=DISK(IP)%CENTER(2)
	!		END IF
	!		D_PW=SQRT( (DISK(IP)%CENTER(1)-R_W(1))**2+(DISK(IP)%CENTER(2)-R_W(2))**2 ) !Particle-Wall Distance
			
	!		IF( D_PW>(DISK(IP)%DIAM+ZITA) )THEN
	!			F_PW=0.0
	!		ELSEIF( D_PW<(DISK(IP)%DIAM) )THEN
	!			F_PW(1)=EPS_W*(DISK(IP)%CENTER(1)-R_W(1))*(DISK(IP)%DIAM-D_PW)
	!			F_PW(2)=EPS_W*(DISK(IP)%CENTER(2)-R_W(2))*(DISK(IP)%DIAM-D_PW)
	!		ELSE
	!			F_PW(1)=EPS_W*(DISK(IP)%CENTER(1)-R_W(1))*(DISK(IP)%DIAM+ZITA-D_PW)**2
	!			F_PW(2)=EPS_W*(DISK(IP)%CENTER(2)-R_W(2))*(DISK(IP)%DIAM+ZITA-D_PW)**2
	!		ENDIF
	!		DISK(IP)%F_COL=DISK(IP)%F_COL+F_PW
	!	END DO
	!END DO
	!END IF
	!*******************GLOWINSKY MODEL**************
	!IF(1)THEN
	ZITA=1.0;CIJ=-(RHO_S-RHO_F)*A_G*PI*0.5*DISK(1)%DIAM*0.5*DISK(1)%DIAM !A_G*PI*0.5*DISK(1)%DIAM*0.5*DISK(1)%DIAM;
	!EPS_P=1.0E2;EPS_W=0.5*EPS_P;E_P=2.0E2;E_W=0.5*E_P
	!EPS_P=0.5;EPS_W=EPS_P;E_P=10.0;E_W=E_P
	EPS_P=2.0;EPS_W=0.5*EPS_P;E_P=20.0;E_W=0.5*E_P
	
	!IF(N_PARTICLE>1)THEN
	!IP=1;JP=2;
		DO IP=1,N_PARTICLE
			IF(IP==N_PARTICLE)THEN
				EXIT
			END IF
			DO JP=IP+1,N_PARTICLE
			D_PP=SQRT( ( DISK(IP)%CENTER(1)-DISK(JP)%CENTER(1) )**2+( DISK(IP)%CENTER(2)-DISK(JP)%CENTER(2) )**2 )	
			!IF( D_PP>(0.5*DISK(IP)%DIAM+0.5*DISK(JP)%DIAM+ZITA) )THEN
			!	F_PP=0.0
			!ELSEIF( D_PP<=(0.5*DISK(IP)%DIAM+0.5*DISK(JP)%DIAM+ZITA) )THEN
			!	F_PP(1)=0.5*CIJ*((D_PP-0.5*DISK(IP)%DIAM-0.5*DISK(JP)%DIAM-ZITA)/ZITA)**2 *(DISK(IP)%CENTER(1)-DISK(JP)%CENTER(1))/D_PP
			!	F_PP(2)=0.5*CIJ*((D_PP-0.5*DISK(IP)%DIAM-0.5*DISK(JP)%DIAM-ZITA)/ZITA)**2 *(DISK(IP)%CENTER(2)-DISK(JP)%CENTER(2))/D_PP
			!ENDIF
			IF( D_PP>(0.5*DISK(IP)%DIAM+0.5*DISK(JP)%DIAM+ZITA) )THEN
				F_PP=0.0
			ELSEIF( D_PP<=(0.5*DISK(IP)%DIAM+0.5*DISK(JP)%DIAM+ZITA) .AND. D_PW>(0.5*DISK(IP)%DIAM+0.5*DISK(JP)%DIAM) )THEN
				F_PP(1)=EPS_P*CIJ*((D_PP-0.5*DISK(IP)%DIAM-0.5*DISK(JP)%DIAM-ZITA)/ZITA)**2 *(DISK(IP)%CENTER(1)-DISK(JP)%CENTER(1))/D_PP
				F_PP(2)=EPS_P*CIJ*((D_PP-0.5*DISK(IP)%DIAM-0.5*DISK(JP)%DIAM-ZITA)/ZITA)**2 *(DISK(IP)%CENTER(2)-DISK(JP)%CENTER(2))/D_PP
			ELSEIF( D_PW<=(0.5*DISK(IP)%DIAM+0.5*DISK(JP)%DIAM) )THEN
				F_PP(1)=( EPS_P*CIJ*((D_PP-0.5*DISK(IP)%DIAM-0.5*DISK(JP)%DIAM-ZITA)/ZITA)**2 +E_P*CIJ*(0.5*DISK(IP)%DIAM+0.5*DISK(JP)%DIAM-D_PP)/ZITA )*( DISK(IP)%CENTER(1)-DISK(JP)%CENTER(1) )/D_PP
				F_PP(2)=( EPS_P*CIJ*((D_PP-0.5*DISK(IP)%DIAM-0.5*DISK(JP)%DIAM-ZITA)/ZITA)**2 +E_P*CIJ*(0.5*DISK(IP)%DIAM+0.5*DISK(JP)%DIAM-D_PP)/ZITA )*( DISK(IP)%CENTER(2)-DISK(JP)%CENTER(2) )/D_PP
			ENDIF
			DISK(IP)%F_COL=DISK(IP)%F_COL+F_PP
			DISK(JP)%F_COL=DISK(JP)%F_COL-F_PP
			END DO
		END DO
	!END IF

	DO IP=1,N_PARTICLE
		DO IW=1,4 !Index of each Wall
			IF(IW==1)THEN !BOTTOM WALL
			R_W(1)=DISK(IP)%CENTER(1);R_W(2)=-0.5*DISK(IP)%DIAM !DISK(IP)%CENTER(2)
			ELSEIF(IW==2)THEN !RIGHT WALL
			R_W(1)=DBLE(NX)+0.5*DISK(IP)%DIAM;R_W(2)=DISK(IP)%CENTER(2)  !2.0*DBLE(NX)-DISK(IP)%CENTER(1)
			ELSEIF(IW==3)THEN !TOP WALL
			R_W(1)=DISK(IP)%CENTER(1);R_W(2)=DBLE(NY)+0.5*DISK(IP)%DIAM
			ELSEIF(IW==4)THEN !LEFT WALL
			R_W(1)=-0.5*DISK(IP)%DIAM;R_W(2)=DISK(IP)%CENTER(2)
			END IF
			D_PW=SQRT( (DISK(IP)%CENTER(1)-R_W(1))**2+(DISK(IP)%CENTER(2)-R_W(2))**2 ) !Particle-Wall Distance

			!IF( D_PW>(DISK(IP)%DIAM+ZITA) )THEN
			!	F_PW=0.0
			!ELSEIF( D_PW<=(DISK(IP)%DIAM+ZITA) )THEN
			!	F_PW(1)=0.25*CIJ*(( (D_PW-DISK(IP)%DIAM-ZITA)/ZITA )**2)*( DISK(IP)%CENTER(1)-R_W(1) )/D_PW
			!	F_PW(2)=0.25*CIJ*(( (D_PW-DISK(IP)%DIAM-ZITA)/ZITA )**2)*( DISK(IP)%CENTER(2)-R_W(2) )/D_PW
			!ENDIF

			IF( D_PW>(DISK(IP)%DIAM+ZITA) )THEN
				F_PW=0.0
			ELSEIF( D_PW<=(DISK(IP)%DIAM+ZITA) .AND. D_PW>(DISK(IP)%DIAM) )THEN
				F_PW(1)=EPS_W*CIJ*(( (D_PW-DISK(IP)%DIAM-ZITA)/ZITA )**2)*( DISK(IP)%CENTER(1)-R_W(1) )/D_PW
				F_PW(2)=EPS_W*CIJ*(( (D_PW-DISK(IP)%DIAM-ZITA)/ZITA )**2)*( DISK(IP)%CENTER(2)-R_W(2) )/D_PW
			ELSEIF( D_PW<=(DISK(IP)%DIAM) )THEN
				F_PW(1)=( EPS_W*CIJ*( (D_PW-DISK(IP)%DIAM-ZITA)/ZITA )**2+E_W*CIJ*(DISK(IP)%DIAM-D_PW)/ZITA )*( DISK(IP)%CENTER(1)-R_W(1) )/D_PW
				F_PW(2)=( EPS_W*CIJ*( (D_PW-DISK(IP)%DIAM-ZITA)/ZITA )**2+E_W*CIJ*(DISK(IP)%DIAM-D_PW)/ZITA )*( DISK(IP)%CENTER(2)-R_W(2) )/D_PW
			ENDIF
			DISK(IP)%F_COL=DISK(IP)%F_COL+F_PW
		END DO
	END DO
	!END IF
END SUBROUTINE ComputeCollisionForce
!===================================================================
!COMPUTING POST-COLLISION DISTRIBUTION FUNCTION AT EACH NODE IN THE PRESENT TIME STEP
!(THIS SUBROUTINE NEVER CHANGES)
!===================================================================
SUBROUTINE COLLISION(UX,UY,RHO,F,FEQ,G,Omega_BN,IS_SOLID_NODE,DOMAIN_IMAGE)
    use param,only : nx,ny,q,OMEGA,CX,CY,Mu_R,Mu_0,Tau_Y,Gamma_C,M_PAP,NON_NEWTONIAN
    IMPLICIT NONE
    DOUBLE PRECISION,intent(in) :: rho(0:nx,0:ny),ux(0:nx,0:ny),uy(0:nx,0:ny),G(1:2,0:NX,0:NY)
	LOGICAL,INTENT(INOUT) :: IS_SOLID_NODE(0:NX,0:NY)
	INTEGER,intent(inout) :: DOMAIN_IMAGE(0:NX,0:NY)
    DOUBLE PRECISION,intent(inout) :: f(0:q-1,0:nx,0:ny),feq(0:q-1,0:nx,0:ny),OMEGA_BN(0:NX,0:NY)

    DOUBLE PRECISION :: FI(0:8),ZARIB,F_1(0:8),Nu_A,D11,D22,D12,Gamma,TAU_11,TAU_22,TAU_12,TAU_MAG
    INTEGER :: i,j

    DO J=0,nY
        DO I=0,nX

			IF(NON_NEWTONIAN)THEN
			ZARIB=-1.5*OMEGA_BN(I,J)/RHO(I,J)
			F_1(:)=F(:,I,J)-FEQ(:,I,J) !Non-equilibrium f
			D11=ZARIB*(F_1(1)+F_1(5)+F_1(6)+F_1(3)+F_1(7)+F_1(8)) !rate of deformation tensor Components d_ij
			D22=ZARIB*(F_1(5)+F_1(2)+F_1(6)+F_1(7)+F_1(4)+F_1(8))
			D12=ZARIB*(F_1(5)-F_1(6)+F_1(7)-F_1(8))
			!GAMMA=SQRT(0.5*ABS(D11**2+D22**2+2.0*D12**2));
			!GAMMA=SQRT(2.0*ABS(D11**2+D22**2+2.0*D12**2)); !Shear Rate (Original Definition)
			GAMMA=SQRT(4.0*ABS(D12**2-D11*D22)); !Shear Rate (Only for Incompressible Flow)
			!***Biviscosity Model**************************************
			IF(GAMMA<GAMMA_C)THEN !un-yeilded Region
				NU_A=MU_R !High reference Viscosity
				DOMAIN_IMAGE(I,J)=5 !SOLID
			ELSEIF(GAMMA>=GAMMA_C)THEN !yeilded Region
				NU_A=MU_0+TAU_Y/GAMMA
				DOMAIN_IMAGE(I,J)=0 !FLUID
			END IF
			!***Papanastasiou Model***********************************
			!IF(GAMMA<1.0E-9)THEN !zero shear rate limit 
			!	Nu_A=Mu_0+M_PAP*TAU_Y !un-yielded Viscosity
				!IS_SOLID_NODE(I,J)=.TRUE.
			!	DOMAIN_IMAGE(I,J)=1 !SOLID
			!ELSE
			!	Nu_A=Mu_0+(TAU_Y/GAMMA)*(1.0-EXP(-M_PAP*GAMMA)) 
				!IS_SOLID_NODE(I,J)=.FALSE.
			!	DOMAIN_IMAGE(I,J)=0 !FLUID
			!END IF
			!TAU_11=NU_A*2.0*D11;TAU_22=NU_A*2.0*D22;TAU_12=NU_A*2.0*D12;
			!TAU_11=TAU_11-0.5*(1.0-0.5*OMEGA_BN(I,J))*(2.0*G(1,I,J)*UX(I,J))
			!TAU_22=TAU_22-0.5*(1.0-0.5*OMEGA_BN(I,J))*(2.0*G(2,I,J)*UY(I,J))
			!TAU_12=TAU_12-0.5*(1.0-0.5*OMEGA_BN(I,J))*(G(1,I,J)*UY(I,J)+G(2,I,J)*UX(I,J))
			!TAU_MAG=SQRT(0.5*ABS(TAU_11**2+TAU_22**2+2.0*TAU_12**2))
			!IF(TAU_MAG<=TAU_Y)THEN !zero shear rate limit 
				!Nu_A=Mu_0+M_PAP*TAU_Y !un-yielded Viscosity
				!IS_SOLID_NODE(I,J)=.TRUE.
			!	DOMAIN_IMAGE(I,J)=5 !SOLID
			!ELSE
				!Nu_A=Mu_0+(TAU_Y/GAMMA)*(1.0-EXP(-M_PAP*GAMMA)) 
				!IS_SOLID_NODE(I,J)=.FALSE.
			!	DOMAIN_IMAGE(I,J)=0 !FLUID
			!END IF
			!*********************************************************
			OMEGA_BN(I,J)=1.0/(3.0*NU_A+0.5) !important:: no density should appear heare except rho0 (if needed)
			CALL ComputeSourceEquilib_Hydro(UX(I,J),UY(I,J),RHO(I,J),G(1,I,J),G(2,I,J),FI,FEQ(:,I,J))				
			!F(:,i,j)=omega*feq(:,i,j)+(1.0-omega)*f(:,i,j)+(1.-0.5*OMEGA)*FI(:)
			F(:,i,j)=omega_BN(I,J)*feq(:,i,j)+(1.0-omega_BN(I,J))*f(:,i,j)+(1.-0.5*OMEGA_BN(I,J))*FI(:)
			ELSE
				CALL ComputeSourceEquilib_Hydro(UX(I,J),UY(I,J),RHO(I,J),G(1,I,J),G(2,I,J),FI,FEQ(:,I,J))				
				F(:,i,j)=omega*feq(:,i,j)+(1.0-omega)*f(:,i,j)+(1.-0.5*OMEGA)*FI(:)
			END IF
        END DO
    END DO

END SUBROUTINE COLLISION
!************************************************************************************************************
SUBROUTINE ComputeSourceEquilib_Hydro(UX,UY,RHO,GX,GY,F,Feq)
	USE PARAM,ONLY:G_X,G_Y
	IMPLICIT NONE
	DOUBLEPRECISION::UX,UY,RHO,GX,GY,F(0:8),FEQ(0:8),W0,W1,W2,U2,INV_CLB=1.0
	
	W0=4.0/9.0*inv_clb*inv_clb;W1=1./9.*inv_clb*inv_clb;W2=1./36.*inv_clb*inv_clb
	F(0)=W0*RHO*(-3.*UX*(GX+G_X)-3*UY*(GY+G_Y))
	F(1)=W1*RHO*((3.+6.*UX)*(GX+G_X)-3.*UY*(GY+G_Y))
	F(2)=W1*RHO*(-3.*UX*(GX+G_X)+(3.+6.*UY)*(GY+G_Y))
	F(3)=W1*RHO*((-3.+6.*UX)*(GX+G_X)-3.*UY*(GY+G_Y))
	F(4)=W1*RHO*(-3.*UX*(GX+G_X)+(-3.+6.*UY)*(GY+G_Y))
	F(5)=W2*RHO*((3.+6.*UX+9.*UY)*(GX+G_X)+(3.+9.*UX+6.*UY)*(GY+G_Y))
	F(6)=W2*RHO*((-3.+6.*UX-9.*UY)*(GX+G_X)+(3.-9.*UX+6.*UY)*(GY+G_Y))
	F(7)=W2*RHO*((-3.+6.*UX+9.*UY)*(GX+G_X)+(-3.+9.*UX+6.*UY)*(GY+G_Y))
	F(8)=W2*RHO*((3.+6.*UX-9.*UY)*(GX+G_X)+(-3.-9.*UX+6.*UY)*(GY+G_Y))
	
	U2=UX*UX+UY*UY
	W0=2./9.;W1=1./18.;W2=1./36.
	FEQ(0)=W0*RHO*(2.              -3.*U2*inv_clb*inv_clb)
	FEQ(1)=W1*RHO*(2.+6.*UX*inv_clb+9.*UX*UX*inv_clb*inv_clb-3.*U2*inv_clb*inv_clb)
	FEQ(2)=W1*RHO*(2.+6.*UY*inv_clb+9.*UY*UY*inv_clb*inv_clb-3.*U2*inv_clb*inv_clb)
	FEQ(3)=W1*RHO*(2.-6.*UX*inv_clb+9.*UX*UX*inv_clb*inv_clb-3.*U2*inv_clb*inv_clb)
	FEQ(4)=W1*RHO*(2.-6.*UY*inv_clb+9.*UY*UY*inv_clb*inv_clb-3.*U2*inv_clb*inv_clb)
	FEQ(5)=W2*RHO*(1.+3.*(UX+UY)*inv_clb+9.*UX*UY*inv_clb*inv_clb+3.*U2*inv_clb*inv_clb)
	FEQ(6)=W2*RHO*(1.-3.*(UX-UY)*inv_clb-9.*UX*UY*inv_clb*inv_clb+3.*U2*inv_clb*inv_clb)
	FEQ(7)=W2*RHO*(1.-3.*(UX+UY)*inv_clb+9.*UX*UY*inv_clb*inv_clb+3.*U2*inv_clb*inv_clb)
	FEQ(8)=W2*RHO*(1.+3.*(UX-UY)*inv_clb-9.*UX*UY*inv_clb*inv_clb+3.*U2*inv_clb*inv_clb)

END SUBROUTINE ComputeSourceEquilib_Hydro
!==========================================================================================
!MRT-Collision Subroutines
!==========================================================================================
SUBROUTINE Collision_MRT(UX,UY,RHO,F,FEQ,G,Omega_BN,Domain_Image,STRESS)
    use param,only : nx,ny,Omega_mrt,MU_0,M_PAP,TAU_Y,CX,CY,mu_r,gamma_c,NON_NEWTONIAN,NU
    IMPLICIT NONE
    DOUBLE PRECISION,intent(in) :: rho(0:nx,0:ny),ux(0:nx,0:ny),uy(0:nx,0:ny)
    DOUBLE PRECISION,intent(inout) :: f(0:8,0:nx,0:ny),feq(0:8,0:nx,0:ny),G(1:2,0:NX,0:NY),Omega_BN(0:NX,0:NY),Stress(1:3,0:NX,0:NY)
	INTEGER,INTENT(INOUT) :: Domain_Image(0:NX,0:NY)
    DOUBLE PRECISION :: M(0:8),MF(0:8),MEQ(0:8),M_NE(0:8),F_1(0:8),ZARIB,Inv_Rho,D11,D22,D12,GAMMA,NU_A,SUM,TAU_11,TAU_22,TAU_12,TAU_MAG,d_mrt=1./36.
    INTEGER :: i,j,K,II,JJ

    DO J=0,nY
        DO I=0,nX
			M(0)=F(0,I,J)+F(1,I,J)+F(2,I,J)+F(3,I,J)+F(4,I,J)+F(5,I,J)+F(6,I,J)+F(7,I,J)+F(8,I,J)
			M(1)=-4.0*F(0,I,J)-F(1,I,J)-F(2,I,J)-F(3,I,J)-F(4,I,J)+2.0*(F(5,I,J)+F(6,I,J)+F(7,I,J)+F(8,I,J))
			M(2)=4.0*F(0,I,J)-2.0*(F(1,I,J)+F(2,I,J)+F(3,I,J)+F(4,I,J))+(F(5,I,J)+F(6,I,J)+F(7,I,J)+F(8,I,J))
			M(3)=F(1,I,J)-F(3,I,J)+F(5,I,J)-F(6,I,J)-F(7,I,J)+F(8,I,J)
			M(4)=-2.0*(F(1,I,J)-F(3,I,J))+F(5,I,J)-F(6,I,J)-F(7,I,J)+F(8,I,J)
			M(5)=F(2,I,J)-F(4,I,J)+F(5,I,J)+F(6,I,J)-F(7,I,J)-F(8,I,J)
			M(6)=-2.0*(F(2,I,J)-F(4,I,J))+F(5,I,J)+F(6,I,J)-F(7,I,J)-F(8,I,J)
			M(7)=F(1,I,J)-F(2,I,J)+F(3,I,J)-F(4,I,J)
			M(8)=F(5,I,J)-F(6,I,J)+F(7,I,J)-F(8,I,J)

			CALL ComputeSourceEquilib_MRT(UX(I,J),UY(I,J),RHO(I,J),G(1,I,J),G(2,I,J),MF,MEQ)
			M_NE(:)=(M(:)-MEQ(:));
			m(:)=omega_mrt(:)*MEQ(:)+(1.-omega_mrt(:))*m(:)+(1.-0.5*Omega_mrt(:))*MF(:)
			
			F(0,I,J)=d_mrt*( 4.0*(m(0)-m(1)+m(2)) )
			F(1,I,J)=d_mrt*( 4.*m(0)-m(1)-2.*m(2)+6.*m(3)-6.*m(4)+9.*m(7) )
			F(2,I,J)=d_mrt*( 4.0*M(0)-M(1)-2.0*M(2)+6.0*M(5)-6.0*M(6)-9.0*M(7) )
			F(3,I,J)=d_mrt*( 4.0*M(0)-M(1)-2.0*M(2)-6.0*M(3)+6.0*M(4)+9.0*M(7) )
			F(4,I,J)=d_mrt*( 4.0*M(0)-M(1)-2.0*M(2)-6.0*M(5)+6.0*M(6)-9.0*M(7) )
			F(5,I,J)=d_mrt*( 4.0*M(0)+2.0*M(1)+M(2)+6.0*M(3)+3.0*M(4)+6.0*M(5)+3.0*M(6)+9.0*M(8) )
			F(6,I,J)=d_mrt*( 4.0*M(0)+2.0*M(1)+M(2)-6.0*M(3)-3.0*M(4)+6.0*M(5)+3.0*M(6)-9.0*M(8) )
			F(7,I,J)=d_mrt*( 4.0*M(0)+2.0*M(1)+M(2)-6.0*M(3)-3.0*M(4)-6.0*M(5)-3.0*M(6)+9.0*M(8) )
			F(8,I,J)=d_mrt*( 4.0*M(0)+2.0*M(1)+M(2)+6.0*M(3)+3.0*M(4)-6.0*M(5)-3.0*M(6)-9.0*M(8) )

			!***Chai Method for Strain Rate Tensor
			F_1(0)=4.0*OMEGA_MRT(0)*M_NE(0)-4.0*OMEGA_MRT(1)*M_NE(1)+4.0*OMEGA_MRT(2)*M_NE(2)
			F_1(1)=4.0*OMEGA_MRT(0)*M_NE(0)-1.0*OMEGA_MRT(1)*M_NE(1)-2.0*OMEGA_MRT(2)*M_NE(2)+6.0*OMEGA_MRT(3)*M_NE(3)-6.0*OMEGA_MRT(4)*M_NE(4)+9.0*OMEGA_MRT(7)*M_NE(7)
			F_1(2)=4.0*OMEGA_MRT(0)*M_NE(0)-1.0*OMEGA_MRT(1)*M_NE(1)-2.0*OMEGA_MRT(2)*M_NE(2)+6.0*OMEGA_MRT(5)*M_NE(5)-6.0*OMEGA_MRT(6)*M_NE(6)-9.0*OMEGA_MRT(7)*M_NE(7)
			F_1(3)=4.0*OMEGA_MRT(0)*M_NE(0)-1.0*OMEGA_MRT(1)*M_NE(1)-2.0*OMEGA_MRT(2)*M_NE(2)-6.0*OMEGA_MRT(3)*M_NE(3)+6.0*OMEGA_MRT(4)*M_NE(4)+9.0*OMEGA_MRT(7)*M_NE(7)
			F_1(4)=4.0*OMEGA_MRT(0)*M_NE(0)-1.0*OMEGA_MRT(1)*M_NE(1)-2.0*OMEGA_MRT(2)*M_NE(2)-6.0*OMEGA_MRT(5)*M_NE(5)+6.0*OMEGA_MRT(6)*M_NE(6)-9.0*OMEGA_MRT(7)*M_NE(7)
			F_1(5)=4.0*OMEGA_MRT(0)*M_NE(0)+2.0*OMEGA_MRT(1)*M_NE(1)+1.0*OMEGA_MRT(2)*M_NE(2)+6.0*OMEGA_MRT(3)*M_NE(3)+3.0*OMEGA_MRT(4)*M_NE(4)+6.0*OMEGA_MRT(5)*M_NE(5)+3.0*OMEGA_MRT(6)*M_NE(6)+9.0*OMEGA_MRT(8)*M_NE(8)
			F_1(6)=4.0*OMEGA_MRT(0)*M_NE(0)+2.0*OMEGA_MRT(1)*M_NE(1)+1.0*OMEGA_MRT(2)*M_NE(2)-6.0*OMEGA_MRT(3)*M_NE(3)-3.0*OMEGA_MRT(4)*M_NE(4)+6.0*OMEGA_MRT(5)*M_NE(5)+3.0*OMEGA_MRT(6)*M_NE(6)-9.0*OMEGA_MRT(8)*M_NE(8)
			F_1(7)=4.0*OMEGA_MRT(0)*M_NE(0)+2.0*OMEGA_MRT(1)*M_NE(1)+1.0*OMEGA_MRT(2)*M_NE(2)-6.0*OMEGA_MRT(3)*M_NE(3)-3.0*OMEGA_MRT(4)*M_NE(4)-6.0*OMEGA_MRT(5)*M_NE(5)-3.0*OMEGA_MRT(6)*M_NE(6)+9.0*OMEGA_MRT(8)*M_NE(8)
			F_1(8)=4.0*OMEGA_MRT(0)*M_NE(0)+2.0*OMEGA_MRT(1)*M_NE(1)+1.0*OMEGA_MRT(2)*M_NE(2)+6.0*OMEGA_MRT(3)*M_NE(3)+3.0*OMEGA_MRT(4)*M_NE(4)-6.0*OMEGA_MRT(5)*M_NE(5)-3.0*OMEGA_MRT(6)*M_NE(6)-9.0*OMEGA_MRT(8)*M_NE(8)
			
			ZARIB=-1.5*d_mrt
			D11=ZARIB*(F_1(1)+F_1(5)+F_1(6)+F_1(3)+F_1(7)+F_1(8))+0.75*( (OMEGA_MRT(7)-OMEGA_MRT(1))*(G(1,I,J)*UX(I,J)+G(2,I,J)*UY(I,J))-OMEGA_MRT(7)*(2.0*G(1,I,J)*UX(I,J)) )
			D22=ZARIB*(F_1(5)+F_1(2)+F_1(6)+F_1(7)+F_1(4)+F_1(8))+0.75*( (OMEGA_MRT(7)-OMEGA_MRT(1))*(G(1,I,J)*UX(I,J)+G(2,I,J)*UY(I,J))-OMEGA_MRT(7)*(2.0*G(2,I,J)*UY(I,J)) )
			D12=ZARIB*(F_1(5)-F_1(6)+F_1(7)-F_1(8))+0.75*( -OMEGA_MRT(7)*(G(1,I,J)*UY(I,J)+G(2,I,J)*UX(I,J)) )

			TAU_11=NU*2.0*D11;TAU_22=NU*2.0*D22;TAU_12=NU*2.0*D12; !Stress Tensor Component
			STRESS(1,I,J)=TAU_11;STRESS(2,I,J)=TAU_22;STRESS(3,I,J)=TAU_12;

			IF(NON_NEWTONIAN)THEN
			INV_RHO=1.0/RHO(I,J);D11=D11*INV_RHO;D22=D22*INV_RHO;D12=D12*INV_RHO

			GAMMA=SQRT(2.0*ABS(D11**2+D22**2+2.0*D12**2)); !Shear Rate (Original Definition)
			!GAMMA=SQRT(4.0*ABS(D12**2-D11*D22)); !Shear Rate (Only for Incompressible Flow)
			!***Power-Law Model**************************************
			!NU_A=(MU_PL*(GAMMA**(N_PL-1.0)))/RHO(I,J);
			!NU_A=(MU_PL*(GAMMA**(N_PL-1.0d0)));
			!***Biviscosity Model**************************************
			!IF(GAMMA<=GAMMA_C)THEN !un-yeilded Region
			!	NU_A=MU_R !High reference Viscosity
			!	DOMAIN_IMAGE(I,J)=5 !SOLID
			!ELSEIF(GAMMA>GAMMA_C)THEN !yeilded Region
			!	NU_A=MU_0+TAU_Y/GAMMA
			!	DOMAIN_IMAGE(I,J)=0 !FLUID
			!END IF
			!NU_A=NU_A/RHO(I,J)
			!***Papanastasiou Model***********************************
			IF(GAMMA<1.0E-15)THEN !zero shear rate limit !HERE NU_A IS DYNAMIC VISCOSITY (MU) 
				!Nu_A=Mu_0+M_PAP*TAU_Y !Bingham
				Nu_A=SQRT(Mu_0)+SQRT(M_PAP*TAU_Y) !Casson Kinematic Viscosity
			ELSE
				!Nu_A=Mu_0+(TAU_Y/GAMMA)*(1.0-EXP(-M_PAP*GAMMA)) !Bingham
				Nu_A=( SQRT(Mu_0)+SQRT(TAU_Y/GAMMA)*(1.0-EXP(-SQRT(M_PAP*GAMMA))) )**2 !Casson Kinematic Viscosity
			END IF
			
			TAU_11=NU_A*2.0*D11;TAU_22=NU_A*2.0*D22;TAU_12=NU_A*2.0*D12; !Stress Tensor Component
			STRESS(1,I,J)=TAU_11;STRESS(2,I,J)=TAU_22;STRESS(3,I,J)=TAU_12;
			TAU_MAG=SQRT(0.5*ABS(TAU_11**2+TAU_22**2+2.0*TAU_12**2)) !Mag of a Tensor=sqrt(0.5*trace(tensor^2))
			IF(TAU_MAG<=TAU_Y)THEN !Finding Yeilded Region
				DOMAIN_IMAGE(I,J)=5 !SOLID
			ELSE
				DOMAIN_IMAGE(I,J)=0 !FLUID
			END IF
			!*********************************************************
			NU_A=NU_A*Inv_Rho !Kinematic Viscosity=Dynamic Viscosity/Density
			OMEGA_BN(I,J)=1.0/(3.0*NU_A+0.5) !important:: no density should appear heare except rho0 (if needed)
			CALL ComputeRelaxation_MRT(OMEGA_BN(I,J))
			END IF
        END DO
    END DO
END SUBROUTINE Collision_MRT

SUBROUTINE ComputeSourceEquilib_MRT(UX,UY,RHO,GX,GY,F,Feq)
	USE PARAM,ONLY:G_X,G_Y
	IMPLICIT NONE
	DOUBLEPRECISION::UX,UY,RHO,GX,GY,F(0:8),FEQ(0:8),F_X,F_Y,U2
	F_X=Rho*(GX+G_X)
	F_Y=Rho*(GY+G_Y)
	F(0)=0.0
	F(1)=6.0*(F_X*UX+F_Y*UY)
	F(2)=-6.0*(F_X*UX+F_Y*UY)
	F(3)=F_X
	F(4)=-F_X
	F(5)=F_Y
	F(6)=-F_Y
	F(7)=2.0*(F_X*UX-F_Y*UY)
	F(8)=F_Y*UX+F_X*UY
	
	U2=UX*UX+UY*UY
	FEQ(0)=RHO
	FEQ(1)=RHO*(-2.0+3.*U2)
	FEQ(2)=RHO*(1.0-3.*U2)
	FEQ(3)=RHO*( UX)
	FEQ(4)=RHO*(-UX)
	FEQ(5)=RHO*( UY)
	FEQ(6)=RHO*(-UY)
	FEQ(7)=RHO*(UX*UX-UY*UY)
	FEQ(8)=RHO*(UX*UY)
	!FEQ(0)=RHO
	!FEQ(1)=RHO*(1.0-3.0*U2)
	!FEQ(2)=RHO*(9.0*(UX*UX*UY*UY)-3.0*U2+1.0)
	!FEQ(3)=RHO*( UX)
	!FEQ(4)=RHO*(3.0*UX*UX*UX-UX)
	!FEQ(5)=RHO*( UY)
	!FEQ(6)=RHO*(3.0*UY*UY*UY-UY)
	!FEQ(7)=RHO*(UX*UX-UY*UY)
	!FEQ(8)=RHO*(UX*UY)

END SUBROUTINE ComputeSourceEquilib_MRT

SUBROUTINE ComputeRelaxation_MRT(OMEGA)
    USE PARAM,ONLY : OMEGA_MRT
    IMPLICIT NONE
    INTEGER :: i,j,K
    DOUBLE PRECISION :: omega_e,omega_eps,omega_q,OMEGA
	
	!omega_e=1.6;omega_eps=1.8;omega_q=8.*(2.-omega)/(8.-omega)
	omega_e=1.1;omega_eps=1.25;omega_q=8.*(2.-omega)/(8.-omega)
	Omega_mrt(0)=0.0;Omega_mrt(1)=omega_e;Omega_mrt(2)=omega_eps;Omega_mrt(3)=0.0
	Omega_mrt(4)=omega_q;Omega_mrt(5)=0.0;Omega_mrt(6)=omega_q;Omega_mrt(7)=omega;Omega_mrt(8)=omega;
	!****Chai (2011) MRT Viscoplastic
	!omega_e=1.1;omega_eps=1.0;omega_q=1.2
	!Omega_mrt(0)=0.0;Omega_mrt(1)=omega_e;Omega_mrt(2)=omega_eps;Omega_mrt(3)=0.0
	!Omega_mrt(4)=omega_q;Omega_mrt(5)=0.0;Omega_mrt(6)=omega_q;Omega_mrt(7)=omega;Omega_mrt(8)=omega;
	!****Mohammad MRT
	!Omega_mrt(0)=1.0;Omega_mrt(1)=1.4;Omega_mrt(2)=1.4;Omega_mrt(3)=1.0
	!Omega_mrt(4)=1.2;Omega_mrt(5)=1.0;Omega_mrt(6)=1.2;Omega_mrt(7)=omega;Omega_mrt(8)=omega;


END SUBROUTINE ComputeRelaxation_MRT

SUBROUTINE ComputeM1SM_MRT(Omega,M1SM)
    USE PARAM,ONLY : OMEGA_MRT
    IMPLICIT NONE
    DOUBLE PRECISION,INTENT(OUT) :: M1SM(0:8,0:8)
	DOUBLE PRECISION,INTENT(IN) ::Omega
	DOUBLE PRECISION :: S1,S2,S3,S4,S5,S6,S7,S8,S9,OMEGA_E,OMEGA_EPS,OMEGA_Q

	!omega_e=1.6;omega_eps=1.8;omega_q=8.*(2.-omega)/(8.-omega)
	!Omega_mrt(0)=0.0;Omega_mrt(1)=omega_e;Omega_mrt(2)=omega_eps;Omega_mrt(3)=0.0
	!Omega_mrt(4)=omega_q;Omega_mrt(5)=0.0;Omega_mrt(6)=omega_q;Omega_mrt(7)=omega;Omega_mrt(8)=omega;
	!****Chai (2011) MRT Viscoplastic
	omega_e=1.1;omega_eps=1.0;omega_q=1.2
	Omega_mrt(0)=0.0;Omega_mrt(1)=omega_e;Omega_mrt(2)=omega_eps;Omega_mrt(3)=0.0
	Omega_mrt(4)=omega_q;Omega_mrt(5)=0.0;Omega_mrt(6)=omega_q;Omega_mrt(7)=omega;Omega_mrt(8)=omega;
	S1=OMEGA_MRT(0);S2=OMEGA_MRT(1);S3=OMEGA_MRT(2);S4=OMEGA_MRT(3);S5=OMEGA_MRT(4);S6=OMEGA_MRT(5);S7=OMEGA_MRT(6);S8=OMEGA_MRT(7);S9=OMEGA_MRT(8);
	M1SM(0,:)=(/4.0*s1+16.0*s2+16.0*s3,4.0*s1+4.0*s2-8.0*s3,4.0*s1+4.0*s2-8.0*s3,4.0*s1+4.0*s2-8.0*s3,4.0*s1+4.0*s2-8.0*s3,4.0*s1-8.0*s2+4.0*s3,4.0*s1-8.0*s2+4.0*s3,4.0*s1-8.0*s2+4.0*s3,4.0*s1-8.0*s2+4.0*s3/);
	M1SM(1,:)=(/4.0*s1+4.0*s2-8.0*s3,4.0*s1+s2+4.0*s3+6.0*s4+12.0*s5+9.0*s8,4.0*s1+s2+4.0*s3-9.0*s8,4.0*s1+s2+4.0*s3-6.0*s4-12.0*s5+9.0*s8,4.0*s1+s2+4.0*s3-9.0*s8,4.0*s1-2.0*s2-2.0*s3+6.0*s4-6.0*s5,4.0*s1-2.0*s2-2.0*s3-6.0*s4+6.0*s5,4.0*s1-2.0*s2-2.0*s3-6.0*s4+6.0*s5,4.0*s1-2.0*s2-2.0*s3+6.0*s4-6.0*s5/);
	M1SM(2,:)=(/4.0*s1+4.0*s2-8.0*s3,4.0*s1+s2+4.0*s3-9.0*s8,4.0*s1+s2+4.0*s3+6.0*s6+12.0*s7+9.0*s8,4.0*s1+s2+4.0*s3-9.0*s8,4.0*s1+s2+4.0*s3-6.0*s6-12.0*s7+9.0*s8,4.0*s1-2.0*s2-2.0*s3+6.0*s6-6.0*s7,4.0*s1-2.0*s2-2.0*s3+6.0*s6-6.0*s7,4.0*s1-2.0*s2-2.0*s3-6.0*s6+6.0*s7,4.0*s1-2.0*s2-2.0*s3-6.0*s6+6.0*s7/);
	M1SM(3,:)=(/4.0*s1+4.0*s2-8.0*s3,4.0*s1+s2+4.0*s3-6.0*s4-12.0*s5+9.0*s8,4.0*s1+s2+4.0*s3-9.0*s8,4.0*s1+s2+4.0*s3+6.0*s4+12.0*s5+9.0*s8,4.0*s1+s2+4.0*s3-9.0*s8,4.0*s1-2.0*s2-2.0*s3-6.0*s4+6.0*s5,4.0*s1-2.0*s2-2.0*s3+6.0*s4-6.0*s5,4.0*s1-2.0*s2-2.0*s3+6.0*s4-6.0*s5,4.0*s1-2.0*s2-2.0*s3-6.0*s4+6.0*s5/);
	M1SM(4,:)=(/4.0*s1+4.0*s2-8.0*s3,4.0*s1+s2+4.0*s3-9.0*s8,4.0*s1+s2+4.0*s3-6.0*s6-12.0*s7+9.0*s8,4.0*s1+s2+4.0*s3-9.0*s8,4.0*s1+s2+4.0*s3+6.0*s6+12.0*s7+9.0*s8,4.0*s1-2.0*s2-2.0*s3-6.0*s6+6.0*s7,4.0*s1-2.0*s2-2.0*s3-6.0*s6+6.0*s7,4.0*s1-2.0*s2-2.0*s3+6.0*s6-6.0*s7,4.0*s1-2.0*s2-2.0*s3+6.0*s6-6.0*s7/);
	M1SM(5,:)=(/4.0*s1-8.0*s2+4.0*s3,4.0*s1-2.0*s2-2.0*s3+6.0*s4-6.0*s5,4.0*s1-2.0*s2-2.0*s3+6.0*s6-6.0*s7,4.0*s1-2.0*s2-2.0*s3-6.0*s4+6.0*s5,4.0*s1-2.0*s2-2.0*s3-6.0*s6+6.0*s7,4.0*s1+4.0*s2+s3+6.0*s4+3.0*s5+6.0*s6+3.0*s7+9.0*s9,4.0*s1+4.0*s2+s3-6.0*s4-3.0*s5+6.0*s6+3.0*s7-9.0*s9,4.0*s1+4.0*s2+s3-6.0*s4-3.0*s5-6.0*s6-3.0*s7+9.0*s9,4.0*s1+4.0*s2+s3+6.0*s4+3.0*s5-6.0*s6-3.0*s7-9.0*s9/);
	M1SM(6,:)=(/4.0*s1-8.0*s2+4.0*s3,4.0*s1-2.0*s2-2.0*s3-6.0*s4+6.0*s5,4.0*s1-2.0*s2-2.0*s3+6.0*s6-6.0*s7,4.0*s1-2.0*s2-2.0*s3+6.0*s4-6.0*s5,4.0*s1-2.0*s2-2.0*s3-6.0*s6+6.0*s7,4.0*s1+4.0*s2+s3-6.0*s4-3.0*s5+6.0*s6+3.0*s7-9.0*s9,4.0*s1+4.0*s2+s3+6.0*s4+3.0*s5+6.0*s6+3.0*s7+9.0*s9,4.0*s1+4.0*s2+s3+6.0*s4+3.0*s5-6.0*s6-3.0*s7-9.0*s9,4.0*s1+4.0*s2+s3-6.0*s4-3.0*s5-6.0*s6-3.0*s7+9.0*s9/);
	M1SM(7,:)=(/4.0*s1-8.0*s2+4.0*s3,4.0*s1-2.0*s2-2.0*s3-6.0*s4+6.0*s5,4.0*s1-2.0*s2-2.0*s3-6.0*s6+6.0*s7,4.0*s1-2.0*s2-2.0*s3+6.0*s4-6.0*s5,4.0*s1-2.0*s2-2.0*s3+6.0*s6-6.0*s7,4.0*s1+4.0*s2+s3-6.0*s4-3.0*s5-6.0*s6-3.0*s7+9.0*s9,4.0*s1+4.0*s2+s3+6.0*s4+3.0*s5-6.0*s6-3.0*s7-9.0*s9,4.0*s1+4.0*s2+s3+6.0*s4+3.0*s5+6.0*s6+3.0*s7+9.0*s9,4.0*s1+4.0*s2+s3-6.0*s4-3.0*s5+6.0*s6+3.0*s7-9.0*s9/);
	M1SM(8,:)=(/4.0*s1-8.0*s2+4.0*s3,4.0*s1-2.0*s2-2.0*s3+6.0*s4-6.0*s5,4.0*s1-2.0*s2-2.0*s3-6.0*s6+6.0*s7,4.0*s1-2.0*s2-2.0*s3-6.0*s4+6.0*s5,4.0*s1-2.0*s2-2.0*s3+6.0*s6-6.0*s7,4.0*s1+4.0*s2+s3+6.0*s4+3.0*s5-6.0*s6-3.0*s7-9.0*s9,4.0*s1+4.0*s2+s3-6.0*s4-3.0*s5-6.0*s6-3.0*s7+9.0*s9,4.0*s1+4.0*s2+s3-6.0*s4-3.0*s5+6.0*s6+3.0*s7-9.0*s9,4.0*s1+4.0*s2+s3+6.0*s4+3.0*s5+6.0*s6+3.0*s7+9.0*s9/);

	M1SM=(1.0/36.0)*M1SM;

END SUBROUTINE ComputeM1SM_MRT

!================================================================
! Computation of Source Term
!================================================================
SUBROUTINE ImplementSource(G,F,UX,UY,RHO)
	USE PARAM
	IMPLICIT NONE
	DOUBLE PRECISION,INTENT(IN) :: G(1:2,0:NX,0:NY),UX(0:NX,0:NY),UY(0:NX,0:NY),RHO(0:NX,0:NY)
	DOUBLE PRECISION,INTENT(INOUT) :: F(0:Q-1,0:NX,0:NY)
	DOUBLE PRECISION :: S,FK,FI(0:8)
	INTEGER :: I,J,K
	
	DO J=0,NY
		DO I=0,NX
			CALL ComputeSourceTermNew(UX(I,J),UY(I,J),RHO(I,J),G(1,I,J),G(2,I,J),FI)
			DO K=0,Q-1
				!***Gou Forcing Scheme (Validated for Gravity Poise)
				S=(1.-0.5*OMEGA)*FI(K)
				!***LGA Forcing Scheme
				!S=3.0*W(K)*RHO(I,J)*DT*(G(1,I,J)*CX(K)+G(2,I,J)*CY(K))
				!**********************
				F(K,I,J)=F(K,I,J)+S
			END DO
		END DO
	END DO
END SUBROUTINE
SUBROUTINE ComputeSourceTermNew(UX,UY,RHO,GX,GY,F)
	USE PARAM,ONLY:G_X,G_Y
	IMPLICIT NONE
	DOUBLEPRECISION::UX,UY,RHO,GX,GY,F(0:8),W0,W1,W2,A

	W0=4.0/9.0
	W1=1./9.
	W2=1./36.
	IF(1)THEN
	F(0)=W0*RHO*(-3.*UX*(GX+G_X)-3*UY*(GY+G_Y))
	F(1)=W1*RHO*((3.+6.*UX)*(GX+G_X)-3.*UY*(GY+G_Y))
	F(2)=W1*RHO*(-3.*UX*(GX+G_X)+(3.+6.*UY)*(GY+G_Y))
	F(3)=W1*RHO*((-3.+6.*UX)*(GX+G_X)-3.*UY*(GY+G_Y))
	F(4)=W1*RHO*(-3.*UX*(GX+G_X)+(-3.+6.*UY)*(GY+G_Y))
	F(5)=W2*RHO*((3.+6.*UX+9.*UY)*(GX+G_X)+(3.+9.*UX+6.*UY)*(GY+G_Y))
	F(6)=W2*RHO*((-3.+6.*UX-9.*UY)*(GX+G_X)+(3.-9.*UX+6.*UY)*(GY+G_Y))
	F(7)=W2*RHO*((-3.+6.*UX+9.*UY)*(GX+G_X)+(-3.+9.*UX+6.*UY)*(GY+G_Y))
	F(8)=W2*RHO*((3.+6.*UX-9.*UY)*(GX+G_X)+(-3.-9.*UX+6.*UY)*(GY+G_Y))
	ENDIF
	A=1./6.
	IF(0)THEN
	F(0)=A*RHO*(0.0)
	F(1)=A*RHO*(GX+G_X)
	F(2)=A*RHO*(GY+G_Y)
	F(3)=A*RHO*-1.*(GX+G_X)
	F(4)=A*RHO*-1.*(GY+G_Y)
	F(5)=A*RHO*((GX+G_X)+(GY+G_Y))
	F(6)=A*RHO*(-1.*(GX+G_X)+(GY+G_Y))
	F(7)=A*RHO*(-1.*(GX+G_X)-1.*(GY+G_Y))
	F(8)=A*RHO*((GX+G_X)-1.*(GY+G_Y))
	ENDIF
END SUBROUTINE
!===============================================================================

SUBROUTINE VelocityCorrector(UX,UY,G)
	USE PARAM
	IMPLICIT NONE
	DOUBLE PRECISION :: UX(0:NX,0:NY),UY(0:NX,0:NY),G(1:2,0:NX,0:NY)

	INTEGER :: I,J

	DO J=0,NY
		DO I=0,NX
			!***High Order Forcing Scheme
			UX(I,J)=UX(I,J)+0.5*DT*(G(1,I,J)+G_X)
			UY(I,J)=UY(I,J)+0.5*DT*(G(2,I,J)+G_Y)
		END DO
	END DO
END SUBROUTINE

!================================================================
! COMPUTE EQULIBRIUM DISTRIBUTION FUNCTION
!================================================================
SUBROUTINE COMPUTEFEQNEW(UX,UY,RHO,FEQ)
	IMPLICIT NONE
	DOUBLEPRECISION :: UX,UY,RHO,FEQ(0:8),W0,W1,W2,U2

	U2=UX*UX+UY*UY
	W0=2./9.
	W1=1./18.
	W2=1./36.
	FEQ(0)=W0*RHO*(2.              -3.*U2)
	FEQ(1)=W1*RHO*(2.+6.*UX+9.*UX*UX-3.*U2)
	FEQ(2)=W1*RHO*(2.+6.*UY+9.*UY*UY-3.*U2)
	FEQ(3)=W1*RHO*(2.-6.*UX+9.*UX*UX-3.*U2)
	FEQ(4)=W1*RHO*(2.-6.*UY+9.*UY*UY-3.*U2)
	FEQ(5)=W2*RHO*(1.+3.*(UX+UY)+9.*UX*UY+3.*U2)
	FEQ(6)=W2*RHO*(1.-3.*(UX-UY)-9.*UX*UY+3.*U2)
	FEQ(7)=W2*RHO*(1.-3.*(UX+UY)+9.*UX*UY+3.*U2)
	FEQ(8)=W2*RHO*(1.+3.*(UX-UY)-9.*UX*UY+3.*U2)

END SUBROUTINE

!================================================================
!COMPUTING DISTRUBUTION FUNCTION OF EACH NODE FOR THE NEXT TIME STEP
!NO   HAPPENS FOR UNKNOWN DIRECTIONS OF (THIS SUBROUTINE NEVER CHANGES)
!================================================================
SUBROUTINE STREAMING(F,F_POST,IS_SOLID_NODE)
    USE PARAM,ONLY : NX,NY,Q,CX,CY
    IMPLICIT NONE

    DOUBLE PRECISION,INTENT(IN) :: F_POST(0:Q-1,0:NX,0:NY)
	LOGICAL,INTENT(IN) :: IS_SOLID_NODE(0:NX,0:NY)
    DOUBLE PRECISION,INTENT(INOUT) :: F(0:Q-1,0:NX,0:NY)

    INTEGER :: I,J,K,IO,JO,OPP(0:8)=(/0,3,4,1,2,7,8,5,6/),F_TEMP(0:8)

    DO J=0,NY
        DO I=0,NX
			F_TEMP=F(:,I,J)
            DO K=0,Q-1
                IO=I-CX(K)
                JO=J-CY(K)
                IF(IO>=0 .AND. IO<=NX .AND. JO>=0 .AND. JO<=NY ) THEN !ONLY ON FLUID NODES WE HAVE STREAMING
                    F(K,I,J)=F_POST(K,IO,JO)
					!IF(.NOT.(IS_SOLID_NODE(IO,JO)))THEN
					!F(K,I,J)=F_POST(K,IO,JO)
					!ELSE
					!F(K,I,J)=F_POST(OPP(K),I,J)					
					!END IF
					
				!IF(.NOT.(IS_Solid_Node(I,J)) .AND. .NOT.(IS_Solid_Node(IO,JO)) )THEN
				!	F(K,I,J)=F_POST(K,IO,JO)
				!ELSE IF( .NOT.(IS_Solid_Node(I,J)) .AND. IS_Solid_Node(IO,JO))THEN
				!	F(K,I,J)=F_TEMP(OPP(K))
				!END IF

                END IF

            END DO
        END DO
    END DO

END SUBROUTINE STREAMING

!===================================================================
!COMPUTING UNKNOWN DISTRIBUTION FUNCTION ON BOUNDARY NODES
!===================================================================
SUBROUTINE BOUNDARYCONDITION(F,F_POST,RHO,UX,UY,KK)
    USE PARAM
    IMPLICIT NONE
    DOUBLE PRECISION,INTENT(INOUT) :: F(0:Q-1,0:NX,0:NY),F_POST(0:Q-1,0:NX,0:NY),RHO(0:NX,0:NY),UX(0:NX,0:NY),UY(0:NX,0:NY)

    INTEGER :: I,J,KK
    DOUBLE PRECISION :: RHON,RHOS,UXN,UXS,UYN,UYS,UX_L,UY_L,RHO_L,UX_R,UY_R,RHO_R,UX_W,UY_W,RHO_W,c_lb,Inv_Six,YY,Y0
	c_lb=1.0;Inv_Six=1.0/6.0;
    
	UXN=0.0;UYN=A_C*OMG_OSC*SIN(OMG_OSC*DBLE(KK))
	UXS=0.0;UYS=A_C*OMG_OSC*SIN(OMG_OSC*DBLE(KK))
	DO I=1,NX-1
		!CALL COMPUTEFEQNEW(UX(I,1),UY(I,1),RHO(I,1),F(:,I,0))
		!CALL COMPUTEFEQNEW(UX(I,NY-1),UY(I,NY-1),RHO(I,NY-1),F(:,I,NY))
		!F(:,I,NY)=FEQ(:)
		!F(:,I,0 )=FEQ(:)
		!***Halfway Bounca-Back
		!***TOP WALL
		F(8,I,NY)=F_POST(6,I,NY)-Inv_six*(c_lb**2)*Rho(I,Ny)*( -UxN+UyN )
        F(4,I,NY)=F_POST(2,I,NY)-4.0*Inv_six*(c_lb**2)*Rho(I,Ny)*(  UyN )
        F(7,I,NY)=F_POST(5,I,NY)-Inv_six*(c_lb**2)*Rho(I,Ny)*(  UxN+UyN )
		!***BOTTOM WALL
		F(5,I,0)=F_POST(7,I,0)-Inv_six*(c_lb**2)*Rho(I,0)*( -UxS-UyS )
        F(2,I,0)=F_POST(4,I,0)-4.0*Inv_six*(c_lb**2)*Rho(I,0)*( -UyS )
        F(6,I,0)=F_POST(8,I,0)-Inv_six*(c_lb**2)*Rho(I,0)*( UxS-UyS )
		!*******************************************************************************
		!***Zou-He Boundary Condition
		!***TOP MOVING WALL
        !Rhon=(f(0,i,ny)+f(1,i,ny)+f(3,i,ny)+2.*(f(2,i,ny)+f(6,i,ny)+f(5,i,ny)))/(1.0+UYN)
        !f(4,i,ny)=f(2,i,ny)-4.*Inv_six*RHON*UYN
        !f(8,i,ny)=f(6,i,ny)-0.5*(f(1,i,ny)-f(3,i,ny))+0.5*RHON*UXN-Inv_six*RHON*UYN
        !f(7,i,ny)=f(5,i,ny)+0.5*(f(1,i,ny)-f(3,i,ny))-0.5*RHON*UXN-Inv_six*RHON*UYN
		!***BOTTOM MOVING WALL
		!RHOS=(f(0,i,0)+f(1,i,0)+f(3,i,0)+2.*(f(4,i,0)+f(7,i,0)+f(8,i,0)))/(1.0-UYS)
        !f(2,i,0)=f(4,i,0)+4.0*Inv_six*RHOS*UYS
        !f(5,i,0)=f(7,i,0)-0.5*(f(1,i,0)-f(3,i,0))+0.5*RHOS*UXS+Inv_six*RHOS*UYS
        !f(6,i,0)=f(8,i,0)+0.5*(f(1,i,0)-f(3,i,0))-0.5*RHOS*UXS+Inv_six*RHOS*UYS
		!**********************************************************************
		!***Top Periodic
		!F(8,I,NY)=F(8,I,0)
		!F(4,I,NY)=F(4,I,0)
        !F(7,I,NY)=F(7,I,0)
        !***Bottom PERIODIC 
        !F(5,I,0)=F(8,I,0)
		!F(2,I,0)=F(4,I,0)
        !F(6,I,0)=F(7,I,0)
		!**********************************************************************
		!***Top Periodic
		!F(8,I,NY)=F_POST(8,I-1,0)
		!F(4,I,NY)=F_POST(4,I  ,0)
        !F(7,I,NY)=F_POST(7,I+1,0)
        !***Bottom PERIODIC 
        !F(5,I,0)=F_POST(5,I-1,NY)
		!F(2,I,0)=F_POST(2,I  ,NY)
        !F(6,I,0)=F_POST(6,I+1,NY)
		!************************
		!F(:,I,NY)=F(:,I,NY-1)
		!F(:,I,0)=F(:,I,1)
		!***Top Periodic
		!F(8,I,NY)=F(8,I,NY-1)
		!F(4,I,NY)=F(4,I,NY-1)
        !F(7,I,NY)=F(7,I,NY-1)
        !***Bottom PERIODIC 
        !F(5,I,0)=F(5,I,1)
		!F(2,I,0)=F(2,I,1)
        !F(6,I,0)=F(6,I,1)
    END DO
    
	UX_L=0.0;UY_L=A_C*OMG_OSC*SIN(OMG_OSC*DBLE(KK));
	Ux_R=0.0;Uy_R=A_C*OMG_OSC*SIN(OMG_OSC*DBLE(KK));
    DO J=1,NY-1
		!CALL COMPUTEFEQNEW(UX(1,J),UY(1,J),RHO(1,J),F(:,0,J))
		!CALL COMPUTEFEQNEW(UX(NX-1,J),UY(NX-1,J),RHO(NX-1,J),F(:,NX,J))
		!***Left
		!F(5,0,J)=F(6,1,J)
		!F(1,0,J)=F(3,1,J)
        !F(8,0,J)=F(7,1,J)
		!***Right
        !F(6,NX,J)=F(8,NX-1,J)
        !F(3,NX,J)=F(1,NX-1,J)
        !F(7,NX,J)=F(5,NX-1,J)
		!F(:,NX,J)=F(:,NX,J-1)
		!F(:,0 ,J)=F(:,0 ,J-1)
		!F(:,NX,J)=FEQ(:)
	    !F(:,0 ,J)=FEQ(:)
	!***HalfWay Bounce-Back
		!***Left
		!F(5,0,J)=F_POST(7,0,J)
		!F(1,0,J)=F_POST(3,0,J)
        !F(8,0,J)=F_POST(6,0,J)
		F(6,NX,J)=F_POST(8,NX,J)-Inv_six*(c_lb**2)*Rho(NX,J)*( Ux_R-Uy_R )
		F(3,NX,J)=F_POST(1,NX,J)-4.0*Inv_six*(c_lb**2)*Rho(NX,J)*(  UX_R )
		F(7,NX,J)=F_POST(5,NX,J)-Inv_six*(c_lb**2)*Rho(NX,J)*( Ux_R+Uy_R )
        !***Right
        !F(6,NX,J)=F_POST(8,NX,J)
        !F(3,NX,J)=F_POST(1,NX,J)
        !F(7,NX,J)=F_POST(5,NX,J)
		F(5,0,J)=F_POST(7,0,J)-Inv_six*(c_lb**2)*Rho(0,J)*( -Ux_L-Uy_L )
        F(1,0,J)=F_POST(3,0,J)-4.0*Inv_six*(c_lb**2)*Rho(0,J)*(  -UX_L )
        F(8,0,J)=F_POST(6,0,J)-Inv_six*(c_lb**2)*Rho(0,J)*( -Ux_L+Uy_L )
	!***RIGHT BOUNDARY = OUTLET (ZERO-GRADIENT)
		!F(1,NX,J)=2.0*F(1,NX-1,J)-F(1,NX-2,J)
        !F(5,NX,J)=2.0*F(5,NX-1,J)-F(5,NX-2,J)
        !F(8,NX,J)=2.0*F(8,NX-1,J)-F(8,NX-2,J)
	!****************************************************	
	!***LEFT BOUNDARY = VELOCITY INLET (EQUILIBRIUM DISTRUBUTION FUNCTION)
		!***Poise Inlet
		!YY=DBLE(J)-0.5*DBLE(NY);Y0=Tau_Y/Grad_P
		!IF(YY<=Y0 .AND. YY>=0.0)THEN
		!	UX_L=(0.5/MU_0)*Grad_P*(0.5*DBLE(NY)-Y0)**2
		!ELSEIF(YY>Y0 .AND. YY<=0.5*DBLE(NY))THEN
		!	UX_L=(0.5/MU_0)*Grad_P*((0.5*Dble(Ny))**2-YY**2)-(Tau_Y/Mu_0)*(0.5*DBLE(NY)-YY)
		!ELSEIF(YY>=-Y0 .AND. YY<0)THEN
		!	UX_L=(0.5/MU_0)*Grad_P*(0.5*DBLE(NY)-Y0)**2
		!ELSEIF(YY<-Y0 .AND. YY>=-0.5*DBLE(NY))THEN
		!	UX_L=(0.5/MU_0)*Grad_P*((0.5*Dble(Ny))**2-YY**2)-(Tau_Y/Mu_0)*(0.5*DBLE(NY)-ABS(YY))
		!END IF
		!******************
		!Rho_L=1.0;
		!RHO_L=(F(0,0,J)+F(2,0,J)+F(4,0,J)+2.0*(F(3,0,J)+F(6,0,J)+F(7,0,J)))/(1.0-Ux_L)
        !UX_L=1.0-( F(0,0,J)+F(2,0,J)+F(4,0,J)+2.0*(F(3,0,J)+F(6,0,J)+F(7,0,J)) )/RHO_L
		!F(5,0,J)=F(7,0,J)-0.5*(F(2,0,J)-F(4,0,J))+Inv_six*RHO_L*Ux_L+0.5*RHO_L*Uy_L
        !F(8,0,J)=F(6,0,J)+0.5*(F(2,0,J)-F(4,0,J))+Inv_six*RHO_L*Ux_L-0.5*RHO_L*Uy_L
        !F(1,0,J)=F(3,0,J)+4.0*Inv_six*RHO_L*Ux_L
	!***Right BOUNDARY = Pressure Outlet
		!YY=DBLE(J)-0.5*DBLE(NY);Y0=Tau_Y/Grad_P
		!IF(YY<=Y0 .AND. YY>=0.0)THEN
		!	UX_R=(0.5/MU_0)*Grad_P*(0.5*DBLE(NY)-Y0)**2
		!ELSEIF(YY>Y0 .AND. YY<=0.5*DBLE(NY))THEN
		!	UX_R=(0.5/MU_0)*Grad_P*((0.5*Dble(Ny))**2-YY**2)-(Tau_Y/Mu_0)*(0.5*DBLE(NY)-YY)
		!ELSEIF(YY>=-Y0 .AND. YY<0)THEN
		!	UX_R=(0.5/MU_0)*Grad_P*(0.5*DBLE(NY)-Y0)**2
		!ELSEIF(YY<-Y0 .AND. YY>=-0.5*DBLE(NY))THEN
		!	UX_R=(0.5/MU_0)*Grad_P*((0.5*Dble(Ny))**2-YY**2)-(Tau_Y/Mu_0)*(0.5*DBLE(NY)-ABS(YY))
		!END IF
        !!******************
		!RHO_R=(F(0,NX,J)+F(2,NX,J)+F(4,NX,J)+2.0*(F(1,NX,J)+F(5,NX,J)+F(8,NX,J)))/(1.0+Ux_R)
        !Rho_R=-3.0*DBLE(NX)*Grad_P+RHO_L
		!UX_R=-1.0+( F(0,NX,J)+F(2,NX,J)+F(4,NX,J)+2.0*( F(1,NX,J)+F(5,NX,J)+F(8,NX,J) ) )/RHO_R
		!F(6,Nx,J)=F(8,Nx,J)-0.5*(F(2,Nx,J)-F(4,Nx,J))-Inv_six*RHO_R*Ux_R+0.5*RHO_R*Uy_R
        !F(7,Nx,J)=F(5,Nx,J)+0.5*(F(2,Nx,J)-F(4,Nx,J))-Inv_six*RHO_R*Ux_R-0.5*RHO_R*Uy_R
        !F(3,Nx,J)=F(1,Nx,J)-4.0*Inv_six*RHO_R*Ux_R
	!*********************************************************************
		!***Inlet (Ghost Nodes)
		!F(5,0,J)=F(5,NX,J)
		!F(1,0,J)=F(1,NX,J)
        !F(8,0,J)=F(8,NX,J)
        !***OUTLET PERIODIC (Ghost Nodes)
        !F(6,NX,J)=F(6,0,J)
        !F(3,NX,J)=F(3,0,J)
        !F(7,NX,J)=F(7,0,J)
		!***INLET PERIODIC (Real Nodes) (i dont know why,but ghost node works better for smaller domains)
        !F(5,0,J)=F_POST(5,NX,J-1)
		!F(1,0,J)=F_POST(1,NX,J)
        !F(8,0,J)=F_POST(8,NX,J+1)
        !***OUTLET PERIODIC 
        !F(6,NX,J)=F_POST(6,0,J-1)
        !F(3,NX,J)=F_POST(3,0,J)
        !F(7,NX,J)=F_POST(7,0,J+1)
    END DO

	!IF(0)THEN
	!	CALL COMPUTEFEQNEW(UX(1,1),UY(1,1),RHO(1,1),F(:,0,0))
	!	CALL COMPUTEFEQNEW(UX(1,NY-1),UY(1,NY-1),RHO(1,NY-1),F(:,0,NY))
	!	CALL COMPUTEFEQNEW(UX(NX-1,NY-1),UY(NX-1,NY-1),RHO(NX-1,NY-1),F(:,NX,NY))
	!	CALL COMPUTEFEQNEW(UX(NX-1,1),UY(NX-1,1),RHO(NX-1,1),F(:,NX,0))
	!END IF

	!IF(1)THEN
	I=0;J=0;
	Ux_w=0.0;Uy_w=A_C*OMG_OSC*SIN(OMG_OSC*DBLE(KK));Rho_w=Rho(I,J);
	F(1,I,J)=F_POST(3,I,J)-4.0*Inv_six*(c_lb**2)*Rho_w*( -Ux_W )
	F(5,I,J)=F_POST(7,I,J)-Inv_six*(c_lb**2)*Rho_w*( -Ux_w-Uy_w )
	F(2,I,J)=F_POST(4,I,J)-4.0*Inv_six*(c_lb**2)*Rho_w*( -Uy_W )
	
	I=0;J=Ny;
	Ux_w=0.0;Uy_w=A_C*OMG_OSC*SIN(OMG_OSC*DBLE(KK));Rho_w=Rho(I,J);
	F(1,I,J)=F_POST(3,I,J)-4.0*Inv_six*(c_lb**2)*Rho_w*( -Ux_W )
	F(8,I,J)=F_POST(6,I,J)-Inv_six*(c_lb**2)*Rho_w*( -Ux_w+Uy_w )
	F(4,I,J)=F_POST(2,I,J)-4.0*Inv_six*(c_lb**2)*Rho_w*( Uy_W )
	
	I=Nx;J=Ny;
	Ux_w=0.0;Uy_w=A_C*OMG_OSC*SIN(OMG_OSC*DBLE(KK));Rho_w=Rho(I,J);
	F(3,I,J)=F_POST(1,I,J)-4.0*Inv_six*(c_lb**2)*Rho_w*( Ux_W )
	F(7,I,J)=F_POST(5,I,J)-Inv_six*(c_lb**2)*Rho_w*( Ux_w+Uy_w )
	F(4,I,J)=F_POST(2,I,J)-4.0*Inv_six*(c_lb**2)*Rho_w*( Uy_W )

	I=Nx;J=0;
	Ux_w=0.0;Uy_w=A_C*OMG_OSC*SIN(OMG_OSC*DBLE(KK));Rho_w=Rho(I,J);
	F(2,I,J)=F_POST(4,I,J)-4.0*Inv_six*(c_lb**2)*Rho_w*( -Uy_W )
	F(6,I,J)=F_POST(8,I,J)-Inv_six*(c_lb**2)*Rho_w*( Ux_w-Uy_w )
	F(3,I,J)=F_POST(1,I,J)-4.0*Inv_six*(c_lb**2)*Rho_w*( Ux_W )
	!END IF

    RETURN
END SUBROUTINE BOUNDARYCONDITION
!===================================================================
!COMPUTING MACROSCOPIC DENSITY AND VELOCITIES ON ALL NODES
!===================================================================
SUBROUTINE MACROSCOPICVARIABLES(f,rho,ux,uy)
    use param,only : nx,ny,q,cx,cy
    IMPLICIT NONE
    DOUBLE PRECISION,intent(in) :: f(0:q-1,0:nx,0:ny)
    DOUBLE PRECISION,intent(inout) :: rho(0:nx,0:ny),ux(0:nx,0:ny),uy(0:nx,0:ny)

    DOUBLE PRECISION :: RHOINV
    INTEGER :: I,J,K

    DO J=0,NY !MACRO VELOSCITY ARE CALCULATED FOR ALL INTERNAL NODES
        DO I=0,NX
			RHO(I,J)=F(0,I,J)+F(1,I,J)+F(2,I,J)+F(3,I,J)+F(4,I,J)+F(5,I,J)+F(6,I,J)+F(7,I,J)+F(8,I,J)
			RHOINV=1.0/RHO(I,J)
			UX(I,J)=RHOINV*(F(1,I,J)+F(5,I,J)+F(8,I,J)-(F(6,I,J)+F(3,I,J)+F(7,I,J)))
			UY(I,J)=RHOINV*(F(5,I,J)+F(2,I,J)+F(6,I,J)-(F(7,I,J)+F(4,I,J)+F(8,I,J)))
        END DO
    END DO
    return
END
!===================================================================
!WRITING RESULTS OF SIMULATION
!===================================================================
SUBROUTINE WriteResults(UX,UY,RHO,DOMAIN_IMAGE,STRESS,TSTEP)
    USE PARAM
    IMPLICIT NONE
    DOUBLE PRECISION,INTENT(IN) :: RHO(0:NX,0:NY),UX(0:NX,0:NY),UY(0:NX,0:NY),STRESS(1:3,0:NX,0:NY)
	INTEGER,INTENT(IN) :: DOMAIN_IMAGE(0:NX,0:NY)
	DOUBLE PRECISION :: LSCALE,USCALE
    INTEGER :: I,J,TSTEP
	CHARACTER (LEN=100):: FILENAME
	
	WRITE(FILENAME,*) TSTEP
	FILENAME=ADJUSTL(FILENAME)
	OPEN(201,FILE='02_VelocityField_'//TRIM(FILENAME)//'.plt')

	LSCALE=1.0 !/L_C 
	USCALE=1.0 !/U_C)
    
    WRITE(201,*)'VARIABLES =X/D, Y/D, UX, UY, RHO,Solid_Nodes,TAU_11,TAU_22,TAU_12'
    WRITE(201,*)'ZONE ','I=',NX+1,'J=',NY+1,',','F=BLOCK'
    DO J=0,NY
        WRITE(201,*)(DBLE(I)*LSCALE,I=0,NX)
    END DO
    DO J=0,NY
        WRITE(201,*)((DBLE(J))*LSCALE,I=0,NX)
    END DO
    DO J=0,NY
        WRITE(201,*)(UX(I,J)*USCALE,I=0,NX)
    END DO
    DO J=0,NY
        WRITE(201,*)(UY(I,J)*USCALE,I=0,NX)
    END DO
    DO J=0,NY
        WRITE(201,*)(RHO(I,J),I=0,NX)
    END DO
	DO J=0,NY
        WRITE(201,*)(DBLE(DOMAIN_IMAGE(I,J)),I=0,NX)
    END DO
	DO J=0,NY
        WRITE(201,*)(STRESS(1,I,J),I=0,NX)
    END DO
	DO J=0,NY
        WRITE(201,*)(STRESS(2,I,J),I=0,NX)
    END DO
	DO J=0,NY
        WRITE(201,*)(STRESS(3,I,J),I=0,NX)
    END DO
    CLOSE(201)

END SUBROUTINE
!=====================================================================================
!=====================================================================================
SUBROUTINE WRITE_RESULTSPROFILE(UX,UY,RHO,TSTEP)
    USE PARAM
    IMPLICIT NONE
    DOUBLE PRECISION,INTENT(IN) :: RHO(0:NX,0:NY),UX(0:NX,0:NY),UY(0:NX,0:NY)
    INTEGER :: I,J,TSTEP
	CHARACTER (LEN=100):: FILENAME
	DOUBLE PRECISION :: LSCALE,USCALE,YY,Y0,UX_R
	LSCALE=1.0 !(1.0/L_C)
	USCALE=1.0 !(1.0/U_C)
	
	WRITE(FILENAME,*) TSTEP
	FILENAME=ADJUSTL(FILENAME)
	OPEN(15,FILE='01_VelocityProfiles_'//TRIM(FILENAME)//'.plt')
    WRITE(15,*)'VARIABLES = y/D,ux_6d,ux_0,ux_6d,ux_12d'
    DO J=0,NY
        WRITE(15,*) DBLE(J)*lscale,UX(NX/2-30,J)*USCALE,UX(NX/2,J)*USCALE,UX(NX/2+30,J)*USCALE,UX(NX/2+60,J)*USCALE
    END DO
    CLOSE(15)
end SUBROUTINE WRITE_RESULTSPROFILE

SUBROUTINE RESULTSPROFILE(UX,UY,RHO)
    USE PARAM
    IMPLICIT NONE
    DOUBLE PRECISION,INTENT(IN) :: RHO(0:NX,0:NY),UX(0:NX,0:NY),UY(0:NX,0:NY)
    INTEGER :: I,J
	DOUBLE PRECISION :: LSCALE,USCALE,YY,Y0,UX_R
	LSCALE=(1.0/L_C)
	USCALE=1.0 !(1.0/U_C)

    OPEN(10,FILE='VelocityProfiles.plt')
    WRITE(10,*)'VARIABLES = UX,Y,UY'
    DO J=0,NY
		!YY=DBLE(J)-0.5*DBLE(NY);Y0=Tau_Y/Grad_P
		!IF(YY<=Y0 .AND. YY>=0.0)THEN
		!	UX_R=(0.5/MU_0)*Grad_P*(0.5*DBLE(NY)-Y0)**2
		!ELSEIF(YY>Y0 .AND. YY<=0.5*DBLE(NY))THEN
		!	UX_R=(0.5/MU_0)*Grad_P*((0.5*Dble(Ny))**2-YY**2)-(Tau_Y/Mu_0)*(0.5*DBLE(NY)-YY)
		!ELSEIF(YY>=-Y0 .AND. YY<0)THEN
		!	UX_R=(0.5/MU_0)*Grad_P*(0.5*DBLE(NY)-Y0)**2
		!ELSEIF(YY<-Y0 .AND. YY>=-0.5*DBLE(NY))THEN
		!	UX_R=(0.5/MU_0)*Grad_P*((0.5*Dble(Ny))**2-YY**2)-(Tau_Y/Mu_0)*(0.5*DBLE(NY)-ABS(YY))
		!END IF
		!YY=(DBLE(J)-2.5)/100.0
		!UX_R=4.0*(G_X*1.0E4/(8.0*NU))*(YY-YY*YY)
        WRITE(10,*) UX(NX/2,J)*USCALE,DBLE(J)*LSCALE,UY(J,NY/2)*USCALE
    END DO
    CLOSE(10)
end SUBROUTINE
!====================================================================================
SUBROUTINE WRITEIMAGE(PHIP,TSTEP)
    USE PARAM,ONLY : NX,NY
    IMPLICIT NONE
    DOUBLE PRECISION,INTENT(IN) :: PHIP(0:NX,0:NY)
    INTEGER :: I,J,TSTEP
	CHARACTER (LEN=100):: FILENAME
	
	WRITE(FILENAME,*) TSTEP
	FILENAME=ADJUSTL(FILENAME)
	OPEN(10,FILE='ParticleImage_'//TRIM(FILENAME)//'.plt')

    WRITE(10,*)'VARIABLES =X, Y, PHIP'
    WRITE(10,*)'ZONE ','I=',NX+1,'J=',NY+1,',','F=BLOCK'
    DO J=0,NY
        WRITE(10,*)(FLOAT(I),I=0,NX)
    END DO
    DO J=0,NY
        WRITE(10,*)(FLOAT(J),I=0,NX)
    END DO
    DO J=0,NY
        WRITE(10,*)(PHIP(I,J),I=0,NX)
    END DO
    
    CLOSE(10)

end SUBROUTINE

SUBROUTINE WRITE_DOMAINIMAGE(DOMAIN_IMAGE,TSTEP,FILE_NUMBER)
    USE PARAM,ONLY : NX,NY
    IMPLICIT NONE
    INTEGER,INTENT(IN) :: DOMAIN_IMAGE(0:NX,0:NY)
    INTEGER :: I,J,TSTEP,FILE_NUMBER

    WRITE(FILE_NUMBER,*) 'ZONE T="',TSTEP,'",I=',NX+1,',J=',NY+1,',F=BLOCK'
    DO J=0,NY
       WRITE(FILE_NUMBER,*)(DBLE(I),I=0,NX)
    END DO
    DO J=0,NY
        WRITE(FILE_NUMBER,*)(DBLE(J),I=0,NX)
    END DO
    DO J=0,NY
        WRITE(FILE_NUMBER,*)(DOMAIN_IMAGE(I,J),I=0,NX)
    END DO

end SUBROUTINE

!=====================================================
! COMPUTE ERROR
!=====================================================
SUBROUTINE COMPUTEERROR(ER,UX,UY,U0,V0)
    USE PARAM,ONLY : NX,NY
    IMPLICIT NONE

    DOUBLE PRECISION,INTENT(IN)::UX(0:NX,0:NY),UY(0:NX,0:NY)
    DOUBLE PRECISION,INTENT(INOUT)::U0(0:NX,0:NY),V0(0:NX,0:NY),ER

    INTEGER :: I,J
    DOUBLE PRECISION :: ERR1,ERR2
    ERR1=0.0
    ERR2=0.0
    DO I=0,NX
        DO J=0,NY
            ERR1=ERR1+((UX(I,J)-U0(I,J))*(UX(I,J)-U0(I,J))+(UY(I,J)-V0(I,J))*(UY(I,J)-V0(I,J)))
            ERR2=ERR2+(UX(I,J)*UX(I,J)+UY(I,J)*UY(I,J))
            U0(I,J)=UX(I,J)
            V0(I,J)=UY(I,J)
        END DO
    END DO
    ER=SQRT(ERR1/ERR2)
END SUBROUTINE
!======================================================
! Poisseile Velocity Profile Caculation
!======================================================
DOUBLE PRECISION FUNCTION PoiProfile(Y)
	USE PARAM
	DOUBLE PRECISION :: Y,U,H
	U=U_C
	H=FLOAT(NY)
	PoiProfile=6.0*U*Y*(H-Y)/(H*H)
END FUNCTION

!======================================"THE END"=========================================================================