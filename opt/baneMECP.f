      PROGRAM MECP_Optimizer
      implicit none

C     sobMECP--bane modified, 2025.9.11
C     feedback: banerxmd@gmail.com
C     Input: [inputname].xyz, [inputname].grad1, [inputname].grad2
C     Output: new.xyz (optimized geometry)
C     
C     Compile: 
C       cp MECP.f MECP_temp.f
C       sed -i 's/NUMATOM/$NATOM/g' MECP_temp.f
C       gfortran -O3 -march=native -ffast-math -ffixed-line-length-none baneMECP.f -o MECP.x
C     

      INTEGER Natom, Nx
      parameter (Natom=3)  ! Will be replaced by sed
      parameter (Nx=3*Natom)

      INTEGER Nstep, Conv, AtNum(Natom), i
      DOUBLE PRECISION Ea_1, Ea_2, Eb_1, Eb_2
      DOUBLE PRECISION IH_1(Nx,Nx), IH_2(Nx,Nx)
      DOUBLE PRECISION ParG(Nx), PerpG(Nx)
      DOUBLE PRECISION Ga_1(Nx), Ga_2(Nx)
      DOUBLE PRECISION Gb_1(Nx), Gb_2(Nx)
      DOUBLE PRECISION X_1(Nx), X_2(Nx), X_3(Nx)
      DOUBLE PRECISION G_1(Nx), G_2(Nx)
      DOUBLE PRECISION XBP_prev(Nx), YBP_prev(Nx)
      DOUBLE PRECISION XBP_curr(Nx), YBP_curr(Nx)
      INTEGER IALGO
      
C     Convergence thresholds (now dynamic)
      DOUBLE PRECISION TDE, TDXMax, TDXRMS, TGMax, TGRMS
      
C     Trust radius parameter (now dynamic)
      DOUBLE PRECISION STPMX
      
      CHARACTER*256 inputname, filename, arg
      CHARACTER*16 algorithm
      LOGICAL file_exists, first_step, has_bp_prev
      INTEGER nargs, iargc, natom_actual, iarg
      
C     Set default convergence thresholds
      TDE = 5.d-5
      TDXMax = 4.d-3
      TDXRMS = 2.5d-3
      TGMax = 7.d-4
      TGRMS = 5.d-4
      
C     Set default trust radius
      STPMX = 0.1d0
      algorithm = 'harvey'
      has_bp_prev = .false.
      
C     Parse command line arguments
      nargs = iargc()
      if (nargs .lt. 1) then
          write(*,*) 'Usage: ./MECP.x inputname [options]'
          write(*,*) 'Input files required:'
          write(*,*) '  inputname.xyz   - geometry file'
          write(*,*) '  inputname.grad1 - gradient file for state 1'
          write(*,*) '  inputname.grad2 - gradient file for state 2'
          write(*,*) ''
          write(*,*) 'Optional convergence parameters:'
          write(*,*) '  --tde VALUE     Energy gap threshold (default: 5e-5)'
          write(*,*) '  --tgmax VALUE   Max gradient threshold (default: 7e-4)'
          write(*,*) '  --tgrms VALUE   RMS gradient threshold (default: 5e-4)'
          write(*,*) '  --tdxmax VALUE  Max displacement threshold (default: 4e-3)'
          write(*,*) '  --tdxrms VALUE  RMS displacement threshold (default: 2.5e-3)'
          write(*,*) ''
          write(*,*) 'Optional optimization parameters:'
          write(*,*) '  --stpmx VALUE   Trust radius/max step (default: 0.1)'
          write(*,*) '                  Reduce to 0.02 for oscillating systems'
          write(*,*) '  --algo NAME     Algorithm: harvey | bpupd (default: harvey)'
          write(*,*) ''
          write(*,*) 'Example:'
          write(*,*) '  ./MECP.x mol1 --algo bpupd --tgmax 1e-3 --stpmx 0.02'
          stop
      endif
      
      call getarg(1, inputname)
      
C     Parse optional parameters
      iarg = 2
      do while (iarg .le. nargs)
          call getarg(iarg, arg)
          if (arg .eq. '--tde') then
              if (iarg .eq. nargs) then
                  write(*,*) 'Error: --tde requires a value'
                  stop
              endif
              call getarg(iarg + 1, arg)
              read(arg, *, err=100) TDE
              iarg = iarg + 2
          else if (arg .eq. '--tgmax') then
              if (iarg .eq. nargs) then
                  write(*,*) 'Error: --tgmax requires a value'
                  stop
              endif
              call getarg(iarg + 1, arg)
              read(arg, *, err=101) TGMax
              iarg = iarg + 2
          else if (arg .eq. '--tgrms') then
              if (iarg .eq. nargs) then
                  write(*,*) 'Error: --tgrms requires a value'
                  stop
              endif
              call getarg(iarg + 1, arg)
              read(arg, *, err=102) TGRMS
              iarg = iarg + 2
          else if (arg .eq. '--tdxmax') then
              if (iarg .eq. nargs) then
                  write(*,*) 'Error: --tdxmax requires a value'
                  stop
              endif
              call getarg(iarg + 1, arg)
              read(arg, *, err=103) TDXMax
              iarg = iarg + 2
          else if (arg .eq. '--tdxrms') then
              if (iarg .eq. nargs) then
                  write(*,*) 'Error: --tdxrms requires a value'
                  stop
              endif
              call getarg(iarg + 1, arg)
              read(arg, *, err=104) TDXRMS
              iarg = iarg + 2
          else if (arg .eq. '--stpmx') then
              if (iarg .eq. nargs) then
                  write(*,*) 'Error: --stpmx requires a value'
                  stop
              endif
              call getarg(iarg + 1, arg)
              read(arg, *, err=105) STPMX
              iarg = iarg + 2
          else if (arg .eq. '--algo') then
              if (iarg .eq. nargs) then
                  write(*,*) 'Error: --algo requires a value'
                  stop
              endif
              call getarg(iarg + 1, arg)
              algorithm = arg
              iarg = iarg + 2
          else
              write(*,*) 'Warning: Unknown argument: ', trim(arg)
              iarg = iarg + 1
          endif
      end do
      
C     Print parameters being used
      write(*,*) '=========================================='
      write(*,*) '=== Convergence Thresholds ==='
      write(*,'(A,E12.5)') 'Energy gap (TDE):        ', TDE
      write(*,'(A,E12.5)') 'Max gradient (TGMax):    ', TGMax
      write(*,'(A,E12.5)') 'RMS gradient (TGRMS):    ', TGRMS
      write(*,'(A,E12.5)') 'Max displacement (TDXMax):', TDXMax
      write(*,'(A,E12.5)') 'RMS displacement (TDXRMS):', TDXRMS
      write(*,*) '=== Optimization Parameters ==='
      write(*,'(A,E12.5)') 'Trust radius (STPMX):    ', STPMX
      write(*,'(A,A)') 'Algorithm:               ', trim(algorithm)
      if (STPMX .lt. 0.05d0) then
          write(*,*) 'Note: Using reduced trust radius for better stability'
      endif
      write(*,*) '=========================================='
      write(*,*) ''
      
C     Check if this is first step (no MECP.state file exists)
      inquire(file='MECP.state', exist=file_exists)
      first_step = .not. file_exists
      
C     Read current geometry and gradients
      CALL ReadXYZandGrad(inputname, Natom, Nx, natom_actual,
     &                    AtNum, X_2, Ea_2, Ga_2, Eb_2, Gb_2)
     
      IF (first_step) THEN
C         First step - initialize
          Nstep = 0
          CALL Initialize(Nx, IH_1, Ea_1, Eb_1, Ga_1, Gb_1)
          DO i = 1, Nx
              X_1(i) = X_2(i)
              G_1(i) = 0.d0
          END DO
          DO i = 1, Nx
              XBP_prev(i) = 0.d0
              YBP_prev(i) = 0.d0
          END DO
          has_bp_prev = .false.
      ELSE
C         Read previous state
          CALL ReadState(Natom, Nx, natom_actual, Nstep,
     &                   X_1, IH_1, Ea_1, Eb_1, Ga_1, Gb_1, G_1,
     &                   IALGO, XBP_prev, YBP_prev)
          has_bp_prev = .false.
          if (trim(algorithm) .eq. 'bpupd') then
              if (IALGO .eq. 1) then
                  has_bp_prev = .true.
              endif
          endif
      END IF
      
C     Compute the Effective Gradient according to selected algorithm
      if (trim(algorithm) .eq. 'bpupd') then
          CALL BPUPD_Gradient(Nx, Ea_2, Eb_2, Ga_2, Gb_2,
     &                        has_bp_prev, XBP_prev, YBP_prev,
     &                        ParG, PerpG, G_2, XBP_curr, YBP_curr)
      else
          CALL Effective_Gradient(Nx, Ea_2, Eb_2, Ga_2, Gb_2, 
     &                            ParG, PerpG, G_2)
          DO i = 1, Nx
              XBP_curr(i) = 0.d0
              YBP_curr(i) = 0.d0
          END DO
      endif
     
C     Perform BFGS update with configurable trust radius
      CALL UpdateX(Nx, Nstep, first_step, X_1, X_2, X_3, 
     &             G_1, G_2, IH_1, IH_2, STPMX)
     
C     Test convergence (now with dynamic thresholds)
      CALL TestConvergence(Nx, natom_actual, Nstep, AtNum, 
     &                     Ea_2, Eb_2, X_2, X_3, ParG, PerpG, 
     &                     G_2, Conv, TDE, TDXMax, TDXRMS, TGMax, TGRMS)
     
      IF (Conv .eq. 1) THEN
          write(*,*) 'MECP optimization CONVERGED!'
          write(*,'(A,F18.10)') 'Final energy difference: ',
     &                          ABS(Ea_2 - Eb_2)
C         Write convergence status and criteria
          CALL WriteConvergenceInfo('CONVERGED', Nx, natom_actual, 
     &                              Nstep, Ea_2, Eb_2, X_2, X_3, 
     &                              ParG, PerpG, G_2, TDE, TDXMax, 
     &                              TDXRMS, TGMax, TGRMS)
C         Write final geometry
          CALL WriteXYZ('final.xyz', natom_actual, AtNum, X_3,
     &                  'MECP Converged')
C         Clean up state file
          OPEN(UNIT=99, FILE='MECP.state')
          CLOSE(99, STATUS='DELETE')
      ELSE
C         Write convergence status and criteria
          CALL WriteConvergenceInfo('NOT_CONVERGED', Nx, natom_actual, 
     &                              Nstep, Ea_2, Eb_2, X_2, X_3, 
     &                              ParG, PerpG, G_2, TDE, TDXMax,
     &                              TDXRMS, TGMax, TGRMS)
C         Write new geometry for next iteration
          CALL WriteXYZ('new.xyz', natom_actual, AtNum, X_3,
     &                  'MECP Step')
C         Save current state
          IALGO = 0
          if (trim(algorithm) .eq. 'bpupd') IALGO = 1
          CALL WriteState(Natom, Nx, natom_actual, Nstep,
     &                    X_2, X_3, IH_2, Ea_2, Eb_2, 
     &                    Ga_2, Gb_2, G_2, IALGO,
     &                    XBP_curr, YBP_curr)
      END IF

      goto 999

C     Error handling for argument parsing
100   write(*,*) 'Error: Invalid value for --tde: ', trim(arg)
      stop
101   write(*,*) 'Error: Invalid value for --tgmax: ', trim(arg)
      stop
102   write(*,*) 'Error: Invalid value for --tgrms: ', trim(arg)
      stop
103   write(*,*) 'Error: Invalid value for --tdxmax: ', trim(arg)
      stop
104   write(*,*) 'Error: Invalid value for --tdxrms: ', trim(arg)
      stop
105   write(*,*) 'Error: Invalid value for --stpmx: ', trim(arg)
      stop

999   END


C=====================================================================
      SUBROUTINE WriteConvergenceInfo(status, Nx, Natom_actual, 
     &                                 Nstep, Ea, Eb, X_2, X_3, 
     &                                 ParG, PerpG, G, TDE, TDXMax,
     &                                 TDXRMS, TGMax, TGRMS)
      implicit none
      
      CHARACTER*(*) status
      INTEGER Nx, Natom_actual, Nstep, i
      DOUBLE PRECISION Ea, Eb, X_2(Nx), X_3(Nx), ParG(Nx)
      DOUBLE PRECISION PerpG(Nx), G(Nx)
      DOUBLE PRECISION TDE, TDXMax, TDXRMS, TGMax, TGRMS
      
      CHARACTER*3 flags(5)
      LOGICAL PConv(5)
      DOUBLE PRECISION DeltaX(Nx), DE, DXMax, DXRMS, GMax, GRMS
      DOUBLE PRECISION PpGRMS, PGRMS

C     Calculate convergence criteria
      DE = ABS(Ea - Eb)
      DXMax = 0.d0
      DXRMS = 0.d0
      GMax = 0.d0
      GRMS = 0.d0
      PpGRMS = 0.d0
      PGRMS = 0.d0
      
      DO I = 1, Nx
          DeltaX(i) = X_3(i) - X_2(i)
          IF (ABS(DeltaX(i)) .gt. DXMax) DXMax = ABS(DeltaX(i))
          DXRMS = DXRMS + DeltaX(i)**2
          IF (ABS(G(i)) .gt. Gmax) Gmax = ABS(G(i))
          GRMS = GRMS + G(i)**2
          PpGRMS = PpGRMS + PerpG(i)**2
          PGRMS = PGRMS + ParG(i)**2
      END DO
      
      DXRMS = SQRT(DXRMS / Nx)
      GRMS = SQRT(GRMS / Nx)
      PpGRMS= SQRT(PpGRMS / Nx)
      PGRMS = SQRT(PGRMS / Nx)

C     Check convergence flags
      do i = 1, 5
          flags(i) = " NO"
          PConv(i) = .false.
      end do

      IF (GMax .lt. TGMax) THEN
          PConv(1) = .true.
          flags(1) = "YES"
      END IF
      IF (GRMS .lt. TGRMS) THEN
          PConv(2) = .true.
          flags(2) = "YES"
      END IF
      IF (DXMax .lt. TDXMax) THEN
          PConv(3) = .true.
          flags(3) = "YES"
      END IF
      IF (DXRMS .lt. TDXRMS) THEN
          PConv(4) = .true.
          flags(4) = "YES"
      END IF
      IF (DE .lt. TDE) THEN
          PConv(5) = .true.
          flags(5) = "YES"
      END IF

C     Write to convg.tmp file
      OPEN(UNIT=98, FILE='convg.tmp')
      WRITE(98,'(A)') trim(status)
      WRITE(98,'(A,I4)') 'Step: ', Nstep
      WRITE(98,'(A,F20.10)') 'Energy_State_1: ', Ea
      WRITE(98,'(A,F20.10)') 'Energy_State_2: ', Eb
      WRITE(98,'(A,F20.10)') 'Energy_Gap: ', DE
      WRITE(98,'(A)') 'Convergence_Criteria:'
      WRITE(98,'(A,F12.6,A,F9.6,A,A3)') 
     &    'Max_Gradient: ', GMax, ' (', TGMax, ') ', flags(1)
      WRITE(98,'(A,F12.6,A,F9.6,A,A3)') 
     &    'RMS_Gradient: ', GRMS, ' (', TGRMS, ') ', flags(2)
      WRITE(98,'(A,F12.6,A,F9.6,A,A3)') 
     &    'Max_Displacement: ', DXMax, ' (', TDXMax, ') ', flags(3)
      WRITE(98,'(A,F12.6,A,F9.6,A,A3)') 
     &    'RMS_Displacement: ', DXRMS, ' (', TDXRMS, ') ', flags(4)
      WRITE(98,'(A,F12.6,A,F9.6,A,A3)') 
     &    'Energy_Gap_Conv: ', DE, ' (', TDE, ') ', flags(5)
      WRITE(98,'(A,F12.6)') 'Parallel_Gradient_RMS: ', PGRMS
      WRITE(98,'(A,F12.6)') 'Perpendicular_Gradient_RMS: ', PpGRMS
      CLOSE(98)
      
      RETURN
      END


C=====================================================================
      SUBROUTINE ReadXYZandGrad(basename, MaxNatom, MaxNx, 
     &                          Natom_actual, AtNum, X, 
     &                          Ea, Ga, Eb, Gb)
      implicit none
      
      CHARACTER*(*) basename
      CHARACTER*256 filename
      CHARACTER*2 element
      CHARACTER*256 comment
      INTEGER MaxNatom, MaxNx, Natom_actual, AtNum(MaxNatom)
      INTEGER i, j, k, nx_actual
      DOUBLE PRECISION X(MaxNx), Ea, Ga(MaxNx)
      DOUBLE PRECISION Eb, Gb(MaxNx), dummy
      DOUBLE PRECISION bohr
      PARAMETER (bohr=0.529177d0)
      
C     Read geometry from xyz file
      filename = trim(basename) // '.xyz'
      OPEN(UNIT=10, FILE=filename, STATUS='OLD', ERR=900)
      read(10,*) Natom_actual
      IF (Natom_actual .gt. MaxNatom) THEN
          write(*,*) 'Error: Too many atoms in file'
          stop
      END IF
      nx_actual = 3 * Natom_actual
      READ(10,'(A)') comment  ! Skip comment line
      
      DO i = 1, Natom_actual
          k = 3*(i-1) + 1
          READ(10,*) element, X(k), X(k+1), X(k+2)
          CALL ElementToNumber(element, AtNum(i))
      END DO
      CLOSE(10)
      
C     Read gradient for state 1
      filename = trim(basename) // '.grad1'
      OPEN(UNIT=11, FILE=filename, STATUS='OLD', ERR=901)
      read(11,*) i  ! Should match Natom
      read(11,*) Ea  ! Energy in comment line
      DO i = 1, Natom_actual
          k = 3*(i-1) + 1
          READ(11,*) element, Ga(k), Ga(k+1), Ga(k+2)
C         Assume gradients are already in correct sign (handled by C++ script)
C         Only convert units from Hartree/Bohr to Hartree/Ang
          Ga(k) = Ga(k) / bohr
          Ga(k+1) = Ga(k+1) / bohr
          Ga(k+2) = Ga(k+2) / bohr
      END DO
      CLOSE(11)
      
C     Read gradient for state 2
      filename = trim(basename) // '.grad2'
      OPEN(UNIT=12, FILE=filename, STATUS='OLD', ERR=902)
      read(12,*) i  ! Should match Natom
      read(12,*) Eb  ! Energy in comment line
      DO i = 1, Natom_actual
          k = 3*(i-1) + 1
          READ(12,*) element, Gb(k), Gb(k+1), Gb(k+2)
C         Assume gradients are already in correct sign (handled by C++ script)
C         Only convert units from Hartree/Bohr to Hartree/Ang
          Gb(k) = Gb(k) / bohr
          Gb(k+1) = Gb(k+1) / bohr
          Gb(k+2) = Gb(k+2) / bohr
      END DO
      CLOSE(12)
      
      RETURN
      
900   write(*,*) 'Error: Cannot open file ', trim(filename)
      stop
901   write(*,*) 'Error: Cannot open file ', trim(filename)
      stop
902   write(*,*) 'Error: Cannot open file ', trim(filename)
      stop
      
      END


C=====================================================================
      SUBROUTINE WriteXYZ(filename, Natom_actual, AtNum, X, comment)
      implicit none
      
      CHARACTER*(*) filename, comment
      INTEGER Natom_actual, AtNum(*), i, k
      DOUBLE PRECISION X(*)
      CHARACTER*2 element
      
      OPEN(UNIT=20, FILE=filename)
      WRITE(20,*) Natom_actual
      WRITE(20,'(A)') trim(comment)
      
      DO i = 1, Natom_actual
          k = 3*(i-1) + 1
          CALL NumberToElement(AtNum(i), element)
          WRITE(20,'(A2,3F16.8)') element, X(k), X(k+1), X(k+2)
      END DO
      
      CLOSE(20)
      RETURN
      END


C=====================================================================
      SUBROUTINE ReadState(MaxNatom, MaxNx, Natom_actual, Nstep,
     &                     X_1, IH, Ea, Eb, Ga, Gb, G,
     &                     IALGO, XBP_prev, YBP_prev)
      implicit none
      
      INTEGER MaxNatom, MaxNx, Natom_actual, Nstep, i, j
      INTEGER nx_actual
      DOUBLE PRECISION X_1(MaxNx), IH(MaxNx,MaxNx)
      DOUBLE PRECISION Ea, Eb, Ga(MaxNx), Gb(MaxNx), G(MaxNx)
      INTEGER IALGO
      DOUBLE PRECISION XBP_prev(MaxNx), YBP_prev(MaxNx)
      
      nx_actual = 3 * Natom_actual
      
      OPEN(UNIT=30, FILE='MECP.state', STATUS='OLD')
      READ(30,*) Natom_actual
      READ(30,*) Nstep
      
C     Read previous geometry
      DO i = 1, nx_actual
          READ(30,*) X_1(i)
      END DO
      
C     Read energies
      READ(30,*) Ea, Eb
      
C     Read gradients
      DO i = 1, nx_actual
          READ(30,*) Ga(i)
      END DO
      DO i = 1, nx_actual
          READ(30,*) Gb(i)
      END DO
      DO i = 1, nx_actual
          READ(30,*) G(i)
      END DO
      
C     Read inverse Hessian
      DO i = 1, nx_actual
          DO j = 1, nx_actual
              READ(30,*) IH(i,j)
          END DO
      END DO

C     Try to read optional BPUPD state
      IALGO = 0
  910 CONTINUE
      READ(30,*,END=920,ERR=920) IALGO
      DO i = 1, nx_actual
          READ(30,*,END=920,ERR=920) XBP_prev(i)
      END DO
      DO i = 1, nx_actual
          READ(30,*,END=920,ERR=920) YBP_prev(i)
      END DO
      GOTO 930
  920 CONTINUE
      IALGO = 0
      DO i = 1, nx_actual
          XBP_prev(i) = 0.d0
          YBP_prev(i) = 0.d0
      END DO
  930 CONTINUE
      
      CLOSE(30)
      RETURN
      END


C=====================================================================
      SUBROUTINE WriteState(MaxNatom, MaxNx, Natom_actual, Nstep,
     &                      X_2, X_3, IH, Ea, Eb, Ga, Gb, G,
     &                      IALGO, XBP_prev, YBP_prev)
      implicit none
      
      INTEGER MaxNatom, MaxNx, Natom_actual, Nstep, i, j
      INTEGER nx_actual
      DOUBLE PRECISION X_2(MaxNx), X_3(MaxNx), IH(MaxNx,MaxNx)
      DOUBLE PRECISION Ea, Eb, Ga(MaxNx), Gb(MaxNx), G(MaxNx)
      INTEGER IALGO
      DOUBLE PRECISION XBP_prev(MaxNx), YBP_prev(MaxNx)
      
      nx_actual = 3 * Natom_actual
      
      OPEN(UNIT=31, FILE='MECP.state')
      WRITE(31,*) Natom_actual
      WRITE(31,*) Nstep
      
C     Write current geometry (will be "previous" in next step)
      DO i = 1, nx_actual
          WRITE(31,'(F20.12)') X_2(i)
      END DO
      
C     Write energies
      WRITE(31,'(2F20.12)') Ea, Eb
      
C     Write gradients
      DO i = 1, nx_actual
          WRITE(31,'(F20.12)') Ga(i)
      END DO
      DO i = 1, nx_actual
          WRITE(31,'(F20.12)') Gb(i)
      END DO
      DO i = 1, nx_actual
          WRITE(31,'(F20.12)') G(i)
      END DO
      
C     Write inverse Hessian
      DO i = 1, nx_actual
          DO j = 1, nx_actual
              WRITE(31,'(F20.12)') IH(i,j)
          END DO
      END DO

C     Write optional BPUPD state
      WRITE(31,*) IALGO
      DO i = 1, nx_actual
          WRITE(31,'(F20.12)') XBP_prev(i)
      END DO
      DO i = 1, nx_actual
          WRITE(31,'(F20.12)') YBP_prev(i)
      END DO
      
      CLOSE(31)
      RETURN
      END


C=====================================================================
      SUBROUTINE ElementToNumber(element, number)
      implicit none
      
      CHARACTER*2 element, normalized_element
      INTEGER number
      CHARACTER*2 symbols(118)
      INTEGER i, ic1, ic2
      
      DATA symbols /'H ','He','Li','Be','B ','C ','N ','O ','F ',
     & 'Ne','Na','Mg','Al','Si','P ','S ','Cl','Ar','K ','Ca',
     & 'Sc','Ti','V ','Cr','Mn','Fe','Co','Ni','Cu','Zn',
     & 'Ga','Ge','As','Se','Br','Kr','Rb','Sr','Y ','Zr',
     & 'Nb','Mo','Tc','Ru','Rh','Pd','Ag','Cd','In','Sn',
     & 'Sb','Te','I ','Xe','Cs','Ba','La','Ce','Pr','Nd',
     & 'Pm','Sm','Eu','Gd','Tb','Dy','Ho','Er','Tm','Yb',
     & 'Lu','Hf','Ta','W ','Re','Os','Ir','Pt','Au','Hg',
     & 'Tl','Pb','Bi','Po','At','Rn','Fr','Ra','Ac','Th',
     & 'Pa','U ','Np','Pu','Am','Cm','Bk','Cf','Es','Fm',
     & 'Md','No','Lr','Rf','Db','Sg','Bh','Hs','Mt','Ds',
     & 'Rg','Cn','Nh','Fl','Mc','Lv','Ts','Og'/
      
C     Normalize element symbol: first letter uppercase, second lowercase
      normalized_element = element
      
C     Convert first character to uppercase
      ic1 = ICHAR(normalized_element(1:1))
      if (ic1 .ge. ICHAR('a') .and. ic1 .le. ICHAR('z')) then
          normalized_element(1:1) = CHAR(ic1 - ICHAR('a') + ICHAR('A'))
      endif
      
C     Convert second character to lowercase (if it exists and is not space)
      if (normalized_element(2:2) .ne. ' ') then
          ic2 = ICHAR(normalized_element(2:2))
          if (ic2 .ge. ICHAR('A') .and. ic2 .le. ICHAR('Z')) then
              normalized_element(2:2) = CHAR(ic2 - ICHAR('A') + ICHAR('a'))
          endif
      endif
      
      number = 0
      DO i = 1, 118
          IF (normalized_element .eq. symbols(i)) THEN
              number = i
              RETURN
          END IF
      END DO
      
C     Error: Unknown element
      IF (number .eq. 0) THEN
          write(*,*) 'Error: Unknown element symbol: ', element
          write(*,*) 'Normalized to: ', normalized_element
          write(*,*) 'Check your input geometry file for typos'
          stop
      END IF
      
      RETURN
      END


C=====================================================================
      SUBROUTINE NumberToElement(number, element)
      implicit none
      
      INTEGER number
      CHARACTER*2 element
      CHARACTER*2 symbols(118)
      
      DATA symbols /'H ','He','Li','Be','B ','C ','N ','O ','F ',
     & 'Ne','Na','Mg','Al','Si','P ','S ','Cl','Ar','K ','Ca',
     & 'Sc','Ti','V ','Cr','Mn','Fe','Co','Ni','Cu','Zn',
     & 'Ga','Ge','As','Se','Br','Kr','Rb','Sr','Y ','Zr',
     & 'Nb','Mo','Tc','Ru','Rh','Pd','Ag','Cd','In','Sn',
     & 'Sb','Te','I ','Xe','Cs','Ba','La','Ce','Pr','Nd',
     & 'Pm','Sm','Eu','Gd','Tb','Dy','Ho','Er','Tm','Yb',
     & 'Lu','Hf','Ta','W ','Re','Os','Ir','Pt','Au','Hg',
     & 'Tl','Pb','Bi','Po','At','Rn','Fr','Ra','Ac','Th',
     & 'Pa','U ','Np','Pu','Am','Cm','Bk','Cf','Es','Fm',
     & 'Md','No','Lr','Rf','Db','Sg','Bh','Hs','Mt','Ds',
     & 'Rg','Cn','Nh','Fl','Mc','Lv','Ts','Og'/
      
      IF (number .ge. 1 .and. number .le. 118) THEN
          element = symbols(number)
      ELSE
          element = 'XX'
      END IF
      
      RETURN
      END


C=====================================================================
C                          Effective_Gradient 
C=====================================================================
       SUBROUTINE Effective_Gradient(N,Ea,Eb,Ga,Gb,ParG,PerpG,G)
       implicit none
       
C      Computes the parallel and perpendicular compenents of the Effective Gradient,
C      As well as the effective gradient itself.
       
       integer n, i
       double precision Ea, Eb, Ga(N), Gb(N), G(N), ParG(N), PerpG(N), npg, pp
       double precision facPP, facP
       parameter (facPP=140.d0,facP=1.d0)
C  These factors are only really important for the first step
C  The "difference gradient" is typically ca. 0.075 Hartree/Bohr.
C      i.e. 0.14 Hartree/Angstrom.
C  Assuming this is constant, this gives a Hessian term for the func (Ea-Eb)**2
C     of ca. 0.01 Hartree**2 / Angstrom**2  (0.14**2 / 2)
C  The Hessian term along "normal" coordinates is empirically about 1.4 Hartree / Angstrom**2
C  Using facPP ~ 140 /Hartree means that along this coordinate, too, the Hessian is about right.
       
       npg = 0.d0
       pp = 0.d0
       do i = 1, n
           PerpG(i) = Ga(i) - Gb(i)
           npg = npg + PerpG(i)**2
           pp = pp + Ga(i) * PerpG(i)
       end do
       npg = sqrt(npg)
       pp = pp / npg
       do i = 1, n
           ParG(i) = Ga(i) - PerpG(i) / npg * pp
           G(i) = (Ea - Eb) * facPP * PerpG(i) + facP * ParG(i)
       end do
       return
       
       END


C=====================================================================
C                          BPUPD_Gradient
C=====================================================================
      SUBROUTINE BPUPD_Gradient(N,Ea,Eb,Ga,Gb,HasPrev,
     &                         Xprev,Yprev,ParG,PerpG,G,Xout,Yout)
      implicit none
      INTEGER N, i, j
      DOUBLE PRECISION Ea, Eb, Ga(N), Gb(N), ParG(N), PerpG(N), G(N)
      DOUBLE PRECISION Xprev(N), Yprev(N), Xout(N), Yout(N)
      LOGICAL HasPrev
      DOUBLE PRECISION gd_norm, dotxx, dotyx, ax, ay, denom
      DOUBLE PRECISION gm_dot_x, gm_dot_y, de
      DOUBLE PRECISION gd_i, gm_i, x_i, y_i, tmp
      DOUBLE PRECISION gm(N), proj_gm(N)
      DOUBLE PRECISION eps
      PARAMETER (eps=1.0d-12)
      
C     Compute difference and mean gradients
      gd_norm = 0.d0
      DO i = 1, N
          PerpG(i) = Gb(i) - Ga(i)  ! g2 - g1
          gm(i) = 0.5d0 * (Ga(i) + Gb(i))  ! mean gradient
          gd_norm = gd_norm + PerpG(i)*PerpG(i)
      END DO
      gd_norm = SQRT(gd_norm)
      
      IF (gd_norm .lt. eps) THEN
C         Fallback to Harvey method if degenerate
          CALL Effective_Gradient(N, Ea, Eb, Ga, Gb, ParG, PerpG, G)
          DO i = 1, N
              Xout(i) = 0.d0
              Yout(i) = 0.d0
          END DO
          RETURN
      END IF
      
C     x = normalized gradient difference
      DO i = 1, N
          Xout(i) = PerpG(i) / gd_norm
      END DO
      
C     y vector update
      IF (.not. HasPrev) THEN
C         First step: y = (gm - (x·gm)x) / |gm - (x·gm)x|
          gm_dot_x = 0.d0
          DO i = 1, N
              gm_dot_x = gm_dot_x + Xout(i) * gm(i)
          END DO
          
          tmp = 0.d0
          DO i = 1, N
              Yout(i) = gm(i) - gm_dot_x * Xout(i)
              tmp = tmp + Yout(i)*Yout(i)
          END DO
          
          IF (tmp .gt. eps) THEN
              tmp = 1.d0 / SQRT(tmp)
              DO i = 1, N
                  Yout(i) = Yout(i) * tmp
              END DO
          ELSE
              DO i = 1, N
                  Yout(i) = 0.d0
              END DO
          END IF
      ELSE
C         后续步骤的y更新：基于前一步的x和y向量
          dotyx = 0.d0
          dotxx = 0.d0
          DO i = 1, N
              dotyx = dotyx + Yprev(i) * Xout(i)
              dotxx = dotxx + Xprev(i) * Xout(i)
          END DO
          
          denom = dotyx*dotyx + dotxx*dotxx
          IF (denom .gt. eps) THEN
              denom = SQRT(denom)
              DO i = 1, N
                  Yout(i) = (dotyx * Xprev(i) - dotxx * Yprev(i)) / denom
              END DO
          ELSE
              DO i = 1, N
                  Yout(i) = 0.d0
              END DO
          END IF
      END IF
      
C     计算投影后的平均梯度: proj_gm = (I - xx^T - yy^T) * gm
C     这等价于: proj_gm = gm - (x·gm)x - (y·gm)y
      gm_dot_x = 0.d0
      gm_dot_y = 0.d0
      DO i = 1, N
          gm_dot_x = gm_dot_x + Xout(i) * gm(i)
          gm_dot_y = gm_dot_y + Yout(i) * gm(i)
      END DO
      
      DO i = 1, N
          proj_gm(i) = gm(i) - gm_dot_x * Xout(i) - gm_dot_y * Yout(i)
      END DO
      
C     构造有效梯度: G = 2*(E2-E1)*x + proj_gm
      de = Eb - Ea
      DO i = 1, N
C         平行梯度分量（用于输出和分析）
          ParG(i) = proj_gm(i)
C         有效梯度
          G(i) = 2.d0 * de * Xout(i) + proj_gm(i)
      END DO
      
      RETURN
      END

C=====================================================================
C                          Initialize
C=====================================================================
      SUBROUTINE Initialize(N,H,Ea,Eb,Ga,Gb)
      implicit none
       
      INTEGER i, j, n
      DOUBLE PRECISION H(N,N), Ea, Eb, Ga(N), Gb(N)
       
      Ea = 0.d0
      Eb = 0.d0
      do i = 1, n
          Ga(i) = 0.d0
          Gb(i) = 0.d0
          H(i,i) = .7d0
C Inverse Hessian values of .7 correspond to Hessians of 1.4 Hartree/Angstrom**2 - about right
          do j = i + 1, n
              H(i,j) = 0.d0
              H(j,i) = 0.d0
          end do
      end do
      return

      END


C=====================================================================
C                           UpdateX
C=====================================================================
       SUBROUTINE UpdateX(N,Nstep,FirstStep,X_1,X_2,X_3,G_1,G_2,
     &                   HI_1,HI_2,STPMX)
       implicit none
       
       ! Specially Adapted BFGS routine from Numerical Recipes
       ! Modified to accept STPMX as parameter instead of using fixed value
       
       integer i, j, k, n, Nstep
       logical FirstStep
       double precision X_1(N), X_2(N), G_1(N), G_2(N)
       double precision HI_1(N,N), X_3(N), HI_2(N,N)
       double precision STPMX  ! Now passed as parameter
       
       double precision stpmax, DelG(N), HDelG(N), ChgeX(N), DelX(N), w(N), fac,
     & fad, fae, sumdg, sumdx, stpl, lgstst
       
       stpmax = STPMX * N
       IF (FirstStep) THEN
           DO i = 1, n
               ChgeX(i) = -.7d0 * G_2(i)
               DO j = 1, n
                   HI_2(i,j) = HI_1(i,j)
               end do
           end do
       ELSE
           DO i = 1, n
               DelG(i) = G_2(i) - G_1(i)
               DelX(i) = X_2(i) - X_1(i)
           end do
           do i = 1, n
               HDelG(i) = 0.d0
               do j = 1, n
                   hdelg(i) = hdelg(i) + HI_1(i,j) * DelG(j)
               end do
           end do
           fac = 0.d0
           fae = 0.d0
           sumdg = 0.d0
           sumdx = 0.d0
           do i = 1, n
               fac = fac + delg(i) * delx(i)
               fae = fae + delg(i) * hdelg(i)
               sumdg = sumdg + delg(i)**2
               sumdx = sumdx + delx(i)**2
           end do
           fac = 1.d0 / fac
           fad = 1.d0 / fae
           do i = 1, n
               w(i) = fac * DelX(i) - fad * HDelG(i)
           end do
           DO I = 1, N
               do j = 1, n
                   HI_2(i,j) = HI_1(i,j) + fac * delx(i) * delx(j) -
     &                fad * HDelG(I) * HDelG(j) + fae * w(i) * w(j)
               end do
           end do
           do i = 1, n
               ChgeX(i) = 0.d0
               do j = 1, n
                   ChgeX(i) = ChgeX(i) - HI_2(i,j) * G_2(j)
               end do
           end do
       END IF
       
       fac = 0.d0
       do i = 1, n
           fac = fac + ChgeX(i)**2
       end do
       stpl = SQRT(fac)
       IF (stpl .gt. stpmax) THEN
           do i = 1, n
               ChgeX(i) = ChgeX(i) / stpl * stpmax
           end do
       END IF
       lgstst = 0.d0
       do i = 1, n
            IF (ABS(ChgeX(i)) .gt. lgstst) lgstst = ABS(ChgeX(i))
       end do
       IF (lgstst .gt. STPMX) THEN
           do i = 1, n
               ChgeX(i) = ChgeX(i) / lgstst * STPMX
           end do
       END IF
       do i = 1, n
           X_3(i) = X_2(i) + ChgeX(i)
       end do
       
       END


C=====================================================================
      SUBROUTINE TestConvergence(Nx,Natom_actual,Nstep,AtNum,Ea,Eb,
     &                           X_2,X_3,ParG,PerpG,G,Conv,TDE,TDXMax,
     &                           TDXRMS,TGMax,TGRMS)
      implicit none

      INTEGER Nx, Natom_actual, AtNum(*), Nstep, Conv, i, j, k
      DOUBLE PRECISION Ea, Eb, X_2(Nx), ParG(Nx), PerpG(Nx)
      DOUBLE PRECISION X_3(Nx), G(Nx)
      DOUBLE PRECISION TDE, TDXMax, TDXRMS, TGMax, TGRMS

      CHARACTER*3 flags(5)
      LOGICAL PConv(5)
      DOUBLE PRECISION DeltaX(Nx), DE, DXMax, DXRMS, GMax, GRMS
      DOUBLE PRECISION PpGRMS, PGRMS
      CHARACTER*2 element

      DE = ABS(Ea - Eb)
      DXMax = 0.d0
      DXRMS = 0.d0
      GMax = 0.d0
      GRMS = 0.d0
      PpGRMS = 0.d0
      PGRMS = 0.d0
      
      DO I = 1, Nx
          DeltaX(i) = X_3(i) - X_2(i)
          IF (ABS(DeltaX(i)) .gt. DXMax) DXMax = ABS(DeltaX(i))
          DXRMS = DXRMS + DeltaX(i)**2
          IF (ABS(G(i)) .gt. Gmax) Gmax = ABS(G(i))
          GRMS = GRMS + G(i)**2
          PpGRMS = PpGRMS + PerpG(i)**2
          PGRMS = PGRMS + ParG(i)**2
      END DO
      
      DXRMS = SQRT(DXRMS / Nx)
      GRMS = SQRT(GRMS / Nx)
      PpGRMS= SQRT(PpGRMS / Nx)
      PGRMS = SQRT(PGRMS / Nx)

      Conv = 0
      do i = 1, 5
          flags(i) = " NO"
          PConv(i) = .false.
      end do

      IF (GMax .lt. TGMax) THEN
          PConv(1) = .true.
          flags(1) = "YES"
      END IF
      IF (GRMS .lt. TGRMS) THEN
          PConv(2) = .true.
          flags(2) = "YES"
      END IF
      IF (DXMax .lt. TDXMax) THEN
          PConv(3) = .true.
          flags(3) = "YES"
      END IF
      IF (DXRMS .lt. TDXRMS) THEN
          PConv(4) = .true.
          flags(4) = "YES"
      END IF
      IF (DE .lt. TDE) THEN
          PConv(5) = .true.
          flags(5) = "YES"
      END IF

      Nstep = Nstep + 1 
      IF (PConv(1) .and. PConv(2) .and. PConv(3) .and. 
     &    PConv(4) .and. PConv(5)) THEN
          Conv = 1
      ELSE
          Conv = 0
      END IF

C     Write initial geometry for first step
      IF (Nstep .eq. 1) THEN
          write(*,*)
          write(*,'(A)') '  Initial Geometry:'
          DO i = 1, Natom_actual
              k = 3*(i-1) + 1
              CALL NumberToElement(AtNum(i), element)
              WRITE(*,'(A,A2,3F16.8)') '    ', element, X_2(k), 
     &                                 X_2(k+1), X_2(k+2)
          END DO
          write(*,*)
      END IF

C     Write convergence info to screen
      write(*,*)
      write(*,'(A,I4,A)') '  ===== MECP Optimization Step', 
     &                     Nstep, ' ====='
      write(*,'(A,F20.10)') '    Energy of State 1:  ', Ea
      write(*,'(A,F20.10)') '    Energy of State 2:  ', Eb
      write(*,'(A,F20.10)') '    Energy Difference:  ', DE
      write(*,*)
      write(*,'(A)') '  Convergence Check:'
      write(*,'(A,F12.6,A,F9.6,A,A3)') 
     &    '    Max Gradient:     ', GMax, ' (', TGMax, ')  ', flags(1)
      write(*,'(A,F12.6,A,F9.6,A,A3)') 
     &    '    RMS Gradient:     ', GRMS, ' (', TGRMS, ')  ', flags(2)
      write(*,'(A,F12.6,A,F9.6,A,A3)') 
     &    '    Max Displacement: ', DXMax, ' (', TDXMax, ')  ', flags(3)
      write(*,'(A,F12.6,A,F9.6,A,A3)') 
     &    '    RMS Displacement: ', DXRMS, ' (', TDXRMS, ')  ', flags(4)
      write(*,'(A,F12.6,A,F9.6,A,A3)') 
     &    '    Energy Gap:       ', DE, ' (', TDE, ')  ', flags(5)
      write(*,*)
      
C     Write gradient details
      write(*,'(A)') '  Overall Effective Gradient:'
      DO i = 1, Natom_actual
          k = 3*(i-1) + 1
          WRITE(*,'(A,I3,3F16.8)') '  ', i, G(k), G(k+1), G(k+2)
      END DO
      write(*,*)
      
      write(*,'(A,F12.6,A)') '  Difference Gradient (RMS * DE: ',
     &                        PpGRMS * DE, '):'
      DO i = 1, Natom_actual
          k = 3*(i-1) + 1
          WRITE(*,'(A,I3,3F16.8)') '  ', i, PerpG(k), PerpG(k+1), 
     &                             PerpG(k+2)
      END DO
      write(*,*)
      
      write(*,'(A,F12.6,A)') '  Parallel Gradient   (RMS     : ', PGRMS, '):'
      DO i = 1, Natom_actual
          k = 3*(i-1) + 1
          WRITE(*,'(A,I3,3F16.8)') '  ', i, ParG(k), ParG(k+1), 
     &                             ParG(k+2)
      END DO
      write(*,*)

      IF (Conv .ne. 1) THEN
          write(*,'(A,I4)') '  Geometry at Step', Nstep
          DO i = 1, Natom_actual
              k = 3*(i-1) + 1
              CALL NumberToElement(AtNum(i), element)
              WRITE(*,'(A,A2,3F16.8)') '    ', element, X_3(k), 
     &                                 X_3(k+1), X_3(k+2)
      END DO
          write(*,*)
      END IF

      return
      END