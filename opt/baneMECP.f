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

      INTEGER MAXDIIS
      parameter (MAXDIIS=8)

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

C     Lagrangian-constraint / penalty-function parameters (XMECP-style)
C     Default values follow XMECP inp_proc.py (pf_alpha=0.02, pf_sigma=3.50)
      DOUBLE PRECISION PF_ALPHA, PF_SIGMA

C     PF-specific convergence criteria (XMECP-style)
C     pf_thresh=[1e-6, 5e-3] in XMECP.
C       PF_TSTEP: objective function change threshold
C       PF_TGRAD: penalty-function gradient threshold
      DOUBLE PRECISION PF_TSTEP, PF_TGRAD
      
      CHARACTER*256 inputname, filename, arg
      CHARACTER*16 algorithm, opt_method, hupd_method
      LOGICAL file_exists, first_step, has_bp_prev
      INTEGER nargs, iargc, natom_actual, iarg, stepnum, NDIIS
      DOUBLE PRECISION EOBJ
C     Set default convergence thresholds
      TDE = 5.d-5
      TDXMax = 4.d-3
      TDXRMS = 2.5d-3
      TGMax = 7.d-4
      TGRMS = 5.d-4
      
C     Set default trust radius
      STPMX = 0.1d0
      algorithm = 'harvey'
      opt_method = 'bfgs'
      hupd_method = 'bfgs'
      NDIIS = 4
      has_bp_prev = .false.

C     Set default Lagrangian/penalty parameters
      PF_ALPHA = 2.0d-2
      PF_SIGMA = 3.5d0

C     Set default PF-specific convergence thresholds (XMECP-style)
      PF_TSTEP = 1.0d-6
      PF_TGRAD = 5.0d-3
      
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
          write(*,*) '  --algo NAME     Algorithm: harvey | bpupd | lagrange (pf)'
          write(*,*) '  --opt NAME      Opt method: bfgs | gdiis | gediis (default: bfgs)'
          write(*,*) '  --hupd NAME     Hessian update: bfgs | psb | sr1 | bofill (default: bfgs)'
          write(*,*) '  --ndiis N       DIIS history size (default: 4, max: 8)'
          write(*,*) ''
          write(*,*) 'Optional Lagrangian constraint parameters (for lagrange/pf):'
          write(*,*) '  --pf_alpha VALUE  Alpha parameter (default: 0.02)'
          write(*,*) '  --pf_sigma VALUE  Sigma parameter (default: 3.5)'
          write(*,*) '  --pf_tstep VALUE  PF step threshold (default: 1e-6)'
          write(*,*) '  --pf_tgrad VALUE  PF gradient threshold (default: 5e-3)'
          write(*,*) ''
          write(*,*) 'Example:'
          write(*,*) '  ./MECP.x mol1 --algo bpupd --tgmax 1e-3 --stpmx 0.02'
          write(*,*) '  ./MECP.x mol1 --algo lagrange --pf_alpha 0.02 --pf_sigma 3.5 --pf_tstep 1e-6 --pf_tgrad 5e-3'
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
          else if (arg .eq. '--pf_alpha') then
              if (iarg .eq. nargs) then
                  write(*,*) 'Error: --pf_alpha requires a value'
                  stop
              endif
              call getarg(iarg + 1, arg)
              read(arg, *, err=106) PF_ALPHA
              iarg = iarg + 2
          else if (arg .eq. '--pf_sigma') then
              if (iarg .eq. nargs) then
                  write(*,*) 'Error: --pf_sigma requires a value'
                  stop
              endif
              call getarg(iarg + 1, arg)
              read(arg, *, err=107) PF_SIGMA
              iarg = iarg + 2
          else if (arg .eq. '--pf_tstep') then
              if (iarg .eq. nargs) then
                  write(*,*) 'Error: --pf_tstep requires a value'
                  stop
              endif
              call getarg(iarg + 1, arg)
              read(arg, *, err=109) PF_TSTEP
              iarg = iarg + 2
          else if (arg .eq. '--pf_tgrad') then
              if (iarg .eq. nargs) then
                  write(*,*) 'Error: --pf_tgrad requires a value'
                  stop
              endif
              call getarg(iarg + 1, arg)
              read(arg, *, err=110) PF_TGRAD
              iarg = iarg + 2
          else if (arg .eq. '--algo') then
              if (iarg .eq. nargs) then
                  write(*,*) 'Error: --algo requires a value'
                  stop
              endif
              call getarg(iarg + 1, arg)
              algorithm = arg
              iarg = iarg + 2
          else if (arg .eq. '--opt' .or. arg .eq. '--opt_method') then
              if (iarg .eq. nargs) then
                  write(*,*) 'Error: --opt requires a value'
                  stop
              endif
              call getarg(iarg + 1, arg)
              opt_method = arg
              iarg = iarg + 2
                    else if (arg .eq. '--hupd' .or.
     &             arg .eq. '--hess' .or.
     &             arg .eq. '--hessupd') then
              if (iarg .eq. nargs) then
                  write(*,*) 'Error: --hupd requires a value'
                  stop
              endif
              call getarg(iarg + 1, arg)
              hupd_method = arg
              iarg = iarg + 2
          else if (arg .eq. '--ndiis') then

              if (iarg .eq. nargs) then
                  write(*,*) 'Error: --ndiis requires a value'
                  stop
              endif
              call getarg(iarg + 1, arg)
              read(arg, *, err=108) NDIIS
              iarg = iarg + 2
                    else
              write(*,*) 'Warning: Unknown argument: ', trim(arg)
              iarg = iarg + 1
          endif
      end do

C     Clamp DIIS subspace size
      if (NDIIS .lt. 1) NDIIS = 1
      if (NDIIS .gt. MAXDIIS) then
          write(*,*) 'Warning: NDIIS too large, using MAXDIIS =', MAXDIIS
          NDIIS = MAXDIIS
      endif

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
      write(*,'(A,A)') 'Opt method:              ', trim(opt_method)
      write(*,'(A,A)') 'Hessian update:          ', trim(hupd_method)
      write(*,'(A,I8)') 'DIIS history (NDIIS):    ', NDIIS
      if (trim(algorithm) .eq. 'lagrange' .or.
     &    trim(algorithm) .eq. 'pf') then
          write(*,*) '=== Lagrangian Constraint Parameters ==='
          write(*,'(A,E12.5)') 'pf_alpha (alpha):        ', PF_ALPHA
          write(*,'(A,E12.5)') 'pf_sigma (sigma):        ', PF_SIGMA
          write(*,*) '=== PF Convergence (XMECP-style) ==='
          write(*,'(A,E12.5)') 'pf_tstep:                ', PF_TSTEP
          write(*,'(A,E12.5)') 'pf_tgrad:                ', PF_TGRAD
      endif
      if (STPMX .lt. 0.05d0) then
          write(*,*) 'Note: Using reduced trust radius for better stability'
      endif
      write(*,*) '=========================================='
      write(*,*) ''

C     Validate algorithm keyword
      IF (trim(algorithm) .ne. 'harvey' .and.
     &    trim(algorithm) .ne. 'bpupd' .and.
     &    trim(algorithm) .ne. 'lagrange' .and.
     &    trim(algorithm) .ne. 'pf') THEN
          write(*,*) 'Error: Unknown algorithm: ', trim(algorithm)
          write(*,*) 'Supported algorithms: harvey, bpupd, lagrange (pf)'
          STOP
      ENDIF

C     Validate optimization method keyword
      IF (trim(opt_method) .ne. 'bfgs' .and.
     &    trim(opt_method) .ne. 'gdiis' .and.
     &    trim(opt_method) .ne. 'gediis') THEN
          write(*,*) 'Error: Unknown opt_method: ', trim(opt_method)
          write(*,*) 'Supported opt methods: bfgs, gdiis, gediis'
          STOP
      ENDIF
C     Validate Hessian update method keyword
      IF (trim(hupd_method) .ne. 'bfgs' .and.
     &    trim(hupd_method) .ne. 'psb' .and.
     &    trim(hupd_method) .ne. 'sr1' .and.
     &    trim(hupd_method) .ne. 'bofill') THEN
          write(*,*) 'Error: Unknown hessian update: ', trim(hupd_method)
          write(*,*) 'Supported hessian updates: bfgs, psb, sr1, bofill'
          STOP
      ENDIF


C     Check if this is first step (no MECP.state file exists)
      inquire(file='MECP.state', exist=file_exists)
      first_step = .not. file_exists

C     If starting fresh, remove old DIIS history file
      if (first_step) then
          inquire(file='MECP.diis', exist=file_exists)
          if (file_exists) then
              OPEN(UNIT=97, FILE='MECP.diis')
              CLOSE(97, STATUS='DELETE')
          endif
      endif
      
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
      else if (trim(algorithm) .eq. 'lagrange' .or.
     &         trim(algorithm) .eq. 'pf') then
          CALL LAGRANGE_Gradient(Nx, Ea_2, Eb_2, Ga_2, Gb_2,
     &                           PF_ALPHA, PF_SIGMA,
     &                           ParG, PerpG, G_2)
          DO i = 1, Nx
              XBP_curr(i) = 0.d0
              YBP_curr(i) = 0.d0
          END DO
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
     &             G_1, G_2, IH_1, IH_2, STPMX, hupd_method)

C     Optional DIIS acceleration (XMECP-style GDIIS/GEDIIS)
      stepnum = Nstep + 1
      if (trim(opt_method) .eq. 'gdiis' .or.
     &    trim(opt_method) .eq. 'gediis') then
          EOBJ = Ea_2
          CALL DIIS_UpdateX(Nx, stepnum, NDIIS, opt_method,
     &                      X_2, G_2, IH_2, EOBJ,
     &                      TGRMS, STPMX, X_3)
      endif

C     Test convergence
C       - Standard algorithms: use Harvey-style thresholds (TDE/TG/TDX)
C       - pf/lagrange: always use PF-specific criteria (XMECP-style)
      if (trim(algorithm) .eq. 'lagrange' .or.
     &    trim(algorithm) .eq. 'pf') then
          CALL TestConvergencePF(Nx, natom_actual, Nstep, AtNum,
     &                            Ea_1, Eb_1, Ea_2, Eb_2,
     &                            X_2, X_3, G_2,
     &                            PF_ALPHA, PF_SIGMA,
     &                            PF_TSTEP, PF_TGRAD,
     &                            Conv)
      else
          CALL TestConvergence(Nx, natom_actual, Nstep, AtNum,
     &                         Ea_2, Eb_2, X_2, X_3, ParG, PerpG,
     &                         G_2, Conv,
     &                         TDE, TDXMax, TDXRMS, TGMax, TGRMS)
      endif
     
      IF (Conv .eq. 1) THEN
          write(*,*) 'MECP optimization CONVERGED!'
          write(*,'(A,F18.10)') 'Final energy difference: ',
     &                          ABS(Ea_2 - Eb_2)
C         Write convergence status and criteria
          CALL WriteConvergenceInfo('CONVERGED', Nx, natom_actual,
     &                              Nstep, algorithm,
     &                              Ea_1, Eb_1, Ea_2, Eb_2,
     &                              X_2, X_3, ParG, PerpG, G_2,
     &                              TDE, TDXMax, TDXRMS, TGMax, TGRMS,
     &                              PF_ALPHA, PF_SIGMA,
     &                              PF_TSTEP, PF_TGRAD)
C         Write final geometry
          CALL WriteXYZ('final.xyz', natom_actual, AtNum, X_3,
     &                  'MECP Converged')
C         Clean up state file
          OPEN(UNIT=99, FILE='MECP.state')
          CLOSE(99, STATUS='DELETE')
C         Clean up DIIS history file
          inquire(file='MECP.diis', exist=file_exists)
          if (file_exists) then
              OPEN(UNIT=97, FILE='MECP.diis')
              CLOSE(97, STATUS='DELETE')
          endif
      ELSE
C         Write convergence status and criteria
          CALL WriteConvergenceInfo('NOT_CONVERGED', Nx, natom_actual,
     &                              Nstep, algorithm,
     &                              Ea_1, Eb_1, Ea_2, Eb_2,
     &                              X_2, X_3, ParG, PerpG, G_2,
     &                              TDE, TDXMax, TDXRMS, TGMax, TGRMS,
     &                              PF_ALPHA, PF_SIGMA,
     &                              PF_TSTEP, PF_TGRAD)
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

106   write(*,*) 'Error: Invalid value for --pf_alpha: ', trim(arg)
      stop
107   write(*,*) 'Error: Invalid value for --pf_sigma: ', trim(arg)
      stop

109   write(*,*) 'Error: Invalid value for --pf_tstep: ', trim(arg)
      stop
110   write(*,*) 'Error: Invalid value for --pf_tgrad: ', trim(arg)
      stop

108   write(*,*) 'Error: Invalid value for --ndiis: ', trim(arg)
      stop

999   END


C=====================================================================
      SUBROUTINE WriteConvergenceInfo(status, Nx, Natom_actual,
     &                                 Nstep, algorithm,
     &                                 Ea_prev, Eb_prev, Ea, Eb,
     &                                 X_2, X_3,
     &                                 ParG, PerpG, G,
     &                                 TDE, TDXMax, TDXRMS, TGMax, TGRMS,
     &                                 PF_ALPHA, PF_SIGMA,
     &                                 PF_TSTEP, PF_TGRAD)
      implicit none
      
      CHARACTER*(*) status
      CHARACTER*(*) algorithm
      INTEGER Nx, Natom_actual, Nstep, i
      DOUBLE PRECISION Ea_prev, Eb_prev
      DOUBLE PRECISION Ea, Eb, X_2(Nx), X_3(Nx), ParG(Nx)
      DOUBLE PRECISION PerpG(Nx), G(Nx)
      DOUBLE PRECISION TDE, TDXMax, TDXRMS, TGMax, TGRMS
      DOUBLE PRECISION PF_ALPHA, PF_SIGMA, PF_TSTEP, PF_TGRAD
      
      CHARACTER*3 flags(5)
      LOGICAL PConv(5)
      DOUBLE PRECISION DeltaX(Nx), DE, DXMax, DXRMS, GMax, GRMS
      DOUBLE PRECISION PpGRMS, PGRMS

C     PF-specific convergence data (XMECP-style)
      DOUBLE PRECISION de_pf, a_pf, s_pf, W_pf
      DOUBLE PRECISION PF_PREV, PF_CURR
      DOUBLE PRECISION func1, func2, func3
      DOUBLE PRECISION gnorm, u(Nx)
      CHARACTER*3 pflg(3)
      LOGICAL PFConv(3)
      DOUBLE PRECISION eps
      PARAMETER (eps=1.0d-12)

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

C     Append PF-specific convergence info for pf/lagrange.
      if (trim(algorithm) .eq. 'lagrange' .or.
     &    trim(algorithm) .eq. 'pf') then
          a_pf = PF_ALPHA
          s_pf = PF_SIGMA
          if (a_pf .le. eps) a_pf = 1.0d-2

          de_pf = Eb - Ea
          W_pf = 1.0d0 + (de_pf*de_pf)/(a_pf*a_pf)
          PF_CURR = 0.5d0*(Ea+Eb) + s_pf*a_pf*a_pf*log(W_pf)

          de_pf = Eb_prev - Ea_prev
          W_pf = 1.0d0 + (de_pf*de_pf)/(a_pf*a_pf)
          PF_PREV = 0.5d0*(Ea_prev+Eb_prev) + s_pf*a_pf*a_pf*log(W_pf)

          func1 = PF_PREV - PF_CURR

          gnorm = 0.d0
          do i = 1, Nx
              gnorm = gnorm + G(i)*G(i)
          end do
          gnorm = sqrt(gnorm)
          if (gnorm .gt. eps) then
              do i = 1, Nx
                  u(i) = G(i) / gnorm
              end do
          else
              do i = 1, Nx
                  u(i) = 0.d0
              end do
          endif

          func2 = 0.d0
          do i = 1, Nx
              func2 = func2 - G(i)*u(i)
          end do

          func3 = 0.d0
          do i = 1, Nx
              func3 = func3 + ( (-G(i) - func2*u(i))**2 )
          end do
          func3 = sqrt(func3)

          do i = 1, 3
              pflg(i) = " NO"
              PFConv(i) = .false.
          end do
          if (abs(func1) .lt. PF_TSTEP) then
              PFConv(1) = .true.
              pflg(1) = "YES"
          endif
          if (abs(func2) .lt. PF_TGRAD) then
              PFConv(2) = .true.
              pflg(2) = "YES"
          endif
          if (abs(func3) .lt. PF_TGRAD) then
              PFConv(3) = .true.
              pflg(3) = "YES"
          endif

          WRITE(98,'(A)') 'PF_Convergence_Criteria:'
          WRITE(98,'(A,F20.10)') 'PF_Objective: ', PF_CURR
          WRITE(98,'(A,F12.6,A,F9.6,A,A3)')
     &        'PF_Function1: ', func1, ' (', PF_TSTEP, ') ', pflg(1)
          WRITE(98,'(A,F12.6,A,F9.6,A,A3)')
     &        'PF_Function2: ', func2, ' (', PF_TGRAD, ') ', pflg(2)
          WRITE(98,'(A,F12.6,A,F9.6,A,A3)')
     &        'PF_Function3: ', func3, ' (', PF_TGRAD, ') ', pflg(3)
      endif
      CLOSE(98)
      
      RETURN
      END


C=====================================================================
C                    Lagrangian Constraint Gradient
C   (XMECP-style penalty-function / augmented-Lagrangian form)
C
C   Objective function (not explicitly used here):
C     G = (E1+E2)/2 + sigma*alpha^2*log(1 + (dE^2)/(alpha^2))
C   Its gradient:
C     grad = 0.5*(g1+g2) + (2*sigma*dE/(1+dE^2/alpha^2))*(g2-g1)
C
C   Here we return:
C     ParG  = mean gradient (0.5*(g1+g2))
C     PerpG = difference gradient (g2-g1)
C     G     = effective gradient for BFGS
C
C   Parameters:
C     PF_ALPHA (alpha) in Hartree
C     PF_SIGMA (sigma) in 1/Hartree
C=====================================================================
      SUBROUTINE LAGRANGE_Gradient(N,Ea,Eb,Ga,Gb,PF_ALPHA,PF_SIGMA,
     &                            ParG,PerpG,G)
      implicit none
      INTEGER N, i
      DOUBLE PRECISION Ea, Eb, PF_ALPHA, PF_SIGMA
      DOUBLE PRECISION Ga(N), Gb(N), ParG(N), PerpG(N), G(N)
      DOUBLE PRECISION de, a, s, W, gscale
      DOUBLE PRECISION eps
      PARAMETER (eps=1.0d-12)

      de = Eb - Ea
      a = PF_ALPHA
      s = PF_SIGMA

      IF (a .le. eps) THEN
          a = 1.0d-2
      END IF

      W = 1.0d0 + (de*de)/(a*a)
      gscale = 2.0d0 * s * de / W

      DO i = 1, N
          ParG(i) = 0.5d0 * (Ga(i) + Gb(i))
          PerpG(i) = Gb(i) - Ga(i)
          G(i) = ParG(i) + gscale * PerpG(i)
      END DO

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
C                          BPUPD_Gradient (修正版)
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
     &                   HI_1,HI_2,STPMX,HUPDMETHOD)
       implicit none

       ! Quasi-Newton step with selectable inverse-Hessian update:
       !   bfgs   - inverse BFGS update (default / legacy)
       !   psb    - inverse Powell-Symmetric-Broyden update (XMECP-style)
       !   sr1    - inverse Symmetric Rank-1 update
       !   bofill - Bofill mix of SR1 and PSB (inverse form)

       integer i, j, n, Nstep
       logical FirstStep
       double precision X_1(N), X_2(N), G_1(N), G_2(N)
       double precision HI_1(N,N), X_3(N), HI_2(N,N)
       double precision STPMX
       CHARACTER*(*) HUPDMETHOD

       double precision stpmax, DelG(N), HDelG(N), ChgeX(N), DelX(N)
       double precision w(N), fac, fad, fae, stpl, lgstst
       double precision Rvec(N)
       double precision sumyy, sumrr, yTr, denom, eps
       double precision HIPSB(N,N), HISR1(N,N), phi, tmp

       eps = 1.0d-12
       stpmax = STPMX * N

       IF (FirstStep) THEN
           DO i = 1, n
               ChgeX(i) = -.7d0 * G_2(i)
               DO j = 1, n
                   HI_2(i,j) = HI_1(i,j)
               end do
           end do
       ELSE
C         Build secant vectors: s = X_2 - X_1, y = G_2 - G_1
           DO i = 1, n
               DelG(i) = G_2(i) - G_1(i)
               DelX(i) = X_2(i) - X_1(i)
           end do
C         HDelG = HI_1 * y
           do i = 1, n
               HDelG(i) = 0.d0
               do j = 1, n
                   HDelG(i) = HDelG(i) + HI_1(i,j) * DelG(j)
               end do
           end do
C         Default: keep previous inverse Hessian
           do i = 1, n
               do j = 1, n
                   HI_2(i,j) = HI_1(i,j)
               end do
           end do

           if (trim(HUPDMETHOD) .eq. 'bfgs') then
C             --- inverse BFGS update (legacy) ---
               denom = 0.d0
               fae   = 0.d0
               do i = 1, n
                   denom = denom + DelG(i) * DelX(i)   ! y^T s
                   fae   = fae   + DelG(i) * HDelG(i)  ! y^T H y
               end do
               if (abs(denom) .gt. eps .and. abs(fae) .gt. eps) then
                   fac = 1.d0 / denom
                   fad = 1.d0 / fae
                   do i = 1, n
                       w(i) = fac * DelX(i) - fad * HDelG(i)
                   end do
                   do i = 1, n
                       do j = 1, n
                           HI_2(i,j) = HI_1(i,j)
     &                       + fac * DelX(i) * DelX(j)
     &                       - fad * HDelG(i) * HDelG(j)
     &                       + fae * w(i) * w(j)
                       end do
                   end do
               endif

           else
C             Residual r = s - H y  (for inverse-Hessian updates)
               do i = 1, n
                   Rvec(i) = DelX(i) - HDelG(i)
               end do
               sumyy = 0.d0
               sumrr = 0.d0
               yTr   = 0.d0
               do i = 1, n
                   sumyy = sumyy + DelG(i) * DelG(i)
                   sumrr = sumrr + Rvec(i) * Rvec(i)
                   yTr   = yTr   + DelG(i) * Rvec(i)
               end do

               if (trim(HUPDMETHOD) .eq. 'psb') then
C                 --- inverse PSB update ---
                   if (sumyy .gt. eps) then
                       do i = 1, n
                           do j = 1, n
                               HI_2(i,j) = HI_1(i,j)
     &                         + (Rvec(i) * DelG(j)
     &                         +  DelG(i) * Rvec(j)) / sumyy
     &                         - (yTr / (sumyy*sumyy))
     &                         *  DelG(i) * DelG(j)
                           end do
                       end do
                   endif

               else if (trim(HUPDMETHOD) .eq. 'sr1') then
C                 --- inverse SR1 update ---
                   denom = yTr   ! r^T y
                   if (abs(denom) .gt. eps .and. sumrr .gt. eps) then
                       do i = 1, n
                           do j = 1, n
                               HI_2(i,j) = HI_1(i,j)
     &                         + Rvec(i) * Rvec(j) / denom
                           end do
                       end do
                   endif

               else if (trim(HUPDMETHOD) .eq. 'bofill') then
C                 --- Bofill mix of SR1 and PSB (inverse form) ---
C                 Build PSB candidate
                   do i = 1, n
                       do j = 1, n
                           HIPSB(i,j) = HI_1(i,j)
                           HISR1(i,j) = HI_1(i,j)
                       end do
                   end do
                   if (sumyy .gt. eps) then
                       do i = 1, n
                           do j = 1, n
                               HIPSB(i,j) = HI_1(i,j)
     &                         + (Rvec(i) * DelG(j)
     &                         +  DelG(i) * Rvec(j)) / sumyy
     &                         - (yTr / (sumyy*sumyy))
     &                         *  DelG(i) * DelG(j)
                           end do
                       end do
                   endif
C                 Build SR1 candidate
                   denom = yTr
                   if (abs(denom) .gt. eps .and. sumrr .gt. eps) then
                       do i = 1, n
                           do j = 1, n
                               HISR1(i,j) = HI_1(i,j)
     &                         + Rvec(i) * Rvec(j) / denom
                           end do
                       end do
                   endif
C                 Mixing parameter phi (0..1)
                   phi = 0.d0
                   if (sumrr .gt. eps .and. sumyy .gt. eps) then
                       phi = (yTr*yTr) / (sumrr*sumyy)
                       if (phi .lt. 0.d0) phi = 0.d0
                       if (phi .gt. 1.d0) phi = 1.d0
                   endif
                   do i = 1, n
                       do j = 1, n
                           HI_2(i,j) = phi * HISR1(i,j)
     &                               + (1.d0-phi) * HIPSB(i,j)
                       end do
                   end do
               endif
           endif

C         Symmetrize to suppress numerical drift
           do i = 1, n
               do j = i + 1, n
                   tmp = 0.5d0 * (HI_2(i,j) + HI_2(j,i))
                   HI_2(i,j) = tmp
                   HI_2(j,i) = tmp
               end do
           end do

C         Compute quasi-Newton step: dx = -HI_2 * g
           do i = 1, n
               ChgeX(i) = 0.d0
               do j = 1, n
                   ChgeX(i) = ChgeX(i) - HI_2(i,j) * G_2(j)
               end do
           end do
       END IF

C     Step-length control (trust radius + max component)
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
C                    PF-specific Convergence Test
C   For pf/lagrange algorithms: always use PF-exclusive criteria.
C   Reference: XMECP optimizer.isconver_pf()
C
C   Criteria:
C     func1 = PF_prev - PF_curr  (objective change)
C     func2 = dot(-grad, u)      (|grad|, u is normalized grad)
C     func3 = norm(-grad - func2*u) (orthogonal component)
C
C   Converged if:
C     |func1| < PF_TSTEP and |func2| < PF_TGRAD and |func3| < PF_TGRAD
C=====================================================================
      SUBROUTINE TestConvergencePF(Nx,Natom_actual,Nstep,AtNum,
     &                             Ea_prev,Eb_prev,Ea,Eb,
     &                             X_2,X_3,G,
     &                             PF_ALPHA,PF_SIGMA,
     &                             PF_TSTEP,PF_TGRAD,
     &                             Conv)
      implicit none

      INTEGER Nx, Natom_actual, AtNum(*), Nstep, Conv
      INTEGER i, k
      DOUBLE PRECISION Ea_prev, Eb_prev, Ea, Eb
      DOUBLE PRECISION X_2(Nx), X_3(Nx), G(Nx)
      DOUBLE PRECISION PF_ALPHA, PF_SIGMA, PF_TSTEP, PF_TGRAD

      DOUBLE PRECISION a, s, de, W
      DOUBLE PRECISION PF_PREV, PF_CURR
      DOUBLE PRECISION func1, func2, func3
      DOUBLE PRECISION gnorm, u(Nx)
      CHARACTER*3 flags(3)
      LOGICAL PFConv(3)
      CHARACTER*2 element
      DOUBLE PRECISION eps
      PARAMETER (eps=1.0d-12)

      a = PF_ALPHA
      s = PF_SIGMA
      if (a .le. eps) a = 1.0d-2

C     Compute PF objective at current and previous steps
      de = Eb - Ea
      W = 1.0d0 + (de*de)/(a*a)
      PF_CURR = 0.5d0*(Ea+Eb) + s*a*a*log(W)

      de = Eb_prev - Ea_prev
      W = 1.0d0 + (de*de)/(a*a)
      PF_PREV = 0.5d0*(Ea_prev+Eb_prev) + s*a*a*log(W)

      func1 = PF_PREV - PF_CURR

C     Compute func2/func3 from effective gradient
      gnorm = 0.d0
      do i = 1, Nx
          gnorm = gnorm + G(i)*G(i)
      end do
      gnorm = sqrt(gnorm)

      if (gnorm .gt. eps) then
          do i = 1, Nx
              u(i) = G(i) / gnorm
          end do
      else
          do i = 1, Nx
              u(i) = 0.d0
          end do
      endif

      func2 = 0.d0
      do i = 1, Nx
          func2 = func2 - G(i)*u(i)
      end do

      func3 = 0.d0
      do i = 1, Nx
          func3 = func3 + ( (-G(i) - func2*u(i))**2 )
      end do
      func3 = sqrt(func3)

C     Evaluate PF convergence flags
      Conv = 0
      do i = 1, 3
          flags(i) = " NO"
          PFConv(i) = .false.
      end do

      if (abs(func1) .lt. PF_TSTEP) then
          PFConv(1) = .true.
          flags(1) = "YES"
      endif
      if (abs(func2) .lt. PF_TGRAD) then
          PFConv(2) = .true.
          flags(2) = "YES"
      endif
      if (abs(func3) .lt. PF_TGRAD) then
          PFConv(3) = .true.
          flags(3) = "YES"
      endif

      Nstep = Nstep + 1
      if (PFConv(1) .and. PFConv(2) .and. PFConv(3)) then
          Conv = 1
      else
          Conv = 0
      endif

C     Write header for first step
      IF (Nstep .eq. 1) THEN
          write(*,*)
          write(*,'(A)') '  =========================================='
          write(*,'(A)') '       Geometry Optimization of an MECP'
          write(*,'(A)') '       Penalty Function / Augmented-Lagrangian'
          write(*,'(A)') '       baneMECP v1.1, Sept. 2025'
          write(*,'(A)') '       modified from sobMECP@sobereva'
          write(*,'(A)') '       feedback: banerxmd@gmail.com'
          write(*,'(A)') '  =========================================='
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

C     Print PF convergence info to screen (XMECP-style)
      write(*,*)
      write(*,'(A,I4,A)') '  ===== MECP Optimization Step',
     &                     Nstep, ' (PF/Lagrange) ====='
      write(*,'(A,F20.10)') '    Energy of State 1:  ', Ea
      write(*,'(A,F20.10)') '    Energy of State 2:  ', Eb
      write(*,'(A,F20.10)') '    Energy Difference:  ', ABS(Ea-Eb)
      write(*,'(A,F20.10)') '    PF Objective:       ', PF_CURR
      write(*,*)
      write(*,'(A)') '  PF Convergence Check:'
      write(*,'(A,F12.6,A,F9.6,A,A3)')
     &    '    PF Function1:     ', func1, ' (', PF_TSTEP, ')  ', flags(1)
      write(*,'(A,F12.6,A,F9.6,A,A3)')
     &    '    PF Function2:     ', func2, ' (', PF_TGRAD, ')  ', flags(2)
      write(*,'(A,F12.6,A,F9.6,A,A3)')
     &    '    PF Function3:     ', func3, ' (', PF_TGRAD, ')  ', flags(3)
      write(*,*)

      return
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

C     Write header for first step
      IF (Nstep .eq. 1) THEN
          write(*,*)
          write(*,'(A)') '  =========================================='
          write(*,'(A)') '       Geometry Optimization of an MECP'
          write(*,'(A)') '       Original: J. N. Harvey, March 1999'
          write(*,'(A)') '       baneMECP v1.1, Sept. 2025'
          write(*,'(A)') '       modified from sobMECP@sobereva'
          write(*,'(A)') '       feedback: banerxmd@gmail.com'
          write(*,'(A)') '  =========================================='
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

C=====================================================================
C                    DIIS acceleration (GDIIS / GEDIIS)
C=====================================================================
      SUBROUTINE DIIS_UpdateX(Nx, Step, NDIIS, OPTMETHOD,
     &                        XCURR, GCURR, HICURR, EOBJ,
     &                        TGRMS, STPMX, XNEW)
      implicit none

      INTEGER Nx, Step, NDIIS
      CHARACTER*(*) OPTMETHOD
      DOUBLE PRECISION XCURR(Nx), GCURR(Nx), HICURR(Nx,Nx)
      DOUBLE PRECISION EOBJ, TGRMS, STPMX
      DOUBLE PRECISION XNEW(Nx)

      INTEGER MAXDIIS
      parameter (MAXDIIS=8)

      INTEGER nHist
      DOUBLE PRECISION XH(MAXDIIS, Nx), GH(MAXDIIS, Nx)
      DOUBLE PRECISION EH(MAXDIIS)
      DOUBLE PRECISION XGDIIS(Nx), XGEDIIS(Nx)
      INTEGER info
      LOGICAL do_diis

C     Read and update history (always)
      call DIIS_ReadHistory('MECP.diis', MAXDIIS, Nx, nHist, XH, GH, EH)
      call DIIS_AppendHistory(MAXDIIS, Nx, nHist, NDIIS,
     &                        XCURR, GCURR, EOBJ, XH, GH, EH)
      call DIIS_WriteHistory('MECP.diis', Nx, nHist, XH, GH, EH)

C     Only start DIIS after a couple of steps (match XMECP behavior)
      do_diis = .false.
      if (Step .ge. 3) do_diis = .true.
      if (NDIIS .le. 0) do_diis = .false.
      if (trim(OPTMETHOD) .eq. 'bfgs') do_diis = .false.
      if (.not. do_diis) return

C     --- GDIIS ---
      call DIIS_GDIIS(Nx, nHist, XH, GH, HICURR, XGDIIS, info)
      if (info .ne. 0) then
          write(*,*) 'Warning: GDIIS failed (singular system).'
          return
      endif

      if (trim(OPTMETHOD) .eq. 'gediis') then
C         --- GEDIIS (optional) ---
          call DIIS_GEDIIS(Nx, nHist, XH, GH, EH, HICURR, XGEDIIS, info)
          if (info .ne. 0) then
              write(*,*) 'Warning: GEDIIS failed; using GDIIS only.'
              call DIIS_CopyVec(Nx, XGDIIS, XNEW)
          else
C             Blend GDIIS and GEDIIS proposals for stability
              call DIIS_BlendVec(Nx, XGDIIS, XGEDIIS, 0.5d0, XNEW)
          endif
          write(*,'(A,I4)') 'Entering GEDIIS Step ', Step
      else
          call DIIS_CopyVec(Nx, XGDIIS, XNEW)
          write(*,'(A,I4)') 'Entering GDIIS Step ', Step
      endif

C     Apply step-size control (trust radius + mild damping near convergence)
      call DIIS_StepControl(Nx, XCURR, GCURR, TGRMS, STPMX, XNEW)

      RETURN
      END


C=====================================================================
      SUBROUTINE DIIS_ReadHistory(fname, MAXDIIS, Nx, nHist, XH, GH, EH)
      implicit none

      CHARACTER*(*) fname
      INTEGER MAXDIIS, Nx, nHist
      DOUBLE PRECISION XH(MAXDIIS, Nx), GH(MAXDIIS, Nx), EH(MAXDIIS)

      LOGICAL exists
      INTEGER nfile, nxfile, i, j, idx, start
      DOUBLE PRECISION etmp
      DOUBLE PRECISION xtmp(Nx), gtmp(Nx)

      inquire(file=fname, exist=exists)
      if (.not. exists) then
          nHist = 0
          return
      endif

      OPEN(UNIT=91, FILE=fname, STATUS='OLD', ERR=900)
      READ(91,*,END=900,ERR=900) nfile
      READ(91,*,END=900,ERR=900) nxfile

      if (nxfile .ne. Nx) then
          nHist = 0
          CLOSE(91)
          return
      endif

      if (nfile .le. 0) then
          nHist = 0
          CLOSE(91)
          return
      endif

      if (nfile .gt. MAXDIIS) then
          start = nfile - MAXDIIS + 1
          nHist = MAXDIIS
      else
          start = 1
          nHist = nfile
      endif

      idx = 0
      do i = 1, nfile
          READ(91,*,END=900,ERR=900) etmp
          do j = 1, Nx
              READ(91,*,END=900,ERR=900) xtmp(j)
          end do
          do j = 1, Nx
              READ(91,*,END=900,ERR=900) gtmp(j)
          end do

          if (i .ge. start) then
              idx = idx + 1
              EH(idx) = etmp
              do j = 1, Nx
                  XH(idx, j) = xtmp(j)
                  GH(idx, j) = gtmp(j)
              end do
          endif
      end do

      CLOSE(91)
      return

  900 continue
      nHist = 0
      CLOSE(91)
      return
      END


C=====================================================================
      SUBROUTINE DIIS_WriteHistory(fname, Nx, nHist, XH, GH, EH)
      implicit none

      CHARACTER*(*) fname
      INTEGER Nx, nHist, i, j
      DOUBLE PRECISION XH(8, Nx), GH(8, Nx), EH(8)

      OPEN(UNIT=92, FILE=fname, STATUS='UNKNOWN')
      WRITE(92,*) nHist
      WRITE(92,*) Nx
      do i = 1, nHist
          WRITE(92,'(F20.12)') EH(i)
          do j = 1, Nx
              WRITE(92,'(F20.12)') XH(i, j)
          end do
          do j = 1, Nx
              WRITE(92,'(F20.12)') GH(i, j)
          end do
      end do
      CLOSE(92)
      return
      END


C=====================================================================
      SUBROUTINE DIIS_AppendHistory(MAXDIIS, Nx, nHist, NDIIS,
     &                              XCURR, GCURR, EOBJ, XH, GH, EH)
      implicit none

      INTEGER MAXDIIS, Nx, nHist, NDIIS
      DOUBLE PRECISION XCURR(Nx), GCURR(Nx), EOBJ
      DOUBLE PRECISION XH(MAXDIIS, Nx), GH(MAXDIIS, Nx), EH(MAXDIIS)

      INTEGER keep, i, j

      keep = NDIIS
      if (keep .le. 0) keep = 1
      if (keep .gt. MAXDIIS) keep = MAXDIIS

      if (nHist .ge. keep) then
C         Drop oldest to keep history length <= keep-1 before append
          do i = 1, keep-1
              EH(i) = EH(i+1)
              do j = 1, Nx
                  XH(i, j) = XH(i+1, j)
                  GH(i, j) = GH(i+1, j)
              end do
          end do
          nHist = keep - 1
      endif

      nHist = nHist + 1
      EH(nHist) = EOBJ
      do j = 1, Nx
          XH(nHist, j) = XCURR(j)
          GH(nHist, j) = GCURR(j)
      end do

      return
      END


C=====================================================================
      SUBROUTINE DIIS_GDIIS(Nx, nHist, XH, GH, HI, XOUT, info)
      implicit none

      INTEGER Nx, nHist, info
      DOUBLE PRECISION XH(8, Nx), GH(8, Nx), HI(Nx, Nx)
      DOUBLE PRECISION XOUT(Nx)

      INTEGER i, j, k, m, nsys
      DOUBLE PRECISION EV(8, Nx)
      DOUBLE PRECISION B(8, 8)
      DOUBLE PRECISION A(10, 10), rhs(10), sol(10)
      DOUBLE PRECISION tmpX(Nx), tmpG(Nx), step(Nx)

      info = 0
      if (nHist .le. 0) then
          info = 1
          return
      endif

C     Error vectors: e_i = HI * g_i
      do i = 1, nHist
          do k = 1, Nx
              EV(i, k) = 0.d0
              do j = 1, Nx
                  EV(i, k) = EV(i, k) + HI(k, j) * GH(i, j)
              end do
          end do
      end do

C     B matrix = e_i · e_j
      do i = 1, nHist
          do j = 1, nHist
              B(i, j) = 0.d0
              do k = 1, Nx
                  B(i, j) = B(i, j) + EV(i, k) * EV(j, k)
              end do
          end do
      end do

C     Build augmented system
      nsys = nHist + 1
      do i = 1, nsys
          rhs(i) = 0.d0
          do j = 1, nsys
              A(i, j) = 0.d0
          end do
      end do

      do i = 1, nHist
          do j = 1, nHist
              A(i, j) = B(i, j)
          end do
          A(i, nsys) = 1.d0
          A(nsys, i) = 1.d0
      end do
      A(nsys, nsys) = 0.d0
      rhs(nsys) = 1.d0

      call DIIS_Solve(nsys, A, rhs, sol, info)
      if (info .ne. 0) return

C     tmpX = Σ c_i X_i ; tmpG = Σ c_i G_i
      do k = 1, Nx
          tmpX(k) = 0.d0
          tmpG(k) = 0.d0
      end do

      do i = 1, nHist
          do k = 1, Nx
              tmpX(k) = tmpX(k) + sol(i) * XH(i, k)
              tmpG(k) = tmpG(k) + sol(i) * GH(i, k)
          end do
      end do

C     Xout = tmpX - HI * tmpG
      do k = 1, Nx
          step(k) = 0.d0
          do j = 1, Nx
              step(k) = step(k) + HI(k, j) * tmpG(j)
          end do
          XOUT(k) = tmpX(k) - step(k)
      end do

      return
      END


C=====================================================================
      SUBROUTINE DIIS_GEDIIS(Nx, nHist, XH, GH, EH, HI, XOUT, info)
      implicit none

      INTEGER Nx, nHist, info
      DOUBLE PRECISION XH(8, Nx), GH(8, Nx), EH(8), HI(Nx, Nx)
      DOUBLE PRECISION XOUT(Nx)

      INTEGER i, j, k, nsys
      DOUBLE PRECISION M(8, 8)
      DOUBLE PRECISION A(10, 10), rhs(10), sol(10)
      DOUBLE PRECISION tmpX(Nx), tmpG(Nx), step(Nx)
      DOUBLE PRECISION val

      info = 0
      if (nHist .le. 0) then
          info = 1
          return
      endif

C     Build GEDIIS matrix M(i,j) = - (g_i - g_j) · (x_i - x_j)
      do i = 1, nHist
          do j = 1, nHist
              M(i, j) = 0.d0
          end do
      end do

      do i = 1, nHist
          M(i, i) = 0.d0
          do j = i+1, nHist
              val = 0.d0
              do k = 1, Nx
                  val = val + (GH(i, k) - GH(j, k)) * (XH(i, k) - XH(j, k))
              end do
              M(i, j) = -val
              M(j, i) = M(i, j)
          end do
      end do

C     Build augmented system
      nsys = nHist + 1
      do i = 1, nsys
          rhs(i) = 0.d0
          do j = 1, nsys
              A(i, j) = 0.d0
          end do
      end do

      do i = 1, nHist
          do j = 1, nHist
              A(i, j) = M(i, j)
          end do
          A(i, nsys) = 1.d0
          A(nsys, i) = 1.d0
          rhs(i) = -EH(i)
      end do
      A(nsys, nsys) = 0.d0
      rhs(nsys) = 1.d0

      call DIIS_Solve(nsys, A, rhs, sol, info)
      if (info .ne. 0) return

C     tmpX = Σ c_i X_i ; tmpG = Σ c_i G_i
      do k = 1, Nx
          tmpX(k) = 0.d0
          tmpG(k) = 0.d0
      end do

      do i = 1, nHist
          do k = 1, Nx
              tmpX(k) = tmpX(k) + sol(i) * XH(i, k)
              tmpG(k) = tmpG(k) + sol(i) * GH(i, k)
          end do
      end do

C     Xout = tmpX - HI * tmpG
      do k = 1, Nx
          step(k) = 0.d0
          do j = 1, Nx
              step(k) = step(k) + HI(k, j) * tmpG(j)
          end do
          XOUT(k) = tmpX(k) - step(k)
      end do

      return
      END


C=====================================================================
      SUBROUTINE DIIS_StepControl(Nx, XCURR, GCURR, TGRMS, STPMX, XNEW)
      implicit none

      INTEGER Nx, i
      DOUBLE PRECISION XCURR(Nx), GCURR(Nx), TGRMS, STPMX, XNEW(Nx)
      DOUBLE PRECISION dX(Nx), g2, grms, factor
      DOUBLE PRECISION stpmax, stpl, lgstst

C     Mild damping near convergence (similar spirit to XMECP)
      g2 = 0.d0
      do i = 1, Nx
          g2 = g2 + GCURR(i) * GCURR(i)
      end do
      grms = SQRT(g2 / Nx)

      factor = 1.d0
      if (grms .lt. 10.d0 * TGRMS) factor = 0.5d0

      do i = 1, Nx
          dX(i) = (XNEW(i) - XCURR(i)) * factor
      end do

C     Apply trust radius limits (same logic as UpdateX)
      stpmax = STPMX * Nx
      stpl = 0.d0
      do i = 1, Nx
          stpl = stpl + dX(i) * dX(i)
      end do
      stpl = SQRT(stpl)

      if (stpl .gt. stpmax .and. stpl .gt. 0.d0) then
          do i = 1, Nx
              dX(i) = dX(i) / stpl * stpmax
          end do
      endif

      lgstst = 0.d0
      do i = 1, Nx
          if (ABS(dX(i)) .gt. lgstst) lgstst = ABS(dX(i))
      end do

      if (lgstst .gt. STPMX .and. lgstst .gt. 0.d0) then
          do i = 1, Nx
              dX(i) = dX(i) / lgstst * STPMX
          end do
      endif

      do i = 1, Nx
          XNEW(i) = XCURR(i) + dX(i)
      end do

      return
      END


C=====================================================================
      SUBROUTINE DIIS_Solve(N, A, b, x, info)
      implicit none

      INTEGER N, info
      DOUBLE PRECISION A(10,10), b(10), x(10)

      INTEGER i, j, k, p
      DOUBLE PRECISION tmp, factor, sum, maxval

      info = 0

      do k = 1, N-1
C         Pivot selection
          p = k
          maxval = ABS(A(k, k))
          do i = k+1, N
              if (ABS(A(i, k)) .gt. maxval) then
                  maxval = ABS(A(i, k))
                  p = i
              endif
          end do

          if (maxval .lt. 1.d-14) then
              info = 1
              return
          endif

C         Swap rows if needed
          if (p .ne. k) then
              do j = k, N
                  tmp = A(k, j)
                  A(k, j) = A(p, j)
                  A(p, j) = tmp
              end do
              tmp = b(k)
              b(k) = b(p)
              b(p) = tmp
          endif

C         Elimination
          do i = k+1, N
              factor = A(i, k) / A(k, k)
              A(i, k) = 0.d0
              do j = k+1, N
                  A(i, j) = A(i, j) - factor * A(k, j)
              end do
              b(i) = b(i) - factor * b(k)
          end do
      end do

      if (ABS(A(N, N)) .lt. 1.d-14) then
          info = 1
          return
      endif

C     Back substitution
      x(N) = b(N) / A(N, N)
      do i = N-1, 1, -1
          sum = b(i)
          do j = i+1, N
              sum = sum - A(i, j) * x(j)
          end do
          x(i) = sum / A(i, i)
      end do

      return
      END


C=====================================================================
      SUBROUTINE DIIS_CopyVec(Nx, src, dst)
      implicit none
      INTEGER Nx, i
      DOUBLE PRECISION src(Nx), dst(Nx)
      do i = 1, Nx
          dst(i) = src(i)
      end do
      return
      END


C=====================================================================
      SUBROUTINE DIIS_BlendVec(Nx, a, b, w, out)
      implicit none
      INTEGER Nx, i
      DOUBLE PRECISION a(Nx), b(Nx), out(Nx), w
      do i = 1, Nx
          out(i) = w * a(i) + (1.d0 - w) * b(i)
      end do
      return
      END
