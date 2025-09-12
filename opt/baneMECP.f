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
C       gfortran -O3 -march=native -ffast-math -ffixed-line-length-none sbaneMECP.f -o MECP.x
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
      
      CHARACTER*256 inputname, filename
      LOGICAL file_exists, first_step
      INTEGER nargs, iargc, natom_actual
      
C     Get input filename from command line
      nargs = iargc()
      if (nargs .lt. 1) then
          write(*,*) 'Usage: ./MECP.x inputname'
          write(*,*) 'Input files required:'
          write(*,*) '  inputname.xyz   - geometry file'
          write(*,*) '  inputname.grad1 - gradient file for state 1'
          write(*,*) '  inputname.grad2 - gradient file for state 2'
          stop
      endif
      
      call getarg(1, inputname)
      
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
      ELSE
C         Read previous state
          CALL ReadState(Natom, Nx, natom_actual, Nstep,
     &                   X_1, IH_1, Ea_1, Eb_1, Ga_1, Gb_1, G_1)
      END IF
      
C     Compute the Effective Gradient (ORIGINAL ALGORITHM)
      CALL Effective_Gradient(Nx, Ea_2, Eb_2, Ga_2, Gb_2, 
     &                        ParG, PerpG, G_2)
     
C     Perform BFGS update (ORIGINAL ALGORITHM)
      CALL UpdateX(Nx, Nstep, first_step, X_1, X_2, X_3, 
     &             G_1, G_2, IH_1, IH_2)
     
C     Test convergence
      CALL TestConvergence(Nx, natom_actual, Nstep, AtNum, 
     &                     Ea_2, Eb_2, X_2, X_3, ParG, PerpG, 
     &                     G_2, Conv)
     
      IF (Conv .eq. 1) THEN
          write(*,*) 'MECP optimization CONVERGED!'
          write(*,'(A,F18.10)') 'Final energy difference: ',
     &                          ABS(Ea_2 - Eb_2)
C         Write convergence status and criteria
          CALL WriteConvergenceInfo('CONVERGED', Nx, natom_actual, 
     &                              Nstep, Ea_2, Eb_2, X_2, X_3, 
     &                              ParG, PerpG, G_2)
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
     &                              ParG, PerpG, G_2)
C         Write new geometry for next iteration
          CALL WriteXYZ('new.xyz', natom_actual, AtNum, X_3,
     &                  'MECP Step')
C         Save current state
          CALL WriteState(Natom, Nx, natom_actual, Nstep,
     &                    X_2, X_3, IH_2, Ea_2, Eb_2, 
     &                    Ga_2, Gb_2, G_2)
      END IF

      END


C=====================================================================
      SUBROUTINE WriteConvergenceInfo(status, Nx, Natom_actual, 
     &                                 Nstep, Ea, Eb, X_2, X_3, 
     &                                 ParG, PerpG, G)
      implicit none
      
      CHARACTER*(*) status
      INTEGER Nx, Natom_actual, Nstep, i
      DOUBLE PRECISION Ea, Eb, X_2(Nx), X_3(Nx), ParG(Nx)
      DOUBLE PRECISION PerpG(Nx), G(Nx)
      
      CHARACTER*3 flags(5)
      LOGICAL PConv(5)
      DOUBLE PRECISION DeltaX(Nx), DE, DXMax, DXRMS, GMax, GRMS
      DOUBLE PRECISION PpGRMS, PGRMS
      DOUBLE PRECISION TDE, TDXMax, TDXRMS, TGMax, TGRMS
      PARAMETER (TDE=5.d-5,TDXMax=4.d-3,TDXRMS=2.5d-3)
      PARAMETER (TGMax=7.d-4,TGRMS=5.d-4)

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
          write(*,*) 'File has', Natom_actual, 'atoms'
          write(*,*) 'Compiled for max', MaxNatom, 'atoms'
          write(*,*) 'Recompile with larger NUMATOM'
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
C         Convert from force to gradient (negate)
C         Assume input is in Hartree/Bohr, convert to Hartree/Ang
          Ga(k) = -Ga(k) / bohr
          Ga(k+1) = -Ga(k+1) / bohr
          Ga(k+2) = -Ga(k+2) / bohr
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
          Gb(k) = -Gb(k) / bohr
          Gb(k+1) = -Gb(k+1) / bohr
          Gb(k+2) = -Gb(k+2) / bohr
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
     &                     X_1, IH, Ea, Eb, Ga, Gb, G)
      implicit none
      
      INTEGER MaxNatom, MaxNx, Natom_actual, Nstep, i, j
      INTEGER nx_actual
      DOUBLE PRECISION X_1(MaxNx), IH(MaxNx,MaxNx)
      DOUBLE PRECISION Ea, Eb, Ga(MaxNx), Gb(MaxNx), G(MaxNx)
      
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
      
      CLOSE(30)
      RETURN
      END


C=====================================================================
      SUBROUTINE WriteState(MaxNatom, MaxNx, Natom_actual, Nstep,
     &                      X_2, X_3, IH, Ea, Eb, Ga, Gb, G)
      implicit none
      
      INTEGER MaxNatom, MaxNx, Natom_actual, Nstep, i, j
      INTEGER nx_actual
      DOUBLE PRECISION X_2(MaxNx), X_3(MaxNx), IH(MaxNx,MaxNx)
      DOUBLE PRECISION Ea, Eb, Ga(MaxNx), Gb(MaxNx), G(MaxNx)
      
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
      
      CLOSE(31)
      RETURN
      END


C=====================================================================
      SUBROUTINE ElementToNumber(element, number)
      implicit none
      
      CHARACTER*2 element
      INTEGER number
      CHARACTER*2 symbols(118)
      INTEGER i
      
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
      
      number = 0
      DO i = 1, 118
          IF (element .eq. symbols(i)) THEN
              number = i
              RETURN
          END IF
      END DO
      
C     If not found, might be lowercase or mixed case
      IF (number .eq. 0) THEN
          write(*,*) 'Warning: Unknown element ', element
          number = 6  ! Default to carbon
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
     &                   HI_1,HI_2)
       implicit none
       
       ! Specially Adapted BFGS routine from Numerical Recipes
       
       integer i, j, k, n, Nstep
       logical FirstStep
       double precision X_1(N), X_2(N), G_1(N), G_2(N)
       double precision HI_1(N,N), X_3(N), HI_2(N,N)
       
       double precision stpmax, DelG(N), HDelG(N), ChgeX(N), DelX(N), w(N), fac,
     & fad, fae, sumdg, sumdx, stpl, lgstst, STPMX
       parameter (STPMX = 0.1d0)
       
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
     &                           X_2,X_3,ParG,PerpG,G,Conv)
      implicit none

      INTEGER Nx, Natom_actual, AtNum(*), Nstep, Conv, i, j, k
      DOUBLE PRECISION Ea, Eb, X_2(Nx), ParG(Nx), PerpG(Nx)
      DOUBLE PRECISION X_3(Nx), G(Nx)

      CHARACTER*3 flags(5)
      LOGICAL PConv(5)
      DOUBLE PRECISION DeltaX(Nx), DE, DXMax, DXRMS, GMax, GRMS
      DOUBLE PRECISION PpGRMS, PGRMS, TDE
      DOUBLE PRECISION TDXMax, TDxRMS, TGMax, TGRMS
      PARAMETER (TDE=5.d-5,TDXMax=4.d-3,TDXRMS=2.5d-3)
      PARAMETER (TGMax=7.d-4,TGRMS=5.d-4)
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
          write(*,'(A)') '       baneMECP v1.0, Sept. 2025'
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