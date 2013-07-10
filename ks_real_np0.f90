MODULE nrtype
  INTEGER, PARAMETER :: DP = KIND(1.0D0)
  REAL(DP), PARAMETER :: PI=3.141592653589793238462643383279502884197_dp
  REAL(DP), PARAMETER :: PIO2=1.57079632679489661923132169163975144209858_dp
  REAL(DP), PARAMETER :: TWOPI=6.283185307179586476925286766559005768394_dp
END MODULE nrtype

Module fft_coeff
  IMPLICIT NONE
  COMPLEX(KIND=8), ALLOCATABLE, DIMENSION(:) :: coef_d1, coef_d2, coef_d3, coef_d4
  INTEGER(KIND=4) :: fft_coeff_flag =0, fft_plan_flag=0, nend

  INTEGER(KIND=8) :: p_fw0, p_bw1, p_bw2, p_bw3, p_bw4
  REAL(KIND=8), ALLOCATABLE, DIMENSION(:) :: hh, d1h, d2h, d3h, d4h
  COMPLEX(KIND=8), ALLOCATABLE, DIMENSION(:) :: h_fft, dummy_c

  INCLUDE 'fftw3.f90'

  CONTAINS

    SUBROUTINE calculate_modes_r2c(nn, ll)
      USE nrtype
      IMPLICIT NONE
      
      INTEGER(KIND=4), INTENT(IN) :: nn
      REAL(KIND=8), INTENT(IN) :: ll
      INTEGER(KIND=4) :: k, kend

      fft_coeff_flag = 1
!
! nend: total number of Fourier coefficients
!
      IF(MOD(nn,2).EQ.0)THEN
         kend = nn/2
         nend = nn/2+1
      ELSE
         kend = (nn+1)/2
         nend = (nn+1)/2
      END IF

      ALLOCATE(coef_d1(1:nend))
      ALLOCATE(coef_d2(1:nend))
      ALLOCATE(coef_d3(1:nend))
      ALLOCATE(coef_d4(1:nend))
!
!  Setup fourier coefficients of derivatives
!
      coef_d1(1)=(0.0d0, 0.0d0)
      coef_d2(1)=(0.0d0, 0.0d0)
      coef_d3(1)=(0.0d0, 0.0d0)
      coef_d4(1)=(0.0d0, 0.0d0)
      DO k = 2, kend
         coef_d1(k) = DCMPLX(0.0d0, TWOPI*DBLE(k-1)/ll)
         coef_d2(k) = DCMPLX(-(TWOPI*DBLE(k-1)/ll)**2, 0.0d0)
         coef_d3(k) = DCMPLX(0.0d0, -(TWOPI*DBLE(k-1)/ll)**3)
         coef_d4(k) = DCMPLX((TWOPI*DBLE(k-1)/ll)**4, 0.0d0)
      END DO
!
!  Coefficients for the Nyquist mode
!
      IF(MOD(nn,2).EQ.0)THEN
         coef_d1(nend) = (0.0d0, 0.0d0)
         coef_d2(nend) = DCMPLX(-(TWOPI*DBLE(nn/2)/ll)**2, 0.0d0)
         coef_d3(nend) = (0.0d0, 0.0d0)
         coef_d4(nend) = DCMPLX((TWOPI*DBLE(nn/2)/ll)**4, 0.0d0)
      END IF

      RETURN
    END SUBROUTINE calculate_modes_r2c

    SUBROUTINE create_fftw_plan(nn)
      IMPLICIT NONE
      INTEGER(KIND=4), INTENT(IN) :: nn
      INTEGER :: kend

      fft_plan_flag = 1

      IF(MOD(nn,2).EQ.0)THEN
         kend = nn/2+1
      ELSE
         kend = (nn+1)/2
      END IF

      ALLOCATE(hh(1:nn))
      ALLOCATE(d1h(1:nn))
      ALLOCATE(d2h(1:nn))
      ALLOCATE(d3h(1:nn))
      ALLOCATE(d4h(1:nn))

      ALLOCATE(h_fft(1:kend))
      ALLOCATE(dummy_c(1:kend))

      CALL dfftw_plan_dft_r2c_1d(p_fw0, nn, hh,        h_fft, FFTW_MEASURE+FFTW_PRESERVE_INPUT)
      CALL dfftw_plan_dft_c2r_1d(p_bw1, nn, dummy_c,   d1h,   FFTW_MEASURE)
      CALL dfftw_plan_dft_c2r_1d(p_bw2, nn, dummy_c,   d2h,   FFTW_MEASURE)
      CALL dfftw_plan_dft_c2r_1d(p_bw3, nn, dummy_c,   d3h,   FFTW_MEASURE)
      CALL dfftw_plan_dft_c2r_1d(p_bw4, nn, dummy_c,   d4h,   FFTW_MEASURE)

      RETURN
    END SUBROUTINE create_fftw_plan

END Module fft_coeff

!----------------------------------------------------------------------
!----------------------------------------------------------------------
!  ks : Kuramoto-Sivashinski equation
!
!   h_t - c*h_x + h*h_x + h_xx + nu*h_xxxx=0, x\in[-pi,pi]
!   with periodic boundary conditions
!
!  Only derivatives are calculated using FFTW in the Fourier space
!
!----------------------------------------------------------------------
SUBROUTINE FUNC(NDIM,U,ICP,PAR,IJAC,F,DFDU,DFDP)
  USE nrtype
  USE fft_coeff
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: NDIM, IJAC, ICP(*)
  REAL(DP), INTENT(IN) :: U(NDIM), PAR(*)
  REAL(DP), INTENT(OUT) :: F(NDIM), DFDU(NDIM,NDIM), DFDP(NDIM,*)

  INTEGER :: nn, kk
  REAL(DP) :: nu, cc, tmp, dx
     
  nn = NDIM

! Initialisation: calculate fourier modes
  IF(fft_coeff_flag.EQ.0) CALL calculate_modes_r2c(nn, TWOPI)

! Initilize FFTW and create plans
  IF(fft_plan_flag.EQ.0) CALL create_fftw_plan(nn)

  cc = U(1)
  DO kk = 2, nn
     hh(kk) = U(kk)
  END DO
  hh(1) = -SUM(hh(2:nn))

! Compute h_fft from hh:
  CALL dfftw_execute_dft_r2c(p_fw0, hh, h_fft)

! Compute d1h (h_x) ~ d4h (h_xxxx):
  DO kk=1, nend
     dummy_c(kk) = coef_d1(kk)*h_fft(kk)/DBLE(nn)
  END DO
  CALL dfftw_execute_dft_c2r(p_bw1, dummy_c, d1h)
  DO kk=1, nend
     dummy_c(kk) = coef_d2(kk)*h_fft(kk)/DBLE(nn)
  END DO
  CALL dfftw_execute_dft_c2r(p_bw2, dummy_c, d2h)
  DO kk=1, nend
     dummy_c(kk) = coef_d3(kk)*h_fft(kk)/DBLE(nn)
  END DO
  CALL dfftw_execute_dft_c2r(p_bw3, dummy_c, d3h)
  DO kk=1, nend
     dummy_c(kk) = coef_d4(kk)*h_fft(kk)/DBLE(nn)
  END DO
  CALL dfftw_execute_dft_c2r(p_bw4, dummy_c, d4h)

  nu   = PAR(1)

  DO kk=2, NDIM
     F(kk) = cc*d1h(kk) - hh(kk)*d1h(kk) - d2h(kk) - nu*d4h(kk)
  END DO

! Condition for fixing the translational mode

!  F(1) = d1h(1)
  tmp = 0.0d0
  dx = TWOPI/DBLE(nn)
  DO kk=1, NDIM
     tmp = tmp + hh(kk)*DCOS(2.0d0*dx*DBLE(kk-1))
  END DO
  F(1) = tmp*dx

  RETURN
END SUBROUTINE FUNC

SUBROUTINE STPNT(NDIM,U,PAR,T)
  USE nrtype
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: NDIM
  REAL(DP), INTENT(INOUT) :: U(NDIM),PAR(*)
  REAL(DP), INTENT(IN) :: T
  INTEGER :: nn, kk
  REAL(DP) :: nu, dx

  nn = NDIM

! Initialise the equation parameters:
  nu  = 1.0d0

  PAR(1) = nu

!  par(11) = 1.0; ! Internal period, set to one and do not touch

! Initial condition
! U(1) is the traveling wave speed
  U(1) = 0.0d0
  dx = TWOPI/DBLE(nn)
  DO kk = 2, nn
     U(kk) = 1.0d-4*DSIN(dx*DBLE(kk-1))
  END DO

  RETURN
END SUBROUTINE STPNT

SUBROUTINE BCND
  RETURN
END SUBROUTINE BCND

SUBROUTINE ICND
  RETURN
END SUBROUTINE ICND

SUBROUTINE PVLS(NDIM,U,PAR)
  USE nrtype
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: NDIM
  REAL(DP), INTENT(IN) :: U(NDIM)
  REAL(DP), INTENT(INOUT) :: PAR(*)

  REAL(DP), EXTERNAL :: GETP,GETU2
  REAL(DP) :: tmp, tmp2
  INTEGER :: kk
!---------------------------------------------------------------------- 
! NOTE : 
! parameters set in this subroutine should be considered as ``solution 
! measures'' and be used for output purposes only.
! 
! they should never be used as `true'' continuation parameters. 
!
! they may, however, be added as ``over-specified parameters'' in the 
! parameter list associated with the auto-constant nicp, in order to 
! print their values on the screen and in the ``p.xxx file.
!
! they may also appear in the list associated with auto-constant nuzr.
!
!---------------------------------------------------------------------- 
! for algebraic problems the argument u is, as usual, the state vector.
! for differential equations the argument u represents the approximate 
! solution on the entire interval [0,1]. in this case its values must 
! be accessed indirectly by calls to getp, as illustrated below.
!---------------------------------------------------------------------- 
!
! set par(2) equal to the l2-norm of u(1)
!       par(2)=getp('nrm',1,u)
!
! set par(3) equal to the minimum of u(2)
!       par(3)=getp('min',2,u)
!
! set par(4) equal to the value of u(2) at the left boundary.
!       par(4)=getp('bv0',2,u)
!
! set par(5) equal to the pseudo-arclength step size used.
!       par(5)=getp('stp',1,u)
!
!---------------------------------------------------------------------- 
! the first argument of getp may be one of the following:
!        'nrm' (l2-norm),     'max' (maximum),
!        'int' (integral),    'bv0 (left boundary value),
!        'min' (minimum),     'bv1' (right boundary value).
!
! also available are
!   'stp' (pseudo-arclength step size used).
!   'fld' (`fold function', which vanishes at folds).
!   'bif' (`bifurcation function', which vanishes at singular points).
!   'hbf' (`hopf function'; which vanishes at hopf points).
!   'spb' ( function which vanishes at secondary periodic bifurcations).
!   'eig' ( obtains eigenvalues and floquet multipliers:
!            index 1, 3, 5,... obtain the real parts of ev's
!              and 2, 4, 6,... the imaginary parts
!            eigenvalues are ordered by real part.
!            floquet multipliers are ordered by distance from |z|=1)
!   'nbc' ( gets the active nbc).
!   'nint' ( gets the active nint).
!   'ntst' ( gets the active ntst).
!   'ncol' ( gets the active ncol).
!   'ndim' ( get the dimension of the extended system).
!   'mxt' ( gets t value for maximum of component).
!   'mnt' ( gets t value for minumum of component).
!   'win' ( gets weight for integral; 0<=i<=ncol).
!   'dtm' ( gets value for delta t array 1<=i<=ntst).
!---------------------------------------------------------------------- 

  PAR(6) = getp('INT', 1, U)

  tmp = 0.0d0
  tmp2 = 0.0d0
  DO kk=2, NDIM
     tmp = tmp + (getp('INT', kk, U))**2
     tmp2 = tmp2 + getp('INT', kk, U)
  END DO
  tmp = tmp + tmp2**2
  PAR(7) = DSQRT(tmp/DBLE(NDIM))

  RETURN
END SUBROUTINE PVLS

SUBROUTINE FOPT
  RETURN
END SUBROUTINE FOPT
