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

  INTEGER(KIND=8) :: p_bw0, p_fw2
  REAL(KIND=8), ALLOCATABLE, DIMENSION(:) :: hh, hh2
  COMPLEX(KIND=8), ALLOCATABLE, DIMENSION(:) :: h_fft, h2_fft

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
      ALLOCATE(hh2(1:nn))

      ALLOCATE(h_fft(1:kend))
      ALLOCATE(h2_fft(1:kend))

      CALL dfftw_plan_dft_c2r_1d(p_bw0, nn, h_fft,  hh,      FFTW_MEASURE+FFTW_PRESERVE_INPUT)
      CALL dfftw_plan_dft_r2c_1d(p_fw2, nn, hh2,    h2_fft,  FFTW_MEASURE)

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
  COMPLEX(KIND=8) :: F_fft(NDIM)

  INTEGER :: nn, kk
  REAL(DP) :: nu, cc
     
  nn = NDIM

! Initialisation: calculate fourier modes
  IF(fft_coeff_flag.EQ.0) CALL calculate_modes_r2c(nn, TWOPI)

! Initilize FFTW and create plans
  IF(fft_plan_flag.EQ.0) CALL create_fftw_plan(nn)

! obtain cc and h_fft
  cc = U(1)
  h_fft(1) = DCMPLX(0.0d0, 0.0d0)
  IF(MOD(nn,2).EQ.0)THEN
     DO kk = 2, nn/2
        h_fft(kk) = DCMPLX(U(2*kk-2), U(2*kk-1))
     END DO
     h_fft(nend) = DCMPLX(U(nn), 0.0d0)
  ELSE
     DO kk = 2, (nn+1)/2
        h_fft(kk) = DCMPLX(U(2*kk-2), U(2*kk-1))
     END DO
  END IF

! Compute hh from h_fft:
  CALL dfftw_execute_dft_c2r(p_bw0, h_fft, hh)
  DO kk=1, nn
     hh(kk) = hh(kk)/DBLE(nn)
  END DO

! Compute h^2
  DO kk=1, nn
     hh2(kk) = hh(kk)**2
  END DO

! Compute h2_fft: Fourier(h^2)
  CALL dfftw_execute_dft_r2c(p_fw2, hh2, h2_fft)

  nu   = PAR(1)

  DO kk=1, nend
     F_fft(kk) = coef_d1(kk)*(cc*h_fft(kk) - 0.5d0*h2_fft(kk)) &
          - (coef_d2(kk)+nu*coef_d4(kk))*h_fft(kk)
  END DO

  F(1) = U(2)
  IF(MOD(nn,2).EQ.0)THEN
     DO kk = 2, nn/2
        F(2*kk-2) = DREAL(F_fft(kk))
        F(2*kk-1) = DIMAG(F_fft(kk))
     END DO
     F(nn) = DREAL(F_fft(nend))
  ELSE
     DO kk = 2, (nn+1)/2
        F(2*kk-2) = DREAL(F_fft(kk))
        F(2*kk-1) = DIMAG(F_fft(kk))
     END DO
  END IF

  RETURN
END SUBROUTINE FUNC

SUBROUTINE STPNT(NDIM,U,PAR,T)
  USE nrtype
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: NDIM
  REAL(DP), INTENT(INOUT) :: U(NDIM),PAR(*)
  REAL(DP), INTENT(IN) :: T
  INTEGER :: nn, kk
  REAL(DP) :: nu

  nn = NDIM

! Initialise the equation parameters:
!  nu = 1.0d0
  nu = 1.0d0/(2.0d0**2)
!  nu = 1.0d0/(4.0d0**2)

  PAR(1) = nu

!  par(11) = 1.0; ! Internal period, set to one and do not touch

! Initial condition
! U(1) is the traveling wave speed
  DO kk = 1, NDIM
     U(kk) = 0.0d0
  END DO

!  U(3) = 1.0d-2
  U(5) = 1.0d-1
!  U(7) = 1.0d-2
!  U(9) = 1.0d-2

  RETURN
END SUBROUTINE STPNT

SUBROUTINE BCND
  RETURN
END SUBROUTINE BCND

SUBROUTINE ICND
  RETURN
END SUBROUTINE ICND

SUBROUTINE PVLS(NDIM,U,PAR)
!     ---------- ----

  USE nrtype
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: NDIM
  REAL(DP), INTENT(IN) :: U(NDIM)
  REAL(DP), INTENT(INOUT) :: PAR(*)

  REAL(DP), EXTERNAL :: GETP
  INTEGER :: k
  REAL(DP) :: sumall

  sumall = 0.0d0
  DO k=2, NDIM
     sumall = sumall + 2.0d0*(getp("NRM", k, U)**2)
  END DO
  PAR(45) = DSQRT(sumall)/DBLE(NDIM)

  RETURN
END SUBROUTINE PVLS

SUBROUTINE FOPT
  RETURN
END SUBROUTINE FOPT
