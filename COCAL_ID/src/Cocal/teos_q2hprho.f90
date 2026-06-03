subroutine COCAL_ID_teos_q2hprho(q,h,pre,rho,ened)
!

use COCAL_ID_def_teos_parameter
use COCAL_ID_phys_constant, only: long
implicit none
!
  real(8), intent(inout) :: q
  real(8), intent(out)   :: h, pre, rho, ened
  integer, save          :: iphase = 1
  integer                :: i0, i, ii
  real(long), external   :: COCAL_ID_lagint_4th
  real(long)             :: x4(4), f4(4), qloc
  real(8)                :: hin, qin, abin, abct, fac1, fac2, fack, small
!
  if (q <= qi(0)) then
    iphase = 1

    hin  = hi(0)
    qin  = qi(0)
    abin = gamma_crust
    abct = kappa_crust
    fac1 = 1.0d0/(abin - 1.0d0)
    fac2 = abin/(abin - 1.0d0)
    fack = abct**(-fac1)
    small = 1.0d-60
    if (q <= small) q = small
    h = hin + fac2*(q - qin)
    if (h <= 1.0d0) h = 1.0d0
    pre = fack*q**fac2
    rho = fack*q**fac1
    ened = rho*h - pre
    return
  end if

  if (iphase < 1 .or. iphase > nphase) iphase = 1

  ii = iphase
  if (q < qi(iphase)) then
    do i = ii, 1, -1
      if ((q >= qi(i-1)) .and. (q <= qi(i))) then
        iphase = i
        exit
      end if
    end do
  else
    do i = ii, nphase
      if ((q >= qi(i-1)) .and. (q <= qi(i))) then
        iphase = i
        exit
      end if
    end do
  end if

  i0 = min0(max0(iphase-2,0),nphase-3)
  qloc = q

  x4(1:4) = qi(i0:i0+3)
  f4(1:4) = hi(i0:i0+3)
  h       = COCAL_ID_lagint_4th(x4,f4,qloc)
  if (h <= 1.0d0) h = 1.0d0

  f4(1:4) = prei(i0:i0+3)
  pre     = COCAL_ID_lagint_4th(x4,f4,qloc)

  f4(1:4) = rhoi(i0:i0+3)
  rho     = COCAL_ID_lagint_4th(x4,f4,qloc)

  f4(1:4) = enei(i0:i0+3)
  ened    = COCAL_ID_lagint_4th(x4,f4,qloc)
!
end subroutine COCAL_ID_teos_q2hprho
