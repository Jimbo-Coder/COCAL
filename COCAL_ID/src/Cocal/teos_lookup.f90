subroutine COCAL_ID_teos_lookup(qp, qpar, iphase)
!

use COCAL_ID_phys_constant
use COCAL_ID_def_teos_parameter
implicit none
!
  real(8), intent(in)  :: qp, qpar(0:nnteos)
  real(8)              :: det
  integer, intent(out) :: iphase
  integer              :: ii
!
! --  Monotonically increasing qpar is assumed.
!
  iphase = 1
  do ii = 1, nphase
    det = (qp-qpar(ii))*(qp-qpar(ii-1))
    if (det <= 0.0d0) then
      iphase = ii
      exit
    end if
  end do
!
end subroutine COCAL_ID_teos_lookup
