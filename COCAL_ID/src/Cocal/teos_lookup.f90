subroutine COCAL_ID_teos_lookup(qp, qpar, iphase)
!

use COCAL_ID_phys_constant
use COCAL_ID_def_teos_parameter
implicit none
!
  real(8), intent(in)  :: qp, qpar(0:nnteos)
  integer, intent(out) :: iphase
  integer              :: ilo, ihi, imid
!
! --  Monotonically increasing qpar is assumed.
!
  if (qp <= qpar(0)) then
    iphase = 1
    return
  else if (qp >= qpar(nphase)) then
    iphase = nphase
    return
  end if

  ilo = 0
  ihi = nphase
  do while (ihi - ilo > 1)
    imid = (ilo + ihi)/2
    if (qp <= qpar(imid)) then
      ihi = imid
    else
      ilo = imid
    end if
  end do
  iphase = ihi
!
end subroutine COCAL_ID_teos_lookup
