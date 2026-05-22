subroutine COCAL_ID_interpo_gr2cgr_4th(fnc,cfn,xc,yc,zc)

use COCAL_ID_phys_constant, only : long
use COCAL_ID_interface_modules_cartesian, ignore_me => COCAL_ID_interpo_gr2cgr_4th
implicit none
  real(long), pointer     :: fnc(:,:,:)
  real(long), intent(out) :: cfn
  real(long) ::  xc, yc, zc
  real(long) ::  wr(4), wth(4), wphi(4)
  integer :: irgex4(4,4,4), itgex4(4,4,4), ipgex4(4,4,4)
!
! --- Interpolation from the central spherical grid 
! --- to a Cartesian grid point using 4th order Lagrange formula.
!
  call COCAL_ID_gr2cgr_4th_setup(xc,yc,zc,irgex4,itgex4,ipgex4,wr,wth,wphi)
  call COCAL_ID_cgr_4th_apply(fnc,cfn,irgex4,itgex4,ipgex4,wr,wth,wphi)
!
end subroutine COCAL_ID_interpo_gr2cgr_4th
