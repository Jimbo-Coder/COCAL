subroutine COCAL_ID_interpo_fl2cgr_4th_export(fnc,cfn,xc,yc,zc,rs)

use COCAL_ID_phys_constant, only : long
use COCAL_ID_interface_modules_cartesian, ignore_me => COCAL_ID_interpo_fl2cgr_4th_export
implicit none
  real(long), pointer     :: fnc(:,:,:), rs(:,:)
  real(long), intent(out) :: cfn
  real(long) ::  xc, yc, zc
  real(long) ::  wr(4), wth(4), wphi(4)
  integer :: irgex4(4,4,4), itgex4(4,4,4), ipgex4(4,4,4)
  logical :: outside
!
! --- Interpolation from the central spherical grid 
! --- to a Cartesian grid point using 4th order Lagrange formula.
!
  cfn = 0.0d0
  call COCAL_ID_fl2cgr_4th_setup(xc,yc,zc,rs,outside,irgex4,itgex4,ipgex4,wr,wth,wphi)
  if (outside) return
  call COCAL_ID_cgr_4th_apply(fnc,cfn,irgex4,itgex4,ipgex4,wr,wth,wphi)
!
end subroutine COCAL_ID_interpo_fl2cgr_4th_export
