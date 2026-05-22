subroutine COCAL_ID_gr2cgr_4th_setup(xc,yc,zc,irgex4,itgex4,ipgex4,wr,wth,wphi)

use COCAL_ID_phys_constant, only : long, pi
use COCAL_ID_grid_parameter, only : nrg
implicit none
  real(long), intent(in)  :: xc, yc, zc
  integer, intent(out)    :: irgex4(4,4,4), itgex4(4,4,4), ipgex4(4,4,4)
  real(long), intent(out) :: wr(4), wth(4), wphi(4)
  real(long) :: rc, thc, phic, varpic
!
  rc     = dsqrt(dabs(xc**2 + yc**2 + zc**2))
  varpic = dsqrt(dabs(xc**2 + yc**2))
  thc  = dmod(2.0d0*pi + datan2(varpic,zc),2.0d0*pi)
  phic = dmod(2.0d0*pi + datan2(    yc,xc),2.0d0*pi)
!
  call COCAL_ID_cgr_4th_setup_from_spherical(rc,thc,phic,nrg,irgex4,itgex4,ipgex4,wr,wth,wphi)
!
end subroutine COCAL_ID_gr2cgr_4th_setup

subroutine COCAL_ID_fl2cgr_4th_setup(xc,yc,zc,rs,outside,irgex4,itgex4,ipgex4,wr,wth,wphi)

use COCAL_ID_phys_constant, only : long, pi
use COCAL_ID_grid_parameter, only : nrf
use COCAL_ID_coordinate_grav_r, only : rg
use COCAL_ID_interface_interpo_lag4th_2Dsurf
implicit none
  real(long), intent(in)  :: xc, yc, zc
  real(long), pointer     :: rs(:,:)
  logical, intent(out)    :: outside
  integer, intent(out)    :: irgex4(4,4,4), itgex4(4,4,4), ipgex4(4,4,4)
  real(long), intent(out) :: wr(4), wth(4), wphi(4)
  real(long) :: rsca, rc_gr, rc, thc, phic, varpic
!
  outside = .false.
  rc_gr  = dsqrt(dabs(xc**2 + yc**2 + zc**2))
  varpic = dsqrt(dabs(xc**2 + yc**2))
  thc  = dmod(2.0d0*pi + datan2(varpic,zc),2.0d0*pi)
  phic = dmod(2.0d0*pi + datan2(    yc,xc),2.0d0*pi)
!
  call COCAL_ID_interpo_lag4th_2Dsurf(rsca,rs,thc,phic)
  rc = rc_gr/rsca
  if (rc.gt.rg(nrf)) then
    outside = .true.
    return
  end if
!
  call COCAL_ID_cgr_4th_setup_from_spherical(rc,thc,phic,nrf,irgex4,itgex4,ipgex4,wr,wth,wphi)
!
end subroutine COCAL_ID_fl2cgr_4th_setup

subroutine COCAL_ID_cgr_4th_setup_from_spherical(rc,thc,phic,nr,irgex4,itgex4,ipgex4,wr,wth,wphi)

use COCAL_ID_phys_constant, only : long
use COCAL_ID_grid_parameter, only : ntg, npg
use COCAL_ID_coordinate_grav_extended
implicit none
  real(long), intent(in)  :: rc, thc, phic
  integer, intent(in)     :: nr
  integer, intent(out)    :: irgex4(4,4,4), itgex4(4,4,4), ipgex4(4,4,4)
  real(long), intent(out) :: wr(4), wth(4), wphi(4)
  real(long) :: r4(4), th4(4), phi4(4)
  integer :: irg, itg, ipg, ir0, it0, ip0
  integer :: ii, jj, kk, irg0, itg0, ipg0
!
  do irg = 0, nr+1
    if (rc.lt.rgex(irg).and.rc.ge.rgex(irg-1)) ir0 = min0(irg-2,nr-3)
  end do
  do itg = 0, ntg+1
    if (thc.lt.thgex(itg).and.thc.ge.thgex(itg-1)) it0 = itg-2
  end do
  do ipg = 0, npg+1
    if (phic.lt.phigex(ipg).and.phic.ge.phigex(ipg-1)) ip0 = ipg-2
  end do
!
  do ii = 1, 4
    irg0 = ir0 + ii - 1
    itg0 = it0 + ii - 1
    ipg0 = ip0 + ii - 1
    r4(ii) = rgex(irg0)
    th4(ii) = thgex(itg0)
    phi4(ii) = phigex(ipg0)
  end do
!
  call COCAL_ID_lagint_4th_weights(r4,rc,wr)
  call COCAL_ID_lagint_4th_weights(th4,thc,wth)
  call COCAL_ID_lagint_4th_weights(phi4,phic,wphi)
!
  do kk = 1, 4
    ipg0 = ip0 + kk - 1
    do jj = 1, 4
      itg0 = it0 + jj - 1
      do ii = 1, 4
        irg0 = ir0 + ii - 1
        irgex4(ii,jj,kk) = irgex_r(irg0)
        itgex4(ii,jj,kk) = itgex_r(itgex_th(itg0),irg0)
        ipgex4(ii,jj,kk) = ipgex_r(ipgex_th(ipgex_phi(ipg0),itg0),irg0)
      end do
    end do
  end do
!
end subroutine COCAL_ID_cgr_4th_setup_from_spherical

subroutine COCAL_ID_cgr_4th_apply(fnc,cfn,irgex4,itgex4,ipgex4,wr,wth,wphi)

use COCAL_ID_phys_constant, only : long
use COCAL_ID_coordinate_grav_extended
implicit none
  real(long), pointer     :: fnc(:,:,:)
  real(long), intent(out) :: cfn
  integer, intent(in)     :: irgex4(4,4,4), itgex4(4,4,4), ipgex4(4,4,4)
  real(long), intent(in)  :: wr(4), wth(4), wphi(4)
  real(long) :: fr4(4), ft4(4), fp4(4)
  integer :: ii, jj, kk
  real(long), external :: COCAL_ID_lagint_4th_apply
!
  cfn = 0.0d0
!
  do kk = 1, 4
    do jj = 1, 4
      do ii = 1, 4
        fr4(ii) = fnc(irgex4(ii,jj,kk),itgex4(ii,jj,kk),ipgex4(ii,jj,kk))
      end do
      ft4(jj) = COCAL_ID_lagint_4th_apply(wr,fr4)
    end do
    fp4(kk) = COCAL_ID_lagint_4th_apply(wth,ft4)
  end do
  cfn = COCAL_ID_lagint_4th_apply(wphi,fp4)
!
end subroutine COCAL_ID_cgr_4th_apply
