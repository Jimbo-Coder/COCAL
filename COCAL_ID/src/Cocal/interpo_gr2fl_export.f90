subroutine COCAL_ID_interpo_gr2fl_export(grv,flv,rs)

use COCAL_ID_phys_constant, only : long
use COCAL_ID_grid_parameter, only : nrg, ntg, npg, nrf, ntf, npf
use COCAL_ID_coordinate_grav_r, only : rg!
!  use def_matter, only : rs
implicit none
  real(long), external :: COCAL_ID_lagint_4th
  real(long), pointer :: grv(:,:,:), flv(:,:,:), rs(:,:)
  real(long) :: x(4), f(4)
  real(long) :: rrff,   small = 1.0d-14
  integer :: irg, irf, itf, ipf, ir0
!
  flv(0:nrf,0:ntf,0:npf) = 0.0d0
!
  do ipf = 0, npf
    do itf = 0, ntf
      do irf = 0, nrf
        rrff = rs(itf,ipf)*rg(irf)
        if (rrff < rg(0) .or. rrff > rg(nrg)) error stop 'COCAL_ID: fluid surface exceeds gravitational grid'
        ir0 = nrg-3
        do irg = 0, nrg
          if (rrff.le.rg(irg)) then 
            ir0 = min0(max0(0,irg-2),nrg-3)
            exit
          end if
        end do
        x(1:4) = rg(ir0:ir0+3)
        f(1:4) = grv(ir0:ir0+3,itf,ipf)
        flv(irf,itf,ipf) = COCAL_ID_lagint_4th(x,f,rrff)
      end do
    end do
  end do
!
end subroutine COCAL_ID_interpo_gr2fl_export
