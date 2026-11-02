subroutine COCAL_ID_coordinate_patch_kit_grav_grid_coc2cac_mpt(igrid,impt)

use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
use COCAL_ID_grid_parameter
use COCAL_ID_coordinate_grav_r
use COCAL_ID_coordinate_grav_phi
use COCAL_ID_coordinate_grav_theta
use COCAL_ID_coordinate_grav_extended
use COCAL_ID_trigonometry_grav_theta
use COCAL_ID_trigonometry_grav_phi
implicit none
  integer, intent(in) :: igrid, impt
  character(len=256) :: message
  real(long), parameter :: grid_tol = 1.0d-10

  ! Regenerate patch impt's grid from its rnspar parameters.
  select case (igrid)
  case (1)
    call grid_r
  case (2)
    call grid_r_bhex('eBH')
  case (3)
    call grid_r_bns
  case (4)
    call grid_r_bqs
  case default
    call invalid_grid('Choose legacy grid 1, 2, 3 or 4.')
  end select
  call grid_theta
  call grid_phi

  if (impt <= 2) then
    if (nrf < 3 .or. nrf > nrg .or. ntf /= ntg .or. npf /= npg) &
      call invalid_grid('Fluid grid must use the gravitational angular grid and at least four radial points.')
    if (.not. ieee_is_finite(r_surf) .or. r_surf <= 0.0d0) &
      call invalid_grid('Invalid r_surf.')
    if (abs(rg(nrf)-r_surf) > grid_tol*max(1.0d0,abs(r_surf))) then
      write(message,'(a,i0,a,es20.12,a,es20.12)') &
        'rg(nrf) differs from r_surf: nrf=', nrf, '; rg(nrf)=', rg(nrf), '; r_surf=', r_surf
      call invalid_grid(trim(message))
    end if
  end if
  call trig_grav_theta
  call trig_grav_phi
  call grid_extended

contains
  subroutine invalid_grid(message)
    character(len=*), intent(in) :: message
    write(*,'(a,i0,a)') 'COCAL_ID: BNS patch ', impt, ': '//message
    error stop 'Invalid BNS grid'
  end subroutine invalid_grid
end subroutine COCAL_ID_coordinate_patch_kit_grav_grid_coc2cac_mpt
