subroutine COCAL_ID_calc_bht_excision_radius

use COCAL_ID_phys_constant, only : long, pi
use COCAL_ID_grid_parameter, only : rgin
use COCAL_ID_def_bh_parameter, only : mass_bh, spin_bh
implicit none
  real(long) :: drh, rh, aom, sqrt_one_minus_aom2

  if (mass_bh <= 0.0d0) stop "Invalid BHT parameters: mass_bh must be positive."

  aom = dabs(spin_bh/mass_bh)
  if (aom > 1.0d0) stop "Invalid BHT parameters: abs(spin_bh) exceeds mass_bh."
  sqrt_one_minus_aom2 = dsqrt(max(1.0d0-aom**2, 0.0d0))

  rh = mass_bh*(1.0d0 + sqrt_one_minus_aom2)
  drh = rh - spin_bh
!  rgin = spin_bh + 0.9*drh

!  rgin = mass_bh*( 1.0d0 + 0.2d0*dsqrt(1.0d0-aom**2) )
  rgin = mass_bh*( 1.0d0 + 0.8d0*sqrt_one_minus_aom2 )

!  open(120,file="rnspar_rgin.dat", status='unknown')
!  write(120,'(a16,1p,3e23.15)') "r horizon      =", rh 
!  write(120,'(a16,1p,3e23.15)') "Kerr parameter =", spin_bh
!  write(120,'(a16,1p,3e23.15)') "rgin           =", rgin 
!  close(120)
!
end subroutine COCAL_ID_calc_bht_excision_radius
