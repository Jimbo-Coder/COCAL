subroutine COCAL_ID_calc_bht_excision_radius

use COCAL_ID_phys_constant, only : long
use COCAL_ID_grid_parameter, only : rgin
use COCAL_ID_def_bh_parameter, only : mass_bh, spin_bh
implicit none
  real(long) :: spin_over_mass

  spin_over_mass = spin_bh/mass_bh
  rgin = mass_bh*(1.0d0 + 0.8d0*dsqrt(max(1.0d0-spin_over_mass**2, 0.0d0)))
end subroutine COCAL_ID_calc_bht_excision_radius
