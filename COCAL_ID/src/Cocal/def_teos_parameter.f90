module COCAL_ID_def_teos_parameter
use COCAL_ID_phys_constant         !nnteos
implicit none
  real(8) :: rhoi(0:nnteos), qi(0:nnteos), hi(0:nnteos)
  real(8) :: prei(0:nnteos), enei(0:nnteos)
  real(8) :: rhocgs(0:nnteos), precgs(0:nnteos), enecgs(0:nnteos)
  real(8) :: kappa_crust, gamma_crust
  real(8) :: rhoini_cgs, rhoini_gcm1, emdini_gcm1  !used in TOV solver
  integer :: nphase
end module COCAL_ID_def_teos_parameter
