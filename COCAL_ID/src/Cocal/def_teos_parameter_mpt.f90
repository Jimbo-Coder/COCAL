module COCAL_ID_def_teos_parameter_mpt
use COCAL_ID_phys_constant         !nnteos, nmpt
implicit none
  real(8) :: rhoi_(0:nnteos,nmpt), qi_(0:nnteos,nmpt), hi_(0:nnteos,nmpt)
  real(8) :: prei_(0:nnteos,nmpt), enei_(0:nnteos,nmpt)
  real(8) :: rhocgs_(0:nnteos,nmpt), precgs_(0:nnteos,nmpt), enecgs_(0:nnteos,nmpt)
  real(8) :: def_teos_param_real_(0:20,nmpt)
  integer :: def_teos_param_int_(0:10,nmpt)
end module COCAL_ID_def_teos_parameter_mpt
