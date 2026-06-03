subroutine COCAL_ID_copy_def_teos_parameter_from_mpt(impt)

use COCAL_ID_def_teos_parameter
use COCAL_ID_def_teos_parameter_mpt
use COCAL_ID_phys_constant, only : nnteos
use COCAL_ID_copy_array_static_2dto1d_mpt
implicit none
  integer :: i, impt
!
  call copy_arraystatic_2dto1d_mpt(impt, rhoi_,   rhoi,   0, nnteos)
  call copy_arraystatic_2dto1d_mpt(impt, qi_,     qi,     0, nnteos)
  call copy_arraystatic_2dto1d_mpt(impt, hi_,     hi,     0, nnteos)
  call copy_arraystatic_2dto1d_mpt(impt, prei_,   prei,   0, nnteos)
  call copy_arraystatic_2dto1d_mpt(impt, enei_,   enei,   0, nnteos)
  call copy_arraystatic_2dto1d_mpt(impt, rhocgs_, rhocgs, 0, nnteos)
  call copy_arraystatic_2dto1d_mpt(impt, precgs_, precgs, 0, nnteos)
  call copy_arraystatic_2dto1d_mpt(impt, enecgs_, enecgs, 0, nnteos)
!
  i=0
  i=i+1; rhoini_cgs  = def_teos_param_real_(i,impt)
  i=i+1; rhoini_gcm1 = def_teos_param_real_(i,impt)
  i=i+1; emdini_gcm1 = def_teos_param_real_(i,impt)
  i=i+1; kappa_crust = def_teos_param_real_(i,impt)
  i=i+1; gamma_crust = def_teos_param_real_(i,impt)
!
  i=0
  i=i+1; nphase = def_teos_param_int_(i,impt)
!
end subroutine COCAL_ID_copy_def_teos_parameter_from_mpt
