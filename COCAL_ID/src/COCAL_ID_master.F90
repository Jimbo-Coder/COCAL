#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Functions.h"
#include "cctk_Parameters.h"

#define module_use use

! COCAL's own routines (Cocal/), as exported by PCOCAL.
#include "Cocal/phys_constant.f90"
#include "Cocal/def_matter_parameter.f90"
#include "Cocal/def_quantities.f90"
#include "Cocal/def_bh_parameter.f90"
#include "Cocal/def_bht_parameter.f90"
#include "Cocal/def_peos_parameter.f90"
#include "Cocal/def_teos_parameter.f90"

#include "Cocal/def_binary_parameter.f90"
#include "Cocal/def_peos_parameter_mpt.f90"
#include "Cocal/def_teos_parameter_mpt.f90"
#include "Cocal/make_char1_array_2d.f90"
#include "Cocal/make_char2_array_2d.f90"
#include "Cocal/make_int_array_2d.f90"

#include "Cocal/make_array_2d.f90"
#include "Cocal/grid_parameter.f90"
#include "Cocal/coordinate_grav_r.f90"
#include "Cocal/coordinate_grav_phi.f90"
#include "Cocal/coordinate_grav_theta.f90"
#include "Cocal/coordinate_grav_extended.f90"
#include "Cocal/trigonometry_grav_theta.f90"
#include "Cocal/trigonometry_grav_phi.f90"

#include "Cocal/grid_parameter_binary_excision.f90"

#include "Cocal/interface_IO_input_CF_grav_export.f90"
#include "Cocal/interface_IO_input_WL_grav_export_hij.f90"
#include "Cocal/interface_IO_input_grav_export_Kij.f90"
#include "Cocal/interface_IO_input_CF_star_export.f90"
#include "Cocal/interface_IO_input_matter_BHT_export.f90"
#include "Cocal/interface_invhij_WL_export.f90"
#include "Cocal/interface_index_vec_down2up_export.f90"
#include "Cocal/interface_IO_input_grav_export_Ai.f90"
#include "Cocal/interface_IO_input_grav_export_Faraday.f90"
#include "Cocal/interface_IO_input_star4ve_export.f90"

#include "Cocal/interface_IO_input_CF_surf_export.f90"
#include "Cocal/interface_IO_input_CF_flco_export.f90"
#include "Cocal/interface_IO_input_CF_flir_export.f90"
#include "Cocal/interface_IO_input_CF_flsp_export.f90"
#include "Cocal/interface_IO_input_gradvep_export.f90"

#include "Cocal/IO_input_CF_grav_export.f90"
#include "Cocal/IO_input_WL_grav_export_hij.f90"
#include "Cocal/IO_input_grav_export_Kij.f90"
#include "Cocal/IO_input_CF_star_export.f90"
#include "Cocal/IO_input_matter_BHT_export.f90"
#include "Cocal/invhij_WL_export.f90"
#include "Cocal/index_vec_down2up_export.f90"
#include "Cocal/IO_input_grav_export_Ai.f90"
#include "Cocal/IO_input_grav_export_Faraday.f90"
#include "Cocal/IO_input_star4ve_export.f90"
#include "Cocal/interpo_lag4th_weights.f90"
#include "Cocal/lagint_4th.f90"
#include "Cocal/peos_initialize.f90"
#include "Cocal/peos_q2hprho.f90"
#include "Cocal/peos_lookup.f90"
#include "Cocal/teos_initialize.f90"
#include "Cocal/teos_q2hprho.f90"
#include "Cocal/teos_lookup.f90"
#include "Cocal/read_bht_parameter.f90"
#include "Cocal/calc_bht_excision_radius.f90"
#include "Cocal/copy_array_static_1dto2d_mpt.f90"
#include "Cocal/copy_array_static_2dto1d_mpt.f90"
#include "Cocal/grid_parameter_mpt.f90"
#include "Cocal/grid_parameter_binary_excision_mpt.f90"
#include "Cocal/allocate_grid_parameter_mpt.f90"
#include "Cocal/allocate_grid_parameter_binary_excision_mpt.f90"
#include "Cocal/read_parameter_mpt.f90"
#include "Cocal/read_parameter_binary_excision_mpt.f90"
#include "Cocal/read_surf_parameter_mpt.f90"
#include "Cocal/copy_grid_parameter_to_mpt.f90"
#include "Cocal/copy_grid_parameter_from_mpt.f90"
#include "Cocal/copy_grid_parameter_binary_excision_from_mpt.f90"
#include "Cocal/copy_grid_parameter_binary_excision_to_mpt.f90"
#include "Cocal/copy_def_peos_parameter_to_mpt.f90"
#include "Cocal/copy_def_peos_parameter_from_mpt.f90"
#include "Cocal/copy_def_teos_parameter_to_mpt.f90"
#include "Cocal/copy_def_teos_parameter_from_mpt.f90"
#include "Cocal/IO_input_CF_surf_export.f90"
#include "Cocal/IO_input_CF_flco_export.f90"
#include "Cocal/IO_input_CF_flir_export.f90"
#include "Cocal/IO_input_CF_flsp_export.f90"
#include "Cocal/peos_initialize_mpt.f90"
#include "Cocal/coordinate_patch_kit_grav_grid_coc2cac_mpt.f90"
#include "Cocal/set_allocate_size_mpt.f90"
#include "Cocal/IO_input_gradvep_export.f90"

! The reader (cocal_id_*): the shared domain, state and holes, one file per
! kind of ID, and the evaluation at a point.
#include "cocal_id_node_shared.f90"
#include "cocal_id_spline.f90"
#include "cocal_id_domain.f90"
#include "cocal_id_state.f90"
#include "cocal_id_holes.f90"
#include "cocal_id_rns.F90"
#include "cocal_id_bns.F90"
#include "cocal_id_bht.F90"
#include "cocal_id_bbh.f90"
#include "cocal_id_eval.f90"


! Read the ID once per MPI rank, collectively over each node.
subroutine COCAL_ID_read_data(CCTK_ARGUMENTS)
  use COCAL_ID_state
  use COCAL_ID_bbh, only : bbh_error, COCAL_ID_read_bbh
  implicit none
  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_PARAMETERS
  DECLARE_CCTK_FUNCTIONS
  character(400) :: dir, msg
  integer :: dir_len
  logical :: ok, three_patch

  if (cocal_kind /= 0) return
  call CCTK_FortranString(dir_len, COCAL_ID_PathToID, dir)
  cocal_tabulated_eos = COCAL_ID_read_tabulated_eos == 1
  if (CCTK_EQUALS(initial_data, "CocalRNS")) then
    cocal_kind = COCAL_RNS
    call COCAL_ID_read_rns(trim(dir))
  else if (CCTK_EQUALS(initial_data, "CocalBNS")) then
    cocal_kind = COCAL_BNS
    call COCAL_ID_read_bns(trim(dir))
  else if (CCTK_EQUALS(initial_data, "CocalBHT")) then
    cocal_kind = COCAL_BHT
    call COCAL_ID_read_bht(trim(dir))
  else
    cocal_kind = merge(COCAL_BBHT, COCAL_BBH, CCTK_EQUALS(initial_data, "CocalBBHT"))
    call COCAL_ID_read_bbh(trim(dir), ok, CCTK_EQUALS(COCAL_ID_bbh_fill_method, "smooth"), cocal_kind == COCAL_BBHT)
    if (.not. ok) call CCTK_WARN(CCTK_WARN_ABORT, trim(bbh_error))
  end if

  ! Patch blend zones of a binary's metric.
  three_patch = cocal_metric%npatches == 3
  if (three_patch .and. (cocal_metric%slab_halfwidth == 0 .or. cocal_metric%shell_hi == cocal_metric%shell_lo)) &
    call CCTK_WARN(CCTK_WARN_ALERT, "COCAL_ID: no usable patch overlap; using a hard patch switch there.")
  if (three_patch .and. COCAL_ID_verbose == 1) then
    write(msg, '(a,es12.4,a,2es12.4,a)') 'COCAL_ID patch blend: midplane |x-mid| <', cocal_metric%slab_halfwidth, &
      ', outer shell', cocal_metric%shell_lo, cocal_metric%shell_hi, ' (COCAL units)'
    call CCTK_INFO(trim(msg))
  end if
  if (COCAL_ID_verbose == 1) then
    write(msg, '(a,es20.12)') 'COCAL_ID: read the COCAL data; radi = ', cocal_radi
    call CCTK_INFO(trim(msg))
  end if
end subroutine COCAL_ID_read_data

! The values on the grid: the driver's own file (Carpet: the grid loop;
! CarpetX: a point evaluator for its C++ loops).
#include "cocal_id_grid.F90"

! Collective over the node.
subroutine COCAL_ID_free_data(CCTK_ARGUMENTS)
  use COCAL_ID_state
  use COCAL_ID_holes
  implicit none
  DECLARE_CCTK_ARGUMENTS

  call COCAL_ID_holes_free()
  call COCAL_ID_state_free()
end subroutine COCAL_ID_free_data
