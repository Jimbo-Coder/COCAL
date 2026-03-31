!_____________________________________________________________________________
!
!    SHARED SUPPORT SOURCES FOR COCAL_ID SPLIT DRIVERS
!_____________________________________________________________________________
!

! Shared modules
#include "Cocal/phys_constant.f90"
#include "Cocal/def_matter_parameter.f90"
#include "Cocal/def_quantities.f90"
#include "Cocal/def_peos_parameter.f90"
#include "Cocal/make_char1_array_2d.f90"
#include "Cocal/make_char2_array_2d.f90"
#include "Cocal/make_int_array_2d.f90"
#include "Cocal/make_array_1d.f90"
#include "Cocal/make_array_2d.f90"
#include "Cocal/make_array_3d.f90"
#include "Cocal/grid_parameter.f90"
#include "Cocal/interface_modules_cartesian.f90"
#include "Cocal/coordinate_grav_r.f90"
#include "Cocal/coordinate_grav_phi.f90"
#include "Cocal/coordinate_grav_theta.f90"
#include "Cocal/coordinate_grav_extended.f90"
#include "Cocal/trigonometry_grav_theta.f90"
#include "Cocal/trigonometry_grav_phi.f90"
#include "Cocal/interface_IO_input_CF_grav_export.f90"
#include "Cocal/interface_IO_input_grav_export_Kij.f90"
#include "Cocal/interface_interpo_gr2fl_metric_CF_export.f90"

! Shared subroutines
#include "Cocal/IO_input_CF_grav_export.f90"
#include "Cocal/IO_input_grav_export_Kij.f90"
#include "Cocal/interpo_gr2fl_metric_CF_export.f90"
#include "Cocal/interpo_gr2cgr_4th.f90"
#include "Cocal/interpo_fl2cgr_4th_export.f90"
#include "Cocal/peos_initialize.f90"
#include "Cocal/peos_q2hprho.f90"
#include "Cocal/peos_lookup.f90"

