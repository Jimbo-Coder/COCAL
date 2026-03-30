!_____________________________________________________________________________
!
!    CACTUS READER OF COCAL ROTATING STAR IN WAVELESS FORMALISM 
!_____________________________________________________________________________
!

#define module_use use

#include "Cocal/phys_constant.f90"!


#include "Cocal/def_binary_parameter.f90"! 
#include "Cocal/def_matter_parameter_mpt.f90"!
#include "Cocal/def_binary_parameter_mpt.f90"!

#include "Cocal/grid_parameter.f90"!
#include "Cocal/interface_modules_cartesian.f90"!
#include "Cocal/coordinate_grav_r.f90"!
#include "Cocal/coordinate_grav_phi.f90"
#include "Cocal/coordinate_grav_theta.f90"
#include "Cocal/coordinate_grav_extended.f90"!
#include "Cocal/trigonometry_grav_theta.f90"
#include "Cocal/trigonometry_grav_phi.f90"!

#include "Cocal/grid_parameter_binary_excision.f90"!

#include "Cocal/interface_IO_input_CF_grav_export.f90"!
#include "Cocal/interface_IO_input_grav_export_Kij.f90"!
#include "Cocal/interface_interpo_gr2fl_metric_CF_export.f90"!
#include "Cocal/interface_interpo_lag4th_2Dsurf.f90"!

#include "Cocal/interface_IO_input_CF_surf_export.f90"!
#include "Cocal/interface_IO_input_CF_flco_export.f90"!
#include "Cocal/interface_IO_input_CF_flir_export.f90"!
#include "Cocal/interface_IO_input_CF_flsp_export.f90"!
#include "Cocal/interface_IO_input_gradvep_export.f90"!


#include "Cocal/grid_parameter_mpt.f90"!
#include "Cocal/grid_parameter_binary_excision_mpt.f90"!
#include "Cocal/allocate_grid_parameter_mpt.f90"!
#include "Cocal/allocate_grid_parameter_binary_excision_mpt.f90"!
#include "Cocal/allocate_def_matter_parameter_mpt.f90"!
#include "Cocal/trigonometry_grav_phi_mpt.f90"!
#include "Cocal/trigonometry_grav_theta_mpt.f90"
#include "Cocal/coordinate_grav_r_mpt.f90"!
#include "Cocal/coordinate_grav_extended_mpt.f90"!
#include "Cocal/copy_grid_parameter_to_mpt.f90"!
#include "Cocal/copy_grid_parameter_from_mpt.f90"!
#include "Cocal/copy_grid_parameter_binary_excision_from_mpt.f90"!
#include "Cocal/copy_grid_parameter_binary_excision_to_mpt.f90"!
#include "Cocal/copy_coordinate_grav_extended_to_mpt.f90"!
#include "Cocal/copy_coordinate_grav_extended_from_mpt.f90"!
#include "Cocal/copy_coordinate_grav_r_to_mpt.f90"!
#include "Cocal/copy_coordinate_grav_r_from_mpt.f90"!
#include "Cocal/copy_trigonometry_grav_phi_to_mpt.f90"!
#include "Cocal/copy_trigonometry_grav_phi_from_mpt.f90"!
#include "Cocal/copy_def_binary_parameter_to_mpt.f90"!
#include "Cocal/copy_def_binary_parameter_from_mpt.f90"!
#include "Cocal/allocate_trigonometry_grav_phi_mpt.f90"!

#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Functions.h"
#include "cctk_Parameters.h"

!!!BNS ROUTINE NOW
!_____________________________________
 !
 !      MASTER CF BNS COCAL ID READER to CACTUS
 !______________________________________________
SUBROUTINE coc2cac_bns(CCTK_ARGUMENTS)
 !
   use grid_parameter_binary_excision
   use phys_constant
   use grid_parameter, eps_cocal => eps
   use coordinate_grav_r
   use coordinate_grav_extended
   use interface_modules_cartesian
 !  use interface_calc_gradvep_export
   use trigonometry_grav_phi
   use def_binary_parameter, only : dis
   use def_matter_parameter_mpt
   use interface_interpo_lag4th_2Dsurf
   use interface_IO_input_CF_grav_export
   use interface_IO_input_CF_flir_export !IR
   use interface_IO_input_CF_flco_export !CO
   use interface_IO_input_CF_flsp_export !SP
   use interface_IO_input_CF_surf_export
 !  use interface_excurve_CF_gridpoint_export
   use interface_IO_input_grav_export_Kij
   use interface_interpo_gr2fl_metric_CF_export
   use interface_IO_input_gradvep_export !unique to IR/SP
   
   
   

   implicit none
   DECLARE_CCTK_ARGUMENTS
   DECLARE_CCTK_FUNCTIONS
   DECLARE_CCTK_PARAMETERS
   integer :: i, j, k, nx, ny, nz
   integer :: impt, impt_ex, ico, irr, isp, igrid

   logical ::  bool_lapse, bool_shift, bool_hydro, bool_Bvec
   
    
   real(8) :: huta, alphfca2  !unique to IR, SP has this + more
   real(8) :: vepxfca, vepyfca, vepzfca      !Unique to IR/SP
   
   real(8) :: confpow, psifcacp, wxspfca, wyspfca, wzspfca, wdvep, w2, wterm,confpow_p1, confpow_p2 !unique to SP
   real(8), pointer :: wxspf_p1(:,:,:), wyspf_p1(:,:,:), wzspf_p1(:,:,:) !unique to SP
   real(8), pointer :: wxspf_p2(:,:,:), wyspf_p2(:,:,:), wzspf_p2(:,:,:) !unique to SP
   real(long) ::  fr4wxspf(4) , ft4wxspf(4) , fp4wxspf(4) !Unique to SP
   real(long) ::  fr4wyspf(4) , ft4wyspf(4) , fp4wyspf(4) !Unique to SP
   real(long) ::  fr4wzspf(4) , ft4wzspf(4) , fp4wzspf(4) !Unique to SP

   real(long) ::  fr4vepxf(4) , ft4vepxf(4) , fp4vepxf(4) !unique to IR/SP
   real(long) ::  fr4vepyf(4) , ft4vepyf(4) , fp4vepyf(4)
   real(long) ::  fr4vepzf(4) , ft4vepzf(4) , fp4vepzf(4)

   real(8), pointer :: vep_p1(:,:,:)                !unique to IR/SP
   real(8), pointer :: vepxf_p1(:,:,:), vepyf_p1(:,:,:), vepzf_p1(:,:,:) !unique to IR/SP, SP has more
   real(8), pointer :: vep_p2(:,:,:)                 !unique to IR/SP
   real(8), pointer :: vepxf_p2(:,:,:), vepyf_p2(:,:,:), vepzf_p2(:,:,:) !unique to IR/SP
   
   character(30) :: char1
   character(400) :: outstr
   character*400 :: dir_path, coc2cac_readformatf
   integer        :: dir_path_len, iAB, coc2cac_bnstypelen, coc2cac_readformatlen
   character(2)   :: id_type
   integer        :: ierr=0
   character(len=1) :: np(5) = (/'1', '2','3', '4', '5'/)
   real(8) :: rr3, rrcm, xc,yc,zc, xc_p1, yc_p1, zc_p1, xc_p2, yc_p2, zc_p2, dis_cm
   real(8) :: xc_p3, yc_p3, zc_p3
   real(8) :: xcac, ycac, zcac
   real(8) :: xcoc, ycoc, zcoc
   real(8) :: emdca, vepca, psica, alphca, bvxdca, bvydca, bvzdca, psi4ca, psif4ca
   real(8) :: hca, preca, rhoca, eneca, epsca
   real(8) :: vxu, vyu, vzu, lam_p1, lam_p2 
   real(8) :: bxcor, bycor, bzcor, bvxdfca, bvydfca, bvzdfca, psifca, alphfca
   real(8) :: gxx1, gxy1, gxz1, gyy1, gyz1, gzz1, kxx1, kxy1, kxz1, kyy1, kyz1, kzz1
   real(8) :: axx, axy, axz, ayy, ayz, azz
   real(8) :: ome_p1, ber_p1, radi_p1, r_surf_p1
   real(8) :: ome_p2, ber_p2, radi_p2, r_surf_p2
   integer :: nrg_p1,  ntg_p1, npg_p1, nrf_p1, ntf_p1, npf_p1
   integer :: nrg_p2,  ntg_p2, npg_p2, nrf_p2, ntf_p2, npf_p2
   integer :: nrg_p3,  ntg_p3, npg_p3, nrf_p3, ntf_p3, npf_p3
 !
   real(long) ::  rc_p1, thc_p1, phic_p1, varpic_p1, rcf_p1, rsca_p1
   real(long) ::  rc_p2, thc_p2, phic_p2, varpic_p2, rcf_p2, rsca_p2
   real(long) ::  rc_p3, thc_p3, phic_p3, varpic_p3
   real(long) ::  r4(4), th4(4), phi4(4), fr4(4), ft4(4), fp4(4)
   real(long) ::  fr4psi(4)   , ft4psi(4)   , fp4psi(4)
   real(long) ::  fr4alph(4)  , ft4alph(4)  , fp4alph(4)
   real(long) ::  fr4bvxd(4)  , ft4bvxd(4)  , fp4bvxd(4)
   real(long) ::  fr4bvyd(4)  , ft4bvyd(4)  , fp4bvyd(4)
   real(long) ::  fr4bvzd(4)  , ft4bvzd(4)  , fp4bvzd(4)
   real(long) ::  fr4axx(4)   , ft4axx(4)   , fp4axx(4)
   real(long) ::  fr4axy(4)   , ft4axy(4)   , fp4axy(4)
   real(long) ::  fr4axz(4)   , ft4axz(4)   , fp4axz(4)
   real(long) ::  fr4ayy(4)   , ft4ayy(4)   , fp4ayy(4)
   real(long) ::  fr4ayz(4)   , ft4ayz(4)   , fp4ayz(4)
   real(long) ::  fr4azz(4)   , ft4azz(4)   , fp4azz(4)
   real(long) ::  fr4emd(4)   , ft4emd(4)   , fp4emd(4)
 
   real(long) ::  fr4psif(4)  , ft4psif(4)  , fp4psif(4)
   real(long) ::  fr4alphf(4) , ft4alphf(4) , fp4alphf(4)
   real(long) ::  fr4bvxdf(4) , ft4bvxdf(4) , fp4bvxdf(4)
   real(long) ::  fr4bvydf(4) , ft4bvydf(4) , fp4bvydf(4)
   real(long) ::  fr4bvzdf(4) , ft4bvzdf(4) , fp4bvzdf(4)
 !
   integer :: irg, itg, ipg, irgex, itgex, ipgex
   integer :: ir0, it0, ip0, irg0 , itg0 , ipg0, ii, jj, kk
   real(long), external :: lagint_4th
 !
   real(8), pointer :: rg_p1(:), rgex_p1(:), thgex_p1(:), phigex_p1(:)
   integer, pointer :: irgex_r_p1(:), itgex_th_p1(:), ipgex_phi_p1(:)
   integer, pointer :: itgex_r_p1(:,:), ipgex_r_p1(:,:), ipgex_th_p1(:,:)
   real(8), pointer :: emd_p1(:,:,:), rs_p1(:,:)
   real(8), pointer :: psif_p1(:,:,:), alphf_p1(:,:,:), bvxdf_p1(:,:,:), bvydf_p1(:,:,:), bvzdf_p1(:,:,:)
   real(8), pointer :: psi_p1(:,:,:), alph_p1(:,:,:), bvxd_p1(:,:,:), bvyd_p1(:,:,:), bvzd_p1(:,:,:)
   real(8), pointer :: axx_p1(:,:,:), axy_p1(:,:,:) , axz_p1(:,:,:) , ayy_p1(:,:,:) , ayz_p1(:,:,:), azz_p1(:,:,:)
 !
   real(8), pointer :: rg_p2(:), rgex_p2(:), thgex_p2(:), phigex_p2(:)
   integer, pointer :: irgex_r_p2(:), itgex_th_p2(:), ipgex_phi_p2(:)
   integer, pointer :: itgex_r_p2(:,:), ipgex_r_p2(:,:), ipgex_th_p2(:,:)
   real(8), pointer :: emd_p2(:,:,:), rs_p2(:,:)
   real(8), pointer :: psif_p2(:,:,:), alphf_p2(:,:,:), bvxdf_p2(:,:,:), bvydf_p2(:,:,:), bvzdf_p2(:,:,:)
   real(8), pointer :: psi_p2(:,:,:), alph_p2(:,:,:), bvxd_p2(:,:,:), bvyd_p2(:,:,:), bvzd_p2(:,:,:)
   real(8), pointer :: axx_p2(:,:,:), axy_p2(:,:,:) , axz_p2(:,:,:) , ayy_p2(:,:,:) , ayz_p2(:,:,:), azz_p2(:,:,:)
 !
   real(8), pointer :: rg_p3(:), rgex_p3(:), thgex_p3(:), phigex_p3(:)
   integer, pointer :: irgex_r_p3(:), itgex_th_p3(:), ipgex_phi_p3(:)
   integer, pointer :: itgex_r_p3(:,:), ipgex_r_p3(:,:), ipgex_th_p3(:,:)
   real(8), pointer :: psi_p3(:,:,:), alph_p3(:,:,:), bvxd_p3(:,:,:), bvyd_p3(:,:,:), bvzd_p3(:,:,:)
   real(8), pointer :: axx_p3(:,:,:), axy_p3(:,:,:) , axz_p3(:,:,:) , ayy_p3(:,:,:) , ayz_p3(:,:,:), azz_p3(:,:,:)
 !
   
 
 
   call CCTK_FortranString(dir_path_len,coc2cac_dir_path_ID,dir_path)  !For reading initial data
 ! -- Read ID type
   call CCTK_INFO("COCAL_ID:: Reading BNS Type")
   
   ierr = read_id_type(trim(dir_path)//"/rnspar_mpt1.dat",id_type)

   write(outstr,'(a25,a2)') "COCAL_ID:: BNS Type: ", id_type
   call CCTK_INFO(outstr)
   call CCTK_INFO(trim(id_type))

   !
   if (ierr.ne.0)     call CCTK_INFO("COCAL_ID:: problem reading file rnspar_mpt1.dat.")
 
   select case (trim(id_type))
     case ("CO")
       call CCTK_INFO("COCAL_ID:: Reading corotating BNS ID")
     case ("IR")
       call CCTK_INFO("COCAL_ID:: Reading irrotational BNS ID")
     case ("SP")
       call CCTK_INFO("COCAL_ID:: Reading spinning BNS ID")
   end select
 !
 ! The coc2cac_ini_sub.F90 subroutine has here option for zero initial shift
 !
   gxx1=0.0d0; gxy1=0.0d0; gxz1=0.0d0; gyy1=0.0d0; gyz1=0.0d0; gzz1=0.0d0
   kxx1=0.0d0; kxy1=0.0d0; kxz1=0.0d0; kyy1=0.0d0; kyz1=0.0d0; kzz1=0.0d0
   axx=0.0d0; axy=0.0d0; axz=0.0d0; ayy=0.0d0; ayz=0.0d0; azz=0.0d0
 
   !TODO remove this
   !dir_path="/home/astro/mundim/tmp/ET_2014_05_wheeler/Cactus/repos/Cocal/standalone/Cocal/ID_BNS"
   !dir_path="../../standalone/Cocal/ID_BNS"
   !dir_path='.'
   call CCTK_INFO("In coc2cac_ir: dir_path="//dir_path)
   call CCTK_FortranString(coc2cac_readformatlen,coc2cac_readformat,coc2cac_readformatf)
 
 !--------------------- Choose gravitational grid -----------------------
 !igrid3  igrid = 3     ! 3:r_surf is used
   igrid = 4     ! 4:r_surf=1.0
 !-----------------------------------------------------------------------
 
   call CCTK_INFO("Cocal: Reading parameters")
 ! -- Read parameters
   call allocate_grid_parameter_mpt
   call allocate_grid_parameter_binary_excision_mpt
   call allocate_def_matter_parameter_mpt
   do impt = 1, nmpt
     call CCTK_INFO("read_parameter_mpt_cactus")
     call read_parameter_mpt_cactus(impt,dir_path)
     indata_type = '3D'

     call CCTK_INFO("read_surf_parameter_mpt_cactus")
     if (impt .le. 2) call read_surf_parameter_mpt_cactus(impt,dir_path)
     call copy_grid_parameter_to_mpt(impt)

     call CCTK_INFO("read_parameter_binary_excision_mpt_cactus")
     call read_parameter_binary_excision_mpt_cactus(impt,dir_path)
     call copy_grid_parameter_binary_excision_to_mpt(impt)
     
     call CCTK_INFO("peos_initialize_mpt_cactus")
     if (impt .le. 2) call peos_initialize_mpt_cactus(impt,dir_path)
     call copy_def_peos_parameter_to_mpt(impt)
   end do
 ! -- Allocate arrays
   call CCTK_INFO("Cocal: Allocating local arrays...")
 
   call set_allocate_size_mpt
   call allocate_trig_grav_mphi
   call allocate_trigonometry_grav_phi_mpt
 !
   do impt = 1, nmpt
     call copy_grid_parameter_from_mpt(impt)
     call copy_grid_parameter_binary_excision_from_mpt(impt)
     call copy_def_peos_parameter_from_mpt(impt)
     call coordinate_patch_kit_grav_grid_coc2cac_mpt(igrid)  ! 3:r_surf is used
     call calc_parameter_binary_excision
     call copy_grid_parameter_to_mpt(impt)
     call copy_grid_parameter_binary_excision_to_mpt(impt)
     call copy_coordinate_grav_extended_to_mpt(impt)
     call copy_coordinate_grav_phi_to_mpt(impt)
     call copy_coordinate_grav_r_to_mpt(impt)
     call copy_coordinate_grav_theta_to_mpt(impt)
     call copy_def_binary_parameter_to_mpt(impt)
     call copy_trigonometry_grav_phi_to_mpt(impt)
     call copy_trigonometry_grav_theta_to_mpt(impt)
   end do
 
   call CCTK_INFO("Cocal: Internal reading info (BEGIN):")
 
   do impt = 1,nmpt
 !=>    call copy_from_mpatch_interpolation_utility(impt)
     call copy_grid_parameter_from_mpt(impt)
     call copy_grid_parameter_binary_excision_from_mpt(impt)
     call copy_coordinate_grav_extended_from_mpt(impt)
     call copy_coordinate_grav_phi_from_mpt(impt)
     call copy_coordinate_grav_r_from_mpt(impt)
     call copy_coordinate_grav_theta_from_mpt(impt)
     call copy_def_binary_parameter_from_mpt(impt)
     call copy_trigonometry_grav_phi_from_mpt(impt)
     call copy_trigonometry_grav_theta_from_mpt(impt)
 
 !    write(6,'(6i5)') nrg, ntg, npg, nrf, ntf, npf
     if (impt==1) then
       nrg_p1=nrg;  ntg_p1=ntg;  npg_p1=npg;  nrf_p1=nrf;  ntf_p1=ntf;  npf_p1=npf
       call CCTK_INFO("Allocating star1...")
       allocate (       rg_p1( 0:nnrg))
       allocate (     rgex_p1(-2:nnrg+2))
       allocate (    thgex_p1(-2:nntg+2))
       allocate (   phigex_p1(-2:nnpg+2))
       allocate (  irgex_r_p1(-2:nnrg+2))
       allocate ( itgex_th_p1(-2:nntg+2))
       allocate (ipgex_phi_p1(-2:nnpg+2))
 
       allocate ( itgex_r_p1(0:nntg,-2:nnrg+2))
       allocate ( ipgex_r_p1(0:nnpg,-2:nnrg+2))
       allocate (ipgex_th_p1(0:nnpg,-2:nntg+2))
 
       rg_p1=rg;    rgex_p1=rgex;    thgex_p1=thgex;    phigex_p1=phigex
       irgex_r_p1=irgex_r;  itgex_th_p1=itgex_th;  ipgex_phi_p1=ipgex_phi
       itgex_r_p1=itgex_r;  ipgex_r_p1=ipgex_r;    ipgex_th_p1=ipgex_th
 
       rr3 = 0.7d0*(rgout - rgmid)
       dis_cm = dis
       r_surf_p1 = r_surf
       
 
       allocate (  emd_p1(0:nrf,0:ntf,0:npf))
       select case (trim(id_type))
          case ("IR","SP")
             allocate (  vep_p1(0:nrf,0:ntf,0:npf)) ! Unique to IR/SP, SP has 3 more
             allocate (vepxf_p1(0:nrf,0:ntf,0:npf)) ! Unique to IR/SP, SP has 3 more
             allocate (vepyf_p1(0:nrf,0:ntf,0:npf)) ! Unique to IR/SP, SP has 3 more
             allocate (vepzf_p1(0:nrf,0:ntf,0:npf)) ! Unique to IR/SP, SP has 3 more
             vep_p1 =0.0d0;  vepxf_p1=0.0d0;  vepyf_p1=0.0d0;  vepzf_p1=0.0d0 ! Unique to IR/SP, SP has 3 more
             select case (trim(id_type))
                case ("SP")
                   allocate (wxspf_p1(0:nrf,0:ntf,0:npf)) ! Unique to SP
                   allocate (wyspf_p1(0:nrf,0:ntf,0:npf)) ! Unique to SP
                   allocate (wzspf_p1(0:nrf,0:ntf,0:npf)) ! Unique to SP
                   wxspf_p1=0.0d0; wyspf_p1=0.0d0; wzspf_p1=0.0d0 ! Unique to SP
             end select
       end select
       allocate ( psif_p1(0:nrf,0:ntf,0:npf))
       allocate (alphf_p1(0:nrf,0:ntf,0:npf))
       allocate (bvxdf_p1(0:nrf,0:ntf,0:npf))
       allocate (bvydf_p1(0:nrf,0:ntf,0:npf))
       allocate (bvzdf_p1(0:nrf,0:ntf,0:npf))
       allocate (   rs_p1(0:ntf,0:npf))
       allocate (  psi_p1(0:nrg,0:ntg,0:npg))
       allocate ( alph_p1(0:nrg,0:ntg,0:npg))
       allocate ( bvxd_p1(0:nrg,0:ntg,0:npg))
       allocate ( bvyd_p1(0:nrg,0:ntg,0:npg))
       allocate ( bvzd_p1(0:nrg,0:ntg,0:npg))
       allocate (  axx_p1(0:nrg,0:ntg,0:npg))
       allocate (  axy_p1(0:nrg,0:ntg,0:npg))
       allocate (  axz_p1(0:nrg,0:ntg,0:npg))
       allocate (  ayy_p1(0:nrg,0:ntg,0:npg))
       allocate (  ayz_p1(0:nrg,0:ntg,0:npg))
       allocate (  azz_p1(0:nrg,0:ntg,0:npg))
       call CCTK_INFO("Done allocating star 1...")
 
       write(outstr,'(a13,i1,a8)') "Reading COCP-", impt, " data..."
       call CCTK_INFO(outstr)

       
       emd_p1=0.0d0;  rs_p1  =0.0d0;  
       psi_p1=0.0d0;  alph_p1=0.0d0;  bvxd_p1=0.0d0;  bvyd_p1=0.0d0;  bvzd_p1=0.0d0
       axx_p1=0.0d0;  axy_p1 =0.0d0;  axz_p1 =0.0d0;   ayy_p1=0.0d0;   ayz_p1=0.0d0;   azz_p1=0.0d0
 
       call IO_input_CF_grav_export(trim(dir_path)//"/bnsgra_3D_mpt1.las",coc2cac_readformatf,psi_p1,alph_p1,bvxd_p1,bvyd_p1,bvzd_p1)
       select case (trim(id_type))
          case("CO")
             call IO_input_CF_flco_export(trim(dir_path)//"/bnsflu_3D_mpt1.las",coc2cac_readformatf,emd_p1,ome_p1,ber_p1,radi_p1)
          case("IR")
             call IO_input_CF_flir_export(trim(dir_path)//"/bnsflu_3D_mpt1.las",coc2cac_readformatf,emd_p1,vep_p1,ome_p1,ber_p1,radi_p1) ! This line changes for IR/CO/SP
          case("SP")
             call IO_input_CF_flsp_export(trim(dir_path)//"/bnsflu_3D_mpt1.las",coc2cac_readformatf,emd_p1,vep_p1,wxspf_p1,wyspf_p1,wzspf_p1,ome_p1,ber_p1,radi_p1,confpow_p1)
       end select
            
 
       call IO_input_CF_surf_export(trim(dir_path)//"/bnssur_3D_mpt1.las",coc2cac_readformatf,rs_p1)
 
 !      call excurve_CF_gridpoint_export(alph_p1,bvxd_p1,bvyd_p1,bvzd_p1, axx_p1, axy_p1, axz_p1, ayy_p1, ayz_p1, azz_p1)
 
       call IO_input_grav_export_Kij(trim(dir_path)//"/bnsgra_Kij_3D_mpt1.las", axx_p1, axy_p1, axz_p1, ayy_p1, ayz_p1, azz_p1)
 
       call interpo_gr2fl_metric_CF_export(alph_p1, psi_p1, bvxd_p1, bvyd_p1, bvzd_p1, &
          &    alphf_p1, psif_p1, bvxdf_p1, bvydf_p1, bvzdf_p1, rs_p1)
 
 !      call calc_gradvep_export(vep_p1,vepxf_p1,vepyf_p1,vepzf_p1,rs_p1)

       select case (trim(id_type))
          case ("IR","SP")
             call IO_input_gradvep_export(trim(dir_path)//"/bnsdvep_3D_mpt1.las",coc2cac_readformatf,  vepxf_p1, vepyf_p1, vepzf_p1) ! Unique to IR/SP
       end select
 
     end if
     if (impt==2) then
       call CCTK_INFO("Allocating star2...")
       nrg_p2=nrg;  ntg_p2=ntg;  npg_p2=npg;  nrf_p2=nrf;  ntf_p2=ntf;  npf_p2=npf
       allocate (       rg_p2( 0:nnrg))
       allocate (     rgex_p2(-2:nnrg+2))
       allocate (    thgex_p2(-2:nntg+2))
       allocate (   phigex_p2(-2:nnpg+2))
       allocate (  irgex_r_p2(-2:nnrg+2))
       allocate ( itgex_th_p2(-2:nntg+2))
       allocate (ipgex_phi_p2(-2:nnpg+2))
 
       allocate ( itgex_r_p2(0:nntg,-2:nnrg+2))
       allocate ( ipgex_r_p2(0:nnpg,-2:nnrg+2))
       allocate (ipgex_th_p2(0:nnpg,-2:nntg+2))
 
       rg_p2=rg;    rgex_p2=rgex;    thgex_p2=thgex;    phigex_p2=phigex
       irgex_r_p2=irgex_r;  itgex_th_p2=itgex_th;  ipgex_phi_p2=ipgex_phi
       itgex_r_p2=itgex_r;  ipgex_r_p2=ipgex_r;    ipgex_th_p2=ipgex_th
 
       r_surf_p2 = r_surf
       
 
       allocate (  emd_p2(0:nrf,0:ntf,0:npf))
       select case (trim(id_type))
          case ("IR","SP")
             allocate (  vep_p2(0:nrf,0:ntf,0:npf)) ! Unique to IR/SP, SP has 3 more
             allocate (vepxf_p2(0:nrf,0:ntf,0:npf)) ! Unique to IR/SP, SP has 3 more
             allocate (vepyf_p2(0:nrf,0:ntf,0:npf)) ! Unique to IR/SP, SP has 3 more
             allocate (vepzf_p2(0:nrf,0:ntf,0:npf)) ! Unique to IR/SP, SP has 3 more
             vep_p2 =0.0d0;  vepxf_p2=0.0d0;  vepyf_p2=0.0d0;  vepzf_p2=0.0d0 ! Unique to IR/SP, SP has 3 more
             select case (trim(id_type))
                case ("SP")
                   allocate (wxspf_p2(0:nrf,0:ntf,0:npf)) ! Unique to SP
                   allocate (wyspf_p2(0:nrf,0:ntf,0:npf)) ! Unique to SP
                   allocate (wzspf_p2(0:nrf,0:ntf,0:npf)) ! Unique to SP
                   wxspf_p2=0.0d0; wyspf_p2=0.0d0; wzspf_p2=0.0d0 ! Unique to SP
             end select
       end select
       allocate ( psif_p2(0:nrf,0:ntf,0:npf))
       allocate (alphf_p2(0:nrf,0:ntf,0:npf))
       allocate (bvxdf_p2(0:nrf,0:ntf,0:npf))
       allocate (bvydf_p2(0:nrf,0:ntf,0:npf))
       allocate (bvzdf_p2(0:nrf,0:ntf,0:npf))
       allocate (   rs_p2(0:ntf,0:npf))
       allocate (  psi_p2(0:nrg,0:ntg,0:npg))
       allocate ( alph_p2(0:nrg,0:ntg,0:npg))
       allocate ( bvxd_p2(0:nrg,0:ntg,0:npg))
       allocate ( bvyd_p2(0:nrg,0:ntg,0:npg))
       allocate ( bvzd_p2(0:nrg,0:ntg,0:npg))
       allocate (  axx_p2(0:nrg,0:ntg,0:npg))
       allocate (  axy_p2(0:nrg,0:ntg,0:npg))
       allocate (  axz_p2(0:nrg,0:ntg,0:npg))
       allocate (  ayy_p2(0:nrg,0:ntg,0:npg))
       allocate (  ayz_p2(0:nrg,0:ntg,0:npg))
       allocate (  azz_p2(0:nrg,0:ntg,0:npg))
       call CCTK_INFO("Done allocating star 2...")
       write(outstr,'(a13,i1,a8)') "Reading COCP-", impt, " data..."
       call CCTK_INFO(outstr)

       emd_p2=0.0d0;  rs_p2=0.0d0
       psi_p2=0.0d0;  alph_p2=0.0d0;  bvxd_p2=0.0d0;   bvyd_p2=0.0d0;  bvzd_p2=0.0d0
       axx_p2=0.0d0;   axy_p2=0.0d0;   axz_p2=0.0d0;    ayy_p2=0.0d0;   ayz_p2=0.0d0;   azz_p2=0.0d0
 
       call IO_input_CF_grav_export(trim(dir_path)//"/bnsgra_3D_mpt2.las",coc2cac_readformatf,psi_p2,alph_p2,bvxd_p2,bvyd_p2,bvzd_p2)
 
       select case (trim(id_type))
       case("CO")
          call IO_input_CF_flco_export(trim(dir_path)//"/bnsflu_3D_mpt2.las",coc2cac_readformatf,emd_p2,ome_p2,ber_p2,radi_p2)
       case("IR")
          call IO_input_CF_flir_export(trim(dir_path)//"/bnsflu_3D_mpt2.las",coc2cac_readformatf,emd_p2,vep_p2,ome_p2,ber_p2,radi_p2) ! This line changes for IR/CO/SP
       case("SP")
          call IO_input_CF_flsp_export(trim(dir_path)//"/bnsflu_3D_mpt2.las",coc2cac_readformatf,emd_p2,vep_p2,wxspf_p2,wyspf_p2,wzspf_p2,ome_p2,ber_p2,radi_p2,confpow_p2)
       end select
 
       call IO_input_CF_surf_export(trim(dir_path)//"/bnssur_3D_mpt2.las",coc2cac_readformatf,rs_p2)
 
 !      call excurve_CF_gridpoint_export(alph_p2,bvxd_p2,bvyd_p2,bvzd_p2,  axx_p2, axy_p2, axz_p2, ayy_p2, ayz_p2, azz_p2)
 
       call IO_input_grav_export_Kij(trim(dir_path)//"/bnsgra_Kij_3D_mpt2.las", axx_p2, axy_p2, axz_p2, ayy_p2, ayz_p2, azz_p2)
 
       call interpo_gr2fl_metric_CF_export(alph_p2, psi_p2, bvxd_p2, bvyd_p2, bvzd_p2, &
          &    alphf_p2, psif_p2, bvxdf_p2, bvydf_p2, bvzdf_p2, rs_p2)
 
 !      call calc_gradvep_export(vep_p2,vepxf_p2,vepyf_p2,vepzf_p2,rs_p2)
 
       select case (trim(id_type))
          case ("IR","SP")
             call IO_input_gradvep_export(trim(dir_path)//"/bnsdvep_3D_mpt2.las",coc2cac_readformatf,  vepxf_p2, vepyf_p2, vepzf_p2) ! Unique to IR/SP
             vepxf_p2 = - vepxf_p2 ! Unique to IR/SP
             vepyf_p2 = - vepyf_p2 ! Unique to IR/SP
            select case (trim(id_type))
               case ("SP")
                  wxspf_p2 = - wxspf_p2 !unique to SP
                  wyspf_p2 = - wyspf_p2 !unique to SP
            end select
       end select
       bvxdf_p2 = - bvxdf_p2
       bvydf_p2 = - bvydf_p2
       bvxd_p2  = - bvxd_p2
       bvyd_p2  = - bvyd_p2
       axz_p2   = - axz_p2
       ayz_p2   = - ayz_p2
     end if
     if (impt==3) then
       nrg_p3=nrg;  ntg_p3=ntg;  npg_p3=npg;  nrf_p3=nrf;  ntf_p3=ntf;  npf_p3=npf
       allocate (       rg_p3( 0:nnrg))
       allocate (     rgex_p3(-2:nnrg+2))
       allocate (    thgex_p3(-2:nntg+2))
       allocate (   phigex_p3(-2:nnpg+2))
       allocate (  irgex_r_p3(-2:nnrg+2))
       allocate ( itgex_th_p3(-2:nntg+2))
       allocate (ipgex_phi_p3(-2:nnpg+2))
 
       allocate ( itgex_r_p3(0:nntg,-2:nnrg+2))
       allocate ( ipgex_r_p3(0:nnpg,-2:nnrg+2))
       allocate (ipgex_th_p3(0:nnpg,-2:nntg+2))
 
       rg_p3=rg;    rgex_p3=rgex;    thgex_p3=thgex;    phigex_p3=phigex
       irgex_r_p3=irgex_r;  itgex_th_p3=itgex_th;  ipgex_phi_p3=ipgex_phi
       itgex_r_p3=itgex_r;  ipgex_r_p3=ipgex_r;    ipgex_th_p3=ipgex_th
 
       write(6,'(a13,i1,a8)') "Reading ARCP-", impt, " data..."
       allocate ( psi_p3(0:nrg,0:ntg,0:npg))
       allocate (alph_p3(0:nrg,0:ntg,0:npg))
       allocate (bvxd_p3(0:nrg,0:ntg,0:npg))
       allocate (bvyd_p3(0:nrg,0:ntg,0:npg))
       allocate (bvzd_p3(0:nrg,0:ntg,0:npg))
       allocate ( axx_p3(0:nrg,0:ntg,0:npg))
       allocate ( axy_p3(0:nrg,0:ntg,0:npg))
       allocate ( axz_p3(0:nrg,0:ntg,0:npg))
       allocate ( ayy_p3(0:nrg,0:ntg,0:npg))
       allocate ( ayz_p3(0:nrg,0:ntg,0:npg))
       allocate ( azz_p3(0:nrg,0:ntg,0:npg))
       psi_p3=0.0d0;  alph_p3=0.0d0;  bvxd_p3=0.0d0;  bvyd_p3=0.0d0; bvzd_p3=0.0d0
       axx_p3=0.0d0;   axy_p3=0.0d0;   axz_p3=0.0d0;   ayy_p3=0.0d0;  ayz_p3=0.0d0; azz_p3=0.0d0
 
       call IO_input_CF_grav_export(trim(dir_path)//"/bnsgra_3D_mpt3.las",coc2cac_readformatf,psi_p3,alph_p3,bvxd_p3,bvyd_p3,bvzd_p3)
       
       call IO_input_grav_export_Kij(trim(dir_path)//"/bnsgra_Kij_3D_mpt3.las", axx_p3, axy_p3, axz_p3, ayy_p3, ayz_p3, azz_p3)
 
 !      call excurve_CF_gridpoint_export(alph_p3,bvxd_p3,bvyd_p3,bvzd_p3, axx_p3, axy_p3, axz_p3, ayy_p3, ayz_p3, azz_p3)
 
     end if
   end do
   write(6,'(2e20.12)') emd_p1(0,0,0), emd_p1(58,0,0)
   write(6,'(3e20.12)') ome_p1, ber_p1, radi_p1
   write(6,'(e20.12)') dis_cm
 !
   call CCTK_INFO("Cocal: Internal reading info (END).")
   call CCTK_INFO("Cocal: Looping over local cartesian grid:")
 
   nx = cctk_lsh(1)
   ny = cctk_lsh(2)
   nz = cctk_lsh(3)
   if (CCTK_EQUALS(initial_hydro, "Cocal")) then
      bool_hydro = .true.
    else
      bool_hydro = .false.
    end if
    if (CCTK_EQUALS(initial_lapse, "Cocal")) then
       bool_lapse = .true.
    else
       bool_lapse = .false.
    end if
    if (CCTK_EQUALS(initial_shift, "Cocal")) then
       bool_shift = .true.
    else
       bool_shift = .false.
    end if
    if (CCTK_EQUALS(initial_Bvec, "Cocal")) then
       bool_Bvec = .true.
    else
       bool_Bvec = .false.
    end if
 
   ! write(6,'(a56)', ADVANCE = "NO") "Give cartesian coordinates (x,y,z) separated by a space:"
   ! read(5,*) xcac,ycac,zcac
   ! write(6,'(a23,3e20.12)') "Point given wrt CACTUS:", xcac,ycac,zcac
   ! write(6,'(a38,2e20.12)') "Cocal radius scale in COCP-1, COCP-2 :", radi_p1, radi_p2
   ! write(6,'(a38,2e20.12)') "Cocal surface scale in COCP-1, COCP-2:", r_surf_p1, r_surf_p2
   do k = 1, nz
     do j = 1, ny
       do i = 1, nx
         xcac = x(i,j,k)
         ycac = y(i,j,k)
         zcac = z(i,j,k)
         !write(6,*)' i, j, k, xcac, ycac, zcac', i, j, k, xcac, ycac, zcac
         
         xcoc = xcac/(radi_p1)
         ycoc = ycac/(radi_p1)
         zcoc = zcac/(radi_p1)
         ! write(6,'(a23,3e20.12)') "Point given wrt COCAL:", xcoc,ycoc,zcoc     
 
         rrcm = sqrt(xcoc**2 + ycoc**2 + zcoc**2)  
         !write(6,*)  rrcm, rr3
         if (rrcm > rr3) then
       !=>    call copy_from_mpatch_interpolation_utility(3)
       !    call copy_grid_parameter_from_mpt(3)
       !    call copy_grid_parameter_binary_excision_from_mpt(3)
       !    call copy_coordinate_grav_extended_from_mpt(3)
       !    call copy_coordinate_grav_phi_from_mpt(3)
       !    call copy_coordinate_grav_r_from_mpt(3)
       !    call copy_coordinate_grav_theta_from_mpt(3)
       !    call copy_def_binary_parameter_from_mpt(3)
       !    call copy_trigonometry_grav_phi_from_mpt(3)
       !    call copy_trigonometry_grav_theta_from_mpt(3)
 
           xc_p3 = xcoc
           yc_p3 = ycoc
           zc_p3 = zcoc
         !   write(outstr,'(a39,3e20.12)') "Point given wrt COCP-3 in Cocal scale :", xc_p3,yc_p3,zc_p3
         !   call CCTK_INFO(outstr)
 
           psica=0.0d0; alphca=0.0d0; bvxdca=0.0d0; bvydca=0.0d0; bvzdca=0.0d0
           axx=0.0d0  ; axy=0.0d0   ; axz=0.0d0   ; ayy=0.0d0   ; ayz=0.0d0   ; azz=0.0d0
       !
           rc_p3     = dsqrt(dabs(xc_p3**2 + yc_p3**2 + zc_p3**2))
           varpic_p3 = dsqrt(dabs(xc_p3**2 + yc_p3**2))
           thc_p3  = dmod(2.0d0*pi + datan2(varpic_p3,zc_p3),2.0d0*pi)
           phic_p3 = dmod(2.0d0*pi + datan2(    yc_p3,xc_p3),2.0d0*pi)
       !
           do irg = 0, nrg_p3+1
             if (rc_p3.lt.rgex_p3(irg).and.rc_p3.ge.rgex_p3(irg-1)) ir0 = min0(irg-2,nrg_p3-3)
           end do
           do itg = 0, ntg_p3+1
             if (thc_p3.lt.thgex_p3(itg).and.thc_p3.ge.thgex_p3(itg-1)) it0 = itg-2
           end do
           do ipg = 0, npg_p3+1
             if (phic_p3.lt.phigex_p3(ipg).and.phic_p3.ge.phigex_p3(ipg-1)) ip0 = ipg-2
           end do
       !
           do ii = 1, 4
             irg0 = ir0 + ii - 1
             itg0 = it0 + ii - 1
             ipg0 = ip0 + ii - 1
             r4(ii) = rgex_p3(irg0)
             th4(ii) = thgex_p3(itg0)
             phi4(ii) = phigex_p3(ipg0)
           end do
       !
           do kk = 1, 4
             ipg0 = ip0 + kk - 1
             do jj = 1, 4
               itg0 = it0 + jj - 1
               do ii = 1, 4
                 irg0 = ir0 + ii - 1
                 irgex = irgex_r_p3(irg0)
                 itgex = itgex_r_p3(itgex_th_p3(itg0),irg0)
                 ipgex = ipgex_r_p3(ipgex_th_p3(ipgex_phi_p3(ipg0),itg0),irg0)
                 fr4psi(ii)  =  psi_p3(irgex,itgex,ipgex)
                 fr4alph(ii) = alph_p3(irgex,itgex,ipgex)
                 fr4bvxd(ii) = bvxd_p3(irgex,itgex,ipgex)
                 fr4bvyd(ii) = bvyd_p3(irgex,itgex,ipgex)
                 fr4bvzd(ii) = bvzd_p3(irgex,itgex,ipgex)
                 fr4axx(ii)  =  axx_p3(irgex,itgex,ipgex)
                 fr4axy(ii)  =  axy_p3(irgex,itgex,ipgex)
                 fr4axz(ii)  =  axz_p3(irgex,itgex,ipgex)
                 fr4ayy(ii)  =  ayy_p3(irgex,itgex,ipgex)
                 fr4ayz(ii)  =  ayz_p3(irgex,itgex,ipgex)
                 fr4azz(ii)  =  azz_p3(irgex,itgex,ipgex)
               end do
               ft4psi(jj)  = lagint_4th(r4,fr4psi ,rc_p3)
               ft4alph(jj) = lagint_4th(r4,fr4alph,rc_p3)
               ft4bvxd(jj) = lagint_4th(r4,fr4bvxd,rc_p3)
               ft4bvyd(jj) = lagint_4th(r4,fr4bvyd,rc_p3)
               ft4bvzd(jj) = lagint_4th(r4,fr4bvzd,rc_p3)
               ft4axx(jj)  = lagint_4th(r4,fr4axx ,rc_p3)
               ft4axy(jj)  = lagint_4th(r4,fr4axy ,rc_p3)
               ft4axz(jj)  = lagint_4th(r4,fr4axz ,rc_p3)
               ft4ayy(jj)  = lagint_4th(r4,fr4ayy ,rc_p3)
               ft4ayz(jj)  = lagint_4th(r4,fr4ayz ,rc_p3)
               ft4azz(jj)  = lagint_4th(r4,fr4azz ,rc_p3)
             end do
             fp4psi(kk)  = lagint_4th(th4,ft4psi ,thc_p3)
             fp4alph(kk) = lagint_4th(th4,ft4alph,thc_p3)
             fp4bvxd(kk) = lagint_4th(th4,ft4bvxd,thc_p3)
             fp4bvyd(kk) = lagint_4th(th4,ft4bvyd,thc_p3)
             fp4bvzd(kk) = lagint_4th(th4,ft4bvzd,thc_p3)
             fp4axx(kk)  = lagint_4th(th4,ft4axx ,thc_p3)
             fp4axy(kk)  = lagint_4th(th4,ft4axy ,thc_p3)
             fp4axz(kk)  = lagint_4th(th4,ft4axz ,thc_p3)
             fp4ayy(kk)  = lagint_4th(th4,ft4ayy ,thc_p3)
             fp4ayz(kk)  = lagint_4th(th4,ft4ayz ,thc_p3)
             fp4azz(kk)  = lagint_4th(th4,ft4azz ,thc_p3)
           end do
           psica  = lagint_4th(phi4,fp4psi ,phic_p3)
           alphca = lagint_4th(phi4,fp4alph,phic_p3)
           bvxdca = lagint_4th(phi4,fp4bvxd,phic_p3)
           bvydca = lagint_4th(phi4,fp4bvyd,phic_p3)
           bvzdca = lagint_4th(phi4,fp4bvzd,phic_p3)
           axx    = lagint_4th(phi4,fp4axx ,phic_p3)
           axy    = lagint_4th(phi4,fp4axy ,phic_p3)
           axz    = lagint_4th(phi4,fp4axz ,phic_p3)
           ayy    = lagint_4th(phi4,fp4ayy ,phic_p3)
           ayz    = lagint_4th(phi4,fp4ayz ,phic_p3)
           azz    = lagint_4th(phi4,fp4azz ,phic_p3)
 
           psi4ca = psica**4
           vxu = 0.0d0
           vyu = 0.0d0
           vzu = 0.0d0
           emdca = 0.0d0
         else
           if (xcoc<=0.0d0) then
       !=>      call copy_from_mpatch_interpolation_utility(1)
       !      call copy_grid_parameter_from_mpt(1)
       !      call copy_grid_parameter_binary_excision_from_mpt(1)
       !      call copy_coordinate_grav_extended_from_mpt(1)
       !      call copy_coordinate_grav_phi_from_mpt(1)
       !      call copy_coordinate_grav_r_from_mpt(1)
       !      call copy_coordinate_grav_theta_from_mpt(1)
       !      call copy_def_binary_parameter_from_mpt(1)
       !      call copy_trigonometry_grav_phi_from_mpt(1)
       !      call copy_trigonometry_grav_theta_from_mpt(1)
 
             xc_p1 = xcoc + dis_cm
             yc_p1 = ycoc
             zc_p1 = zcoc
            !  write(6,'(a39,3e20.12)') "Point given wrt COCP-1 in Cocal scale :", xc_p1,yc_p1,zc_p1  
       !
             psica=0.0d0; alphca=0.0d0; bvxdca=0.0d0; bvydca=0.0d0; bvzdca=0.0d0
             axx=0.0d0  ; axy=0.0d0   ; axz=0.0d0   ; ayy=0.0d0   ; ayz=0.0d0   ; azz=0.0d0
       !
             rc_p1     = dsqrt(dabs(xc_p1**2 + yc_p1**2 + zc_p1**2))
             varpic_p1 = dsqrt(dabs(xc_p1**2 + yc_p1**2))
             thc_p1  = dmod(2.0d0*pi + datan2(varpic_p1,zc_p1),2.0d0*pi)
             phic_p1 = dmod(2.0d0*pi + datan2(    yc_p1,xc_p1),2.0d0*pi)
       !
             do irg = 0, nrg_p1+1
               if (rc_p1.lt.rgex_p1(irg).and.rc_p1.ge.rgex_p1(irg-1)) ir0 = min0(irg-2,nrg_p1-3)
             end do
             do itg = 0, ntg_p1+1
               if (thc_p1.lt.thgex_p1(itg).and.thc_p1.ge.thgex_p1(itg-1)) it0 = itg-2
             end do
             do ipg = 0, npg_p1+1
               if (phic_p1.lt.phigex_p1(ipg).and.phic_p1.ge.phigex_p1(ipg-1)) ip0 = ipg-2
             end do
       !
       !      write(6,*) "(ir0,it0,ip0)=", ir0,it0,ip0
 
             do ii = 1, 4
               irg0 = ir0 + ii - 1
               itg0 = it0 + ii - 1
               ipg0 = ip0 + ii - 1
               r4(ii) = rgex_p1(irg0)
               th4(ii) = thgex_p1(itg0)
               phi4(ii) = phigex_p1(ipg0)
             end do
       !
             do kk = 1, 4
               ipg0 = ip0 + kk - 1
               do jj = 1, 4
                 itg0 = it0 + jj - 1
                 do ii = 1, 4
                   irg0 = ir0 + ii - 1
                   irgex = irgex_r_p1(irg0)
                   itgex = itgex_r_p1(itgex_th_p1(itg0),irg0)
                   ipgex = ipgex_r_p1(ipgex_th_p1(ipgex_phi_p1(ipg0),itg0),irg0)
                   fr4psi(ii)  =  psi_p1(irgex,itgex,ipgex)
                   fr4alph(ii) = alph_p1(irgex,itgex,ipgex)
                   fr4bvxd(ii) = bvxd_p1(irgex,itgex,ipgex)
                   fr4bvyd(ii) = bvyd_p1(irgex,itgex,ipgex)
                   fr4bvzd(ii) = bvzd_p1(irgex,itgex,ipgex)
                   fr4axx(ii)  =  axx_p1(irgex,itgex,ipgex)
                   fr4axy(ii)  =  axy_p1(irgex,itgex,ipgex)
                   fr4axz(ii)  =  axz_p1(irgex,itgex,ipgex)
                   fr4ayy(ii)  =  ayy_p1(irgex,itgex,ipgex)
                   fr4ayz(ii)  =  ayz_p1(irgex,itgex,ipgex)
                   fr4azz(ii)  =  azz_p1(irgex,itgex,ipgex)
                 end do
                 ft4psi(jj)  = lagint_4th(r4,fr4psi ,rc_p1)
                 ft4alph(jj) = lagint_4th(r4,fr4alph,rc_p1)
                 ft4bvxd(jj) = lagint_4th(r4,fr4bvxd,rc_p1)
                 ft4bvyd(jj) = lagint_4th(r4,fr4bvyd,rc_p1)
                 ft4bvzd(jj) = lagint_4th(r4,fr4bvzd,rc_p1)
                 ft4axx(jj)  = lagint_4th(r4,fr4axx ,rc_p1)
                 ft4axy(jj)  = lagint_4th(r4,fr4axy ,rc_p1)
                 ft4axz(jj)  = lagint_4th(r4,fr4axz ,rc_p1)
                 ft4ayy(jj)  = lagint_4th(r4,fr4ayy ,rc_p1)
                 ft4ayz(jj)  = lagint_4th(r4,fr4ayz ,rc_p1)
                 ft4azz(jj)  = lagint_4th(r4,fr4azz ,rc_p1)
               end do
               fp4psi(kk)  = lagint_4th(th4,ft4psi ,thc_p1)
               fp4alph(kk) = lagint_4th(th4,ft4alph,thc_p1)
               fp4bvxd(kk) = lagint_4th(th4,ft4bvxd,thc_p1)
               fp4bvyd(kk) = lagint_4th(th4,ft4bvyd,thc_p1)
               fp4bvzd(kk) = lagint_4th(th4,ft4bvzd,thc_p1)
               fp4axx(kk)  = lagint_4th(th4,ft4axx ,thc_p1)
               fp4axy(kk)  = lagint_4th(th4,ft4axy ,thc_p1)
               fp4axz(kk)  = lagint_4th(th4,ft4axz ,thc_p1)
               fp4ayy(kk)  = lagint_4th(th4,ft4ayy ,thc_p1)
               fp4ayz(kk)  = lagint_4th(th4,ft4ayz ,thc_p1)
               fp4azz(kk)  = lagint_4th(th4,ft4azz ,thc_p1)
             end do
             psica  = lagint_4th(phi4,fp4psi ,phic_p1)
             alphca = lagint_4th(phi4,fp4alph,phic_p1)
             bvxdca = lagint_4th(phi4,fp4bvxd,phic_p1)
             bvydca = lagint_4th(phi4,fp4bvyd,phic_p1)
             bvzdca = lagint_4th(phi4,fp4bvzd,phic_p1)
             axx    = lagint_4th(phi4,fp4axx ,phic_p1)
             axy    = lagint_4th(phi4,fp4axy ,phic_p1)
             axz    = lagint_4th(phi4,fp4axz ,phic_p1)
             ayy    = lagint_4th(phi4,fp4ayy ,phic_p1)
             ayz    = lagint_4th(phi4,fp4ayz ,phic_p1)
             azz    = lagint_4th(phi4,fp4azz ,phic_p1)
 
             psi4ca = psica**4
       !      write(6,*) axx,axy,axz,ayy,ayz,azz
 
             call interpo_lag4th_2Dsurf(rsca_p1,rs_p1,thc_p1,phic_p1)
             rcf_p1 = rc_p1/rsca_p1
       !
             if (rcf_p1.le.rg_p1(nrf_p1)) then
               do irg = 0, nrf_p1+1
                 if (rcf_p1.lt.rgex_p1(irg).and.rcf_p1.ge.rgex_p1(irg-1)) ir0 = min0(irg-2,nrf_p1-3)
               end do
       !
               do ii = 1, 4
                 irg0 = ir0 + ii - 1
                 r4(ii) = rgex_p1(irg0)
               end do
       !
               do kk = 1, 4
                 ipg0 = ip0 + kk - 1
                 do jj = 1, 4
                   itg0 = it0 + jj - 1
                   do ii = 1, 4
                     irg0 = ir0 + ii - 1
                     irgex = irgex_r_p1(irg0)
                     itgex = itgex_r_p1(itgex_th_p1(itg0),irg0)
                     ipgex = ipgex_r_p1(ipgex_th_p1(ipgex_phi_p1(ipg0),itg0),irg0)
                     fr4emd(ii)   =   emd_p1(irgex,itgex,ipgex)
                     select case (trim(id_type))
                        case ("IR","SP")
                          fr4vepxf(ii) = vepxf_p1(irgex,itgex,ipgex) !Unique to IR/SP, SP has 3 more
                          fr4vepyf(ii) = vepyf_p1(irgex,itgex,ipgex) !Unique to IR/SP, SP has 3 more
                          fr4vepzf(ii) = vepzf_p1(irgex,itgex,ipgex) !Unique to IR/SP, SP has 3 more
                          select case (trim(id_type))
                          case ("SP")
                              fr4wxspf(ii) = wxspf_p1(irgex,itgex,ipgex)
                              fr4wyspf(ii) = wyspf_p1(irgex,itgex,ipgex)
                              fr4wzspf(ii) = wzspf_p1(irgex,itgex,ipgex)
                          end select
                     end select
                     
 
                     fr4psif(ii)  =  psif_p1(irgex,itgex,ipgex)
                     fr4alphf(ii) = alphf_p1(irgex,itgex,ipgex)
                     fr4bvxdf(ii) = bvxdf_p1(irgex,itgex,ipgex)
                     fr4bvydf(ii) = bvydf_p1(irgex,itgex,ipgex)
                     fr4bvzdf(ii) = bvzdf_p1(irgex,itgex,ipgex)
                   end do
                   ft4emd(jj)   = lagint_4th(r4,fr4emd  ,rcf_p1)
                   select case (trim(id_type))
                      case ("IR","SP")
                        ft4vepxf(jj) = lagint_4th(r4,fr4vepxf,rcf_p1) !Unique to IR, SP has 3 more
                        ft4vepyf(jj) = lagint_4th(r4,fr4vepyf,rcf_p1) !Unique to IR, SP has 3 more
                        ft4vepzf(jj) = lagint_4th(r4,fr4vepzf,rcf_p1) !Unique to IR, SP has 3 more
                      select case (trim(id_type))
                        case ("SP")
                          ft4wxspf(jj) = lagint_4th(r4,fr4wxspf,rcf_p1)
                          ft4wyspf(jj) = lagint_4th(r4,fr4wyspf,rcf_p1)
                          ft4wzspf(jj) = lagint_4th(r4,fr4wzspf,rcf_p1)
                      end select
                   end select
                   
 
                   ft4psif(jj)  = lagint_4th(r4,fr4psif ,rcf_p1)
                   ft4alphf(jj) = lagint_4th(r4,fr4alphf,rcf_p1)
                   ft4bvxdf(jj) = lagint_4th(r4,fr4bvxdf,rcf_p1)
                   ft4bvydf(jj) = lagint_4th(r4,fr4bvydf,rcf_p1)
                   ft4bvzdf(jj) = lagint_4th(r4,fr4bvzdf,rcf_p1)
                 end do
                 fp4emd(kk)   = lagint_4th(th4,ft4emd  ,thc_p1)
                 select case (trim(id_type))
                    case ("IR","SP")
                      fp4vepxf(kk) = lagint_4th(th4,ft4vepxf,thc_p1) !Unique to iR, SP has 3 more
                      fp4vepyf(kk) = lagint_4th(th4,ft4vepyf,thc_p1) !Unique to iR, SP has 3 more
                      fp4vepzf(kk) = lagint_4th(th4,ft4vepzf,thc_p1) !Unique to iR, SP has 3 more
                      select case (trim(id_type))
                        case ("SP")
                          fp4wxspf(kk) = lagint_4th(th4,ft4wxspf,thc_p1)
                          fp4wyspf(kk) = lagint_4th(th4,ft4wyspf,thc_p1)
                          fp4wzspf(kk) = lagint_4th(th4,ft4wzspf,thc_p1)
                      end select
                  end select
                 
 
                 fp4psif(kk)  = lagint_4th(th4,ft4psif ,thc_p1)
                 fp4alphf(kk) = lagint_4th(th4,ft4alphf,thc_p1)
                 fp4bvxdf(kk) = lagint_4th(th4,ft4bvxdf,thc_p1)
                 fp4bvydf(kk) = lagint_4th(th4,ft4bvydf,thc_p1)
                 fp4bvzdf(kk) = lagint_4th(th4,ft4bvzdf,thc_p1)
               end do
               emdca   = lagint_4th(phi4,fp4emd  ,phic_p1)

               select case (trim(id_type))
                  case ("IR","SP")
                    vepxfca = lagint_4th(phi4,fp4vepxf,phic_p1) !Unique to IR, SP has 3 more
                    vepyfca = lagint_4th(phi4,fp4vepyf,phic_p1) !Unique to IR, SP has 3 more
                    vepzfca = lagint_4th(phi4,fp4vepzf,phic_p1) !Unique to IR, SP has 3 more
                    select case (trim(id_type))
                      case ("SP")
                        wxspfca = lagint_4th(phi4,fp4wxspf,phic_p1)
                        wyspfca = lagint_4th(phi4,fp4wyspf,phic_p1)
                        wzspfca = lagint_4th(phi4,fp4wzspf,phic_p1)
                    end select
               end select
 
               psifca  = lagint_4th(phi4,fp4psif ,phic_p1)
               alphfca = lagint_4th(phi4,fp4alphf,phic_p1)
               bvxdfca = lagint_4th(phi4,fp4bvxdf,phic_p1)
               bvydfca = lagint_4th(phi4,fp4bvydf,phic_p1)
               bvzdfca = lagint_4th(phi4,fp4bvzdf,phic_p1)
       !
               psif4ca  = psifca**4
               select case (trim(id_type))
                  case ("IR","SP")
                    alphfca2 = alphfca**2 !Unique to IR, SP has 1 more line here
                    select case (trim(id_type))
                      case ("SP")
                        psifcacp = psifca**confpow_p1
                    end select
               end select
            
 
               bxcor   = bvxdfca + ome_p1*(-ycoc) + coc2cac_ecc_cor_velx * xcoc
               bycor   = bvydfca + ome_p1*(xcoc)  + coc2cac_ecc_cor_velx * ycoc
               bzcor   = bvzdfca                  + coc2cac_ecc_cor_velx * zcoc
               

               select case (trim(id_type))
                  case ("CO")
                    vxu = bxcor/alphfca
                    vyu = bycor/alphfca
                    vzu = bzcor/alphfca
                  case ("IR")
                    lam_p1  = ber_p1 + bxcor*vepxfca + bycor*vepyfca + bzcor*vepzfca !Unique to IR/SP, SP has 3 more lines here
                    huta = lam_p1/alphfca !Unique to IR/SP, differening expressions

                    vxu = (vepxfca/psif4ca )/huta !IR/CO/SP differernt expressions
                    vyu = (vepyfca/psif4ca )/huta !IR/CO/SP differernt expressions
                    vzu = (vepzfca/psif4ca )/huta !IR/CO/SP differernt expressions
                  case ("SP")
                    lam_p1  = ber_p1 + bxcor*vepxfca + bycor*vepyfca + bzcor*vepzfca !Unique to IR/SP, SP has 3 more lines here
                    wdvep   = (wxspfca*vepxfca + wyspfca*vepyfca + wzspfca*vepzfca)*psifcacp
                    w2      = psif4ca*(wxspfca*wxspfca + wyspfca*wyspfca + wzspfca*wzspfca)*psifcacp**2.0d0
                    wterm   = wdvep + w2 !wdvep = 0 then huta_sp = huta_ir
                    huta = (lam_p1 + sqrt(lam_p1*lam_p1 + 4.0d0*alphfca2*wterm))/(2.0d0*alphfca)

                    vxu = (vepxfca/psif4ca + psifcacp*wxspfca)/huta
                    vyu = (vepyfca/psif4ca + psifcacp*wyspfca)/huta
                    vzu = (vepzfca/psif4ca + psifcacp*wzspfca)/huta
               end select
 
            
             else
               vepxfca=0.0d0; vepyfca=0.0d0; vepzfca=0.0d0 !Unique to IR/SP
               emdca=0.0d0;
               vxu=0.0d0; vyu=0.0d0; vzu=0.0d0
             end if
           else
       !=>      call copy_from_mpatch_interpolation_utility(2)
       !      call copy_grid_parameter_from_mpt(2)
       !      call copy_grid_parameter_binary_excision_from_mpt(2)
       !      call copy_coordinate_grav_extended_from_mpt(2)
       !      call copy_coordinate_grav_phi_from_mpt(2)
       !      call copy_coordinate_grav_r_from_mpt(2)
       !      call copy_coordinate_grav_theta_from_mpt(2)
       !      call copy_def_binary_parameter_from_mpt(2)
       !      call copy_trigonometry_grav_phi_from_mpt(2)
       !      call copy_trigonometry_grav_theta_from_mpt(2)
 
             xc_p2 = -(xcoc - dis_cm)
             yc_p2 = -ycoc
             zc_p2 =  zcoc
            !  write(6,'(a39,3e20.12)') "Point given wrt COCP-2 in Cocal scale :", xc_p2,yc_p2,zc_p2
       !
             psica=0.0d0; alphca=0.0d0; bvxdca=0.0d0; bvydca=0.0d0; bvzdca=0.0d0
             axx=0.0d0  ; axy=0.0d0   ; axz=0.0d0   ; ayy=0.0d0   ; ayz=0.0d0   ; azz=0.0d0
       !
             rc_p2     = dsqrt(dabs(xc_p2**2 + yc_p2**2 + zc_p2**2))
             varpic_p2 = dsqrt(dabs(xc_p2**2 + yc_p2**2))
             thc_p2  = dmod(2.0d0*pi + datan2(varpic_p2,zc_p2),2.0d0*pi)
             phic_p2 = dmod(2.0d0*pi + datan2(    yc_p2,xc_p2),2.0d0*pi)
       !
             do irg = 0, nrg_p2+1
               if (rc_p2.lt.rgex_p2(irg).and.rc_p2.ge.rgex_p2(irg-1)) ir0 = min0(irg-2,nrg_p2-3)
             end do
             do itg = 0, ntg_p2+1
               if (thc_p2.lt.thgex_p2(itg).and.thc_p2.ge.thgex_p2(itg-1)) it0 = itg-2
             end do
             do ipg = 0, npg_p2+1
               if (phic_p2.lt.phigex_p2(ipg).and.phic_p2.ge.phigex_p2(ipg-1)) ip0 = ipg-2
             end do
       !
             do ii = 1, 4
               irg0 = ir0 + ii - 1
               itg0 = it0 + ii - 1
               ipg0 = ip0 + ii - 1
               r4(ii) = rgex_p2(irg0)
               th4(ii) = thgex_p2(itg0)
               phi4(ii) = phigex_p2(ipg0)
             end do
       !
             do kk = 1, 4
               ipg0 = ip0 + kk - 1
               do jj = 1, 4
                 itg0 = it0 + jj - 1
                 do ii = 1, 4
                   irg0 = ir0 + ii - 1
                   irgex = irgex_r_p2(irg0)
                   itgex = itgex_r_p2(itgex_th_p2(itg0),irg0)
                   ipgex = ipgex_r_p2(ipgex_th_p2(ipgex_phi_p2(ipg0),itg0),irg0)
                   fr4psi(ii)  =  psi_p2(irgex,itgex,ipgex)
                   fr4alph(ii) = alph_p2(irgex,itgex,ipgex)
                   fr4bvxd(ii) = bvxd_p2(irgex,itgex,ipgex)
                   fr4bvyd(ii) = bvyd_p2(irgex,itgex,ipgex)
                   fr4bvzd(ii) = bvzd_p2(irgex,itgex,ipgex)
                   fr4axx(ii)  =  axx_p2(irgex,itgex,ipgex)
                   fr4axy(ii)  =  axy_p2(irgex,itgex,ipgex)
                   fr4axz(ii)  =  axz_p2(irgex,itgex,ipgex)
                   fr4ayy(ii)  =  ayy_p2(irgex,itgex,ipgex)
                   fr4ayz(ii)  =  ayz_p2(irgex,itgex,ipgex)
                   fr4azz(ii)  =  azz_p2(irgex,itgex,ipgex)
                 end do
                 ft4psi(jj)  = lagint_4th(r4,fr4psi ,rc_p2)
                 ft4alph(jj) = lagint_4th(r4,fr4alph,rc_p2)
                 ft4bvxd(jj) = lagint_4th(r4,fr4bvxd,rc_p2)
                 ft4bvyd(jj) = lagint_4th(r4,fr4bvyd,rc_p2)
                 ft4bvzd(jj) = lagint_4th(r4,fr4bvzd,rc_p2)
                 ft4axx(jj)  = lagint_4th(r4,fr4axx ,rc_p2)
                 ft4axy(jj)  = lagint_4th(r4,fr4axy ,rc_p2)
                 ft4axz(jj)  = lagint_4th(r4,fr4axz ,rc_p2)
                 ft4ayy(jj)  = lagint_4th(r4,fr4ayy ,rc_p2)
                 ft4ayz(jj)  = lagint_4th(r4,fr4ayz ,rc_p2)
                 ft4azz(jj)  = lagint_4th(r4,fr4azz ,rc_p2)
               end do
               fp4psi(kk)  = lagint_4th(th4,ft4psi ,thc_p2)
               fp4alph(kk) = lagint_4th(th4,ft4alph,thc_p2)
               fp4bvxd(kk) = lagint_4th(th4,ft4bvxd,thc_p2)
               fp4bvyd(kk) = lagint_4th(th4,ft4bvyd,thc_p2)
               fp4bvzd(kk) = lagint_4th(th4,ft4bvzd,thc_p2)
               fp4axx(kk)  = lagint_4th(th4,ft4axx ,thc_p2)
               fp4axy(kk)  = lagint_4th(th4,ft4axy ,thc_p2)
               fp4axz(kk)  = lagint_4th(th4,ft4axz ,thc_p2)
               fp4ayy(kk)  = lagint_4th(th4,ft4ayy ,thc_p2)
               fp4ayz(kk)  = lagint_4th(th4,ft4ayz ,thc_p2)
               fp4azz(kk)  = lagint_4th(th4,ft4azz ,thc_p2)
             end do
             psica  = lagint_4th(phi4,fp4psi ,phic_p2)
             alphca = lagint_4th(phi4,fp4alph,phic_p2)
             bvxdca = lagint_4th(phi4,fp4bvxd,phic_p2)
             bvydca = lagint_4th(phi4,fp4bvyd,phic_p2)
             bvzdca = lagint_4th(phi4,fp4bvzd,phic_p2)
             axx    = lagint_4th(phi4,fp4axx ,phic_p2)
             axy    = lagint_4th(phi4,fp4axy ,phic_p2)
             axz    = lagint_4th(phi4,fp4axz ,phic_p2)
             ayy    = lagint_4th(phi4,fp4ayy ,phic_p2)
             ayz    = lagint_4th(phi4,fp4ayz ,phic_p2)
             azz    = lagint_4th(phi4,fp4azz ,phic_p2)
 
             psi4ca = psica**4
             call interpo_lag4th_2Dsurf(rsca_p2,rs_p2,thc_p2,phic_p2)
             rcf_p2 = rc_p2/rsca_p2
       !
             if (rcf_p2.le.rg_p2(nrf_p2)) then
               do irg = 0, nrf_p2+1
                 if (rcf_p2.lt.rgex_p2(irg).and.rcf_p2.ge.rgex_p2(irg-1)) ir0 = min0(irg-2,nrf_p2-3)
               end do
       !
               do ii = 1, 4
                 irg0 = ir0 + ii - 1
                 r4(ii) = rgex_p2(irg0)
               end do
       !
               do kk = 1, 4
                 ipg0 = ip0 + kk - 1
                 do jj = 1, 4
                   itg0 = it0 + jj - 1
                   do ii = 1, 4
                     irg0 = ir0 + ii - 1
                     irgex = irgex_r_p2(irg0)
                     itgex = itgex_r_p2(itgex_th_p2(itg0),irg0)
                     ipgex = ipgex_r_p2(ipgex_th_p2(ipgex_phi_p2(ipg0),itg0),irg0)
                     fr4emd(ii)   =   emd_p2(irgex,itgex,ipgex)
                     select case (trim(id_type))
                        case ("IR","SP")
                          fr4vepxf(ii) = vepxf_p2(irgex,itgex,ipgex) !Unique to IR/SP, SP has 3 more
                          fr4vepyf(ii) = vepyf_p2(irgex,itgex,ipgex) !Unique to IR/SP, SP has 3 more
                          fr4vepzf(ii) = vepzf_p2(irgex,itgex,ipgex) !Unique to IR/SP, SP has 3 more
                          select case (trim(id_type))
                          case ("SP")
                            fr4wxspf(ii) = wxspf_p2(irgex,itgex,ipgex)
                            fr4wyspf(ii) = wyspf_p2(irgex,itgex,ipgex)
                            fr4wzspf(ii) = wzspf_p2(irgex,itgex,ipgex)
                          end select
                     end select
 
                     fr4psif(ii)  =  psif_p2(irgex,itgex,ipgex)
                     fr4alphf(ii) = alphf_p2(irgex,itgex,ipgex)
                     fr4bvxdf(ii) = bvxdf_p2(irgex,itgex,ipgex)
                     fr4bvydf(ii) = bvydf_p2(irgex,itgex,ipgex)
                     fr4bvzdf(ii) = bvzdf_p2(irgex,itgex,ipgex)
                   end do
                   ft4emd(jj)   = lagint_4th(r4,fr4emd  ,rcf_p2)
                   select case (trim(id_type))
                     case ("IR","SP")
                        ft4vepxf(jj) = lagint_4th(r4,fr4vepxf,rcf_p2) !Unique to IR, SP has 3 more
                        ft4vepyf(jj) = lagint_4th(r4,fr4vepyf,rcf_p2) !Unique to IR, SP has 3 more
                        ft4vepzf(jj) = lagint_4th(r4,fr4vepzf,rcf_p2) !Unique to IR, SP has 3 more
                        select case (trim(id_type))
                          case ("SP")
                            ft4wxspf(jj) = lagint_4th(r4,fr4wxspf,rcf_p2)
                            ft4wyspf(jj) = lagint_4th(r4,fr4wyspf,rcf_p2)
                            ft4wzspf(jj) = lagint_4th(r4,fr4wzspf,rcf_p2)
                        end select
                   end select
 
                   ft4psif(jj)  = lagint_4th(r4,fr4psif ,rcf_p2)
                   ft4alphf(jj) = lagint_4th(r4,fr4alphf,rcf_p2)
                   ft4bvxdf(jj) = lagint_4th(r4,fr4bvxdf,rcf_p2)
                   ft4bvydf(jj) = lagint_4th(r4,fr4bvydf,rcf_p2)
                   ft4bvzdf(jj) = lagint_4th(r4,fr4bvzdf,rcf_p2)
                 end do
                 fp4emd(kk)   = lagint_4th(th4,ft4emd  ,thc_p2)
                 select case (trim(id_type))
                    case ("IR","SP")
                      fp4vepxf(kk) = lagint_4th(th4,ft4vepxf,thc_p2) !Unique to iR, SP has 3 more
                      fp4vepyf(kk) = lagint_4th(th4,ft4vepyf,thc_p2) !Unique to iR, SP has 3 more
                      fp4vepzf(kk) = lagint_4th(th4,ft4vepzf,thc_p2) !Unique to iR, SP has 3 more
                      select case (trim(id_type))
                        case ("SP")
                          fp4wxspf(kk) = lagint_4th(th4,ft4wxspf,thc_p2)
                          fp4wyspf(kk) = lagint_4th(th4,ft4wyspf,thc_p2)
                          fp4wzspf(kk) = lagint_4th(th4,ft4wzspf,thc_p2)
                      end select
                  end select
 
                 fp4psif(kk)  = lagint_4th(th4,ft4psif ,thc_p2)
                 fp4alphf(kk) = lagint_4th(th4,ft4alphf,thc_p2)
                 fp4bvxdf(kk) = lagint_4th(th4,ft4bvxdf,thc_p2)
                 fp4bvydf(kk) = lagint_4th(th4,ft4bvydf,thc_p2)
                 fp4bvzdf(kk) = lagint_4th(th4,ft4bvzdf,thc_p2)
               end do
               emdca   = lagint_4th(phi4,fp4emd  ,phic_p2)
               select case (trim(id_type))
                  case ("IR","SP")
                    vepxfca = lagint_4th(phi4,fp4vepxf,phic_p2) !Unique to IR, SP has 3 more
                    vepyfca = lagint_4th(phi4,fp4vepyf,phic_p2) !Unique to IR, SP has 3 more
                    vepzfca = lagint_4th(phi4,fp4vepzf,phic_p2) !Unique to IR, SP has 3 more
                    select case (trim(id_type))
                      case ("SP")
                        wxspfca = lagint_4th(phi4,fp4wxspf,phic_p2)
                        wyspfca = lagint_4th(phi4,fp4wyspf,phic_p2)
                        wzspfca = lagint_4th(phi4,fp4wzspf,phic_p2)
                    end select
               end select
 
               psifca  = lagint_4th(phi4,fp4psif ,phic_p2)
               alphfca = lagint_4th(phi4,fp4alphf,phic_p2)
               bvxdfca = lagint_4th(phi4,fp4bvxdf,phic_p2)
               bvydfca = lagint_4th(phi4,fp4bvydf,phic_p2)
               bvzdfca = lagint_4th(phi4,fp4bvzdf,phic_p2)
       !
               psif4ca = psifca**4
               select case (trim(id_type))
                  case ("IR","SP")
                    alphfca2 = alphfca**2 !Unique to IR, SP has 1 more line here
                    select case (trim(id_type))
                      case ("SP")
                        psifcacp = psifca**confpow_p2
                    end select
               end select
            
 
               bxcor   = bvxdfca + ome_p2*(-ycoc) + coc2cac_ecc_cor_velx * xcoc
               bycor   = bvydfca + ome_p2*(xcoc)  + coc2cac_ecc_cor_velx * ycoc
               bzcor   = bvzdfca                  + coc2cac_ecc_cor_velx * zcoc
               
       !
              select case (trim(id_type))
                case ("CO")
                  vxu = bxcor/alphfca
                  vyu = bycor/alphfca
                  vzu = bzcor/alphfca
                case ("IR")
                  lam_p2  = ber_p2 + bxcor*vepxfca + bycor*vepyfca + bzcor*vepzfca !Unique to IR, SP has 3 more lines here
                  huta = lam_p2/alphfca !Unique to IR/SP, differening expressions

                  vxu = (vepxfca/psif4ca )/huta !IR/CO/SP differernt expressions
                  vyu = (vepyfca/psif4ca )/huta !IR/CO/SP differernt expressions
                  vzu = (vepzfca/psif4ca )/huta !IR/CO/SP differernt expressions
                case ("SP")
                  lam_p2  = ber_p2 + bxcor*vepxfca + bycor*vepyfca + bzcor*vepzfca !Unique to IR, SP has 3 more lines here
                  wdvep   = (wxspfca*vepxfca + wyspfca*vepyfca + wzspfca*vepzfca)*psifcacp
                  w2      = psif4ca*(wxspfca*wxspfca + wyspfca*wyspfca + wzspfca*wzspfca)*psifcacp**2.0d0
                  wterm   = wdvep + w2
                  huta = (lam_p2 + sqrt(lam_p2*lam_p2 + 4.0d0*alphfca2*wterm))/(2.0d0*alphfca)

                  vxu = (vepxfca/psif4ca + psifcacp*wxspfca)/huta
                  vyu = (vepyfca/psif4ca + psifcacp*wyspfca)/huta
                  vzu = (vepzfca/psif4ca + psifcacp*wzspfca)/huta
              end select
             else
               vepxfca=0.0d0; vepyfca=0.0d0; vepzfca=0.0d0 !Unique to IR
               emdca=0.0d0;
               vxu=0.0d0; vyu=0.0d0; vzu=0.0d0
             end if
       !
           end if
         end if
         
         gxx1 = psi4ca
         gxy1 = 0.0d0
         gxz1 = 0.0d0
         gyy1 = psi4ca
         gyz1 = 0.0d0
         gzz1 = psi4ca
 
         kxx1 = psi4ca*axx/(radi_p1)
         kxy1 = psi4ca*axy/(radi_p1)
         kxz1 = psi4ca*axz/(radi_p1)
         kyy1 = psi4ca*ayy/(radi_p1)
         kyz1 = psi4ca*ayz/(radi_p1)
         kzz1 = psi4ca*azz/(radi_p1)
 
         gxx(i,j,k) = gxx1
         gxy(i,j,k) = gxy1
         gxz(i,j,k) = gxz1
         gyy(i,j,k) = gyy1
         gyz(i,j,k) = gyz1
         gzz(i,j,k) = gzz1
 
         kxx(i,j,k) = kxx1
         kxy(i,j,k) = kxy1
         kxz(i,j,k) = kxz1
         kyy(i,j,k) = kyy1
         kyz(i,j,k) = kyz1
         kzz(i,j,k) = kzz1
 
         if (gxx1 == gxx1) then
           continue
         else 
           call CCTK_WARN(CCTK_WARN_ABORT, "NaN in gxx")
         end if 
 
         call peos_q2hprho(emdca, hca, preca, rhoca, eneca)
         if (bool_lapse) then
            alp(i,j,k) = alphca
         end if
         if (bool_shift) then
            betax(i,j,k) = bvxdca
            betay(i,j,k) = bvydca
            betaz(i,j,k) = bvzdca
         end if
         if (bool_hydro) then
            rho(i,j,k) = rhoca
            press(i,j,k) = preca
            eps(i,j,k) = eneca/rhoca - 1.0d0
            vel(i,j,k,1) = vxu
            vel(i,j,k,2) = vyu
            vel(i,j,k,3) = vzu
         end if
         if (bool_Bvec) then
            Bvec(i,j,k,1) = 0.0d0
            Bvec(i,j,k,2) = 0.0d0
            Bvec(i,j,k,3) = 0.0d0

            !Avec(i,j,k,1) = 0.0d0
            !Avec(i,j,k,2) = 0.0d0
            !Avec(i,j,k,3) = 0.0d0

            !Aphi(i,j,k) = 0.0d0
         end if
         
 
       !   write(6,'(a6,e20.12)') "psi  =", psica
       !   write(6,'(a6,e20.12)') "alph =", alphca
       !   write(6,'(a6,e20.12)') "bvxd =", bvxdca
       !   write(6,'(a6,e20.12)') "bvyd =", bvydca
       !   write(6,'(a6,e20.12)') "bvzd =", bvzdca
       !   write(6,'(a6,e20.12)') "Radi =", r_surf_p1*radi_p1
       !   write(6,'(a6,e20.12)') "Omeg =", ome_p1/radi_p1
       !   write(6,'(a6,e20.12)') "emd  =", emdca
       !   write(6,'(a6,e20.12)') "h    =", hca
       !   write(6,'(a6,e20.12)') "pre  =", preca
       !   write(6,'(a6,e20.12)') "rho  =", rhoca
       !   write(6,'(a6,e20.12)') "ene  =", eneca
       !   write(6,'(a6,e20.12)') "eps  =", epsca
       !   write(6,'(a6,e20.12)') "vepx =", vepxfca
       !   write(6,'(a6,e20.12)') "vepy =", vepyfca
       !   write(6,'(a6,e20.12)') "vepz =", vepzfca
       ! !
       !   write(6,'(a18)') "Kij at gridpoints:"
       !   write(6,'(3e20.12)') kxx, kxy, kxz
       !   write(6,'(3e20.12)') kxy, kyy, kyz
       !   write(6,'(3e20.12)') kxz, kyz, kzz
 
       !   write(6,'(a13)') "v^i eulerian:"
       !   write(6,'(a6,e20.12)') "vxu  =", vxu
       !   write(6,'(a6,e20.12)') "vyu  =", vyu
       !   write(6,'(a6,e20.12)') "vzu  =", vzu
       end do
     end do
   end do
 
   call CCTK_INFO("Finished looping over Cactus Grid")
 

 
   call CCTK_INFO("Deallocating....")
   deallocate(       rg_p1);  deallocate(       rg_p2);  deallocate(       rg_p3)
   deallocate(     rgex_p1);  deallocate(     rgex_p2);  deallocate(     rgex_p3)
   deallocate(    thgex_p1);  deallocate(    thgex_p2);  deallocate(    thgex_p3)
   deallocate(   phigex_p1);  deallocate(   phigex_p2);  deallocate(   phigex_p3)
   deallocate(  irgex_r_p1);  deallocate(  irgex_r_p2);  deallocate(  irgex_r_p3)
   deallocate( itgex_th_p1);  deallocate( itgex_th_p2);  deallocate( itgex_th_p3)
   deallocate(ipgex_phi_p1);  deallocate(ipgex_phi_p2);  deallocate(ipgex_phi_p3)
   deallocate(  itgex_r_p1);  deallocate(  itgex_r_p2);  deallocate(  itgex_r_p3)
   deallocate(  ipgex_r_p1);  deallocate(  ipgex_r_p2);  deallocate(  ipgex_r_p3)
   deallocate( ipgex_th_p1);  deallocate( ipgex_th_p2);  deallocate( ipgex_th_p3)
   
   deallocate(  emd_p1);  deallocate(  emd_p2);
   select case (trim(id_type))
      case ("IR","SP")
         deallocate(vep_p1);    deallocate(vep_p2);
         deallocate(vepxf_p1);  deallocate(vepxf_p2); 
         deallocate(vepyf_p1);  deallocate(vepyf_p2); 
         deallocate(vepzf_p1);  deallocate(vepzf_p2); 
         select case (trim(id_type))
            case ("SP")
               deallocate(wxspf_p1);  deallocate(wxspf_p2); 
               deallocate(wyspf_p1);  deallocate(wyspf_p2); 
               deallocate(wzspf_p1);  deallocate(wzspf_p2); 
         end select
   end select
   deallocate( psif_p1);  deallocate( psif_p2);
   deallocate(alphf_p1);  deallocate(alphf_p2);
   deallocate(bvxdf_p1);  deallocate(bvxdf_p2);
   deallocate(bvydf_p1);  deallocate(bvydf_p2);
   deallocate(bvzdf_p1);  deallocate(bvzdf_p2);
   deallocate(   rs_p1);  deallocate(   rs_p2);
   deallocate(  psi_p1);  deallocate(  psi_p2);  deallocate(psi_p3)
   deallocate( alph_p1);  deallocate( alph_p2);  deallocate(alph_p3)
   deallocate( bvxd_p1);  deallocate( bvxd_p2);  deallocate(bvxd_p3)
   deallocate( bvyd_p1);  deallocate( bvyd_p2);  deallocate(bvyd_p3)
   deallocate( bvzd_p1);  deallocate( bvzd_p2);  deallocate(bvzd_p3)
   deallocate(  axx_p1);  deallocate(  axx_p2);  deallocate(axx_p3)
   deallocate(  axy_p1);  deallocate(  axy_p2);  deallocate(axy_p3)
   deallocate(  axz_p1);  deallocate(  axz_p2);  deallocate(axz_p3)
   deallocate(  ayy_p1);  deallocate(  ayy_p2);  deallocate(ayy_p3)
   deallocate(  ayz_p1);  deallocate(  ayz_p2);  deallocate(ayz_p3)
   deallocate(  azz_p1);  deallocate(  azz_p2);  deallocate(azz_p3)
   call CCTK_INFO("Finished Deallocating...")
   call CCTK_INFO("Finished coc2cac_bns master subroutine")
 !
   contains
   integer function read_id_type(filename,id_type)
     implicit none
     integer :: nrg, nrf, nrf_deform, nrgin
     character(len=*) :: filename
     character(2) :: id_type, EQ_point
   
     open(unit = 1,file = trim(filename), status='old')
     read(1,'(4i5)') nrg
     read(1,'(4i5)') nrf
     read(1,'(2i5,2(3x,a2))') nrf_deform, nrgin, id_type, EQ_point
     close(1)
     read_id_type = 0
   end function read_id_type
 END SUBROUTINE coc2cac_bns