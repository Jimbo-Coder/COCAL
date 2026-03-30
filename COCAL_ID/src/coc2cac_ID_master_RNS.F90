!_____________________________________________________________________________
!
!    CACTUS READER OF COCAL ROTATING STAR IN WAVELESS FORMALISM 
!_____________________________________________________________________________
!

#define module_use use

#include "Cocal/phys_constant.f90"
#include "Cocal/def_matter_parameter.f90"
#include "Cocal/def_quantities.f90"
#include "Cocal/def_bh_parameter.f90"

#include "Cocal/make_array_2d.f90"
#include "Cocal/grid_parameter.f90"
#include "Cocal/interface_modules_cartesian.f90"
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
#include "Cocal/interface_invhij_WL_export.f90"
#include "Cocal/interface_index_vec_down2up_export.f90"
#include "Cocal/interface_interpo_gr2fl_metric_CF_export.f90"
#include "Cocal/interface_IO_input_grav_export_Ai.f90"
#include "Cocal/interface_IO_input_grav_export_Faraday.f90"
#include "Cocal/interface_IO_input_star4ve_export.f90"

#include "Cocal/IO_input_CF_grav_export.f90"
#include "Cocal/IO_input_WL_grav_export_hij.f90"
#include "Cocal/IO_input_grav_export_Kij.f90"


#include "Cocal/interface_IO_input_matter_BHT_export.f90"


#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Functions.h"
#include "cctk_Parameters.h"

!RNS Routine
subroutine coc2cac_rns(CCTK_ARGUMENTS)
  use phys_constant
  use grid_parameter, eps_cocal => eps
  use interface_modules_cartesian
  use coordinate_grav_r
  use coordinate_grav_phi
  use coordinate_grav_theta
  use coordinate_grav_extended
  use trigonometry_grav_theta
  use trigonometry_grav_phi
  use interface_IO_input_CF_grav_export
  use interface_IO_input_WL_grav_export_hij
  use interface_IO_input_grav_export_Kij
  use interface_IO_input_CF_star_export
  use interface_invhij_WL_export
  use interface_index_vec_down2up_export
  use interface_interpo_gr2fl_metric_CF_export

  use interface_IO_input_matter_BHT_export ! Unique to BHT
 

  use interface_IO_input_grav_export_Ai            ! only used in MRNS_WL
  use interface_IO_input_grav_export_Faraday       !
  use interface_IO_input_star4ve_export            !
  
  

  implicit none
  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_FUNCTIONS
  DECLARE_CCTK_PARAMETERS
  !call CCTK_INFO("Imported modules and CCTK stuff")
  integer :: iAB 
  !character(30) :: char1
  character*400 :: dir_path, outstr, coc2cac_readformatf !,cocout_dir
  

  integer :: dir_path_len,coc2cac_,coc2cac_readformatlen !, cocout_dir_len
  integer :: i, j, k, nx, ny, nz, l, ns

  logical :: bool_hydro, bool_lapse, bool_shift, bool_Bvec
  
  real(8) :: xcac, ycac, zcac
  real(8) :: xcoc, ycoc, zcoc, rcoc, expnr, expnv
  real(8) :: emdca,  omefca, psica,  alphca, psi4ca, psif4ca
  real(8) :: bvxdca, bvydca, bvzdca, bvxuca, bvyuca, bvzuca
  real(8) :: hxxdca, hxydca, hxzdca, hyydca, hyzdca, hzzdca
  real(8) :: hxxuca, hxyuca, hxzuca, hyyuca, hyzuca, hzzuca
  real(8) :: hca, preca, rhoca, eneca, epsca
  real(8) :: kxxca, kxyca, kxzca, kyyca, kyzca, kzzca
  real(8) :: vxu, vyu, vzu
  real(8) :: bxcor, bycor, bzcor, bvxufca, bvyufca, bvzufca, psifca, alphfca
  real(8) :: gxx1, gxy1, gxz1, gyy1, gyz1, gzz1, kxx1, kxy1, kxz1, kyy1, kyz1, kzz1
  real(8) :: ome, ber, radi, rexc

  real(8) :: va1, vaxd1, vayd1, vazd1, fxd1, fyd1, fzd1, fxyd1, fxzd1, fyzd1!
  real(8) :: vaca, vaxdca, vaydca, vazdca, fxdca, fydca, fzdca, fxydca, fxzdca, fyzdca!
  real(8) :: utfca, uxfca, uyfca, uzfca ! Fluid 4 vel
  real(8), pointer :: utf(:,:,:)  ,  uxf(:,:,:) ,  uyf(:,:,:) ,  uzf(:,:,:)!
!
  real(8), pointer :: emd(:,:,:), omef(:,:,:), rs(:,:)
  real(8), pointer :: psif(:,:,:), alphf(:,:,:), bvxuf(:,:,:), bvyuf(:,:,:), bvzuf(:,:,:)
  real(8), pointer :: psi(:,:,:)  , alph(:,:,:)
  real(8), pointer :: bvxd(:,:,:) , bvyd(:,:,:) , bvzd(:,:,:) , bvxu(:,:,:) , bvyu(:,:,:), bvzu(:,:,:)
  real(8), pointer :: hxxd(:,:,:) , hxyd(:,:,:) , hxzd(:,:,:) , hyyd(:,:,:) , hyzd(:,:,:), hzzd(:,:,:)
  real(8), pointer :: hxxu(:,:,:) , hxyu(:,:,:) , hxzu(:,:,:) , hyyu(:,:,:) , hyzu(:,:,:), hzzu(:,:,:)
  real(8), pointer :: kxxa(:,:,:)  , kxya(:,:,:)  , kxza(:,:,:)  , kyya(:,:,:)  , kyza(:,:,:) , kzza(:,:,:)

  real(8), pointer :: alphN(:,:,:,:), betaxdN(:,:,:,:), betaydN(:,:,:,:), betazdN(:,:,:,:) ! New stuff
  real(8), pointer :: gxxN(:,:,:,:), gxyN(:,:,:,:), gxzN(:,:,:,:), gyyN(:,:,:,:), gyzN(:,:,:,:), gzzN(:,:,:,:)
  real(8), pointer :: kxxN(:,:,:,:), kxyN(:,:,:,:), kxzN(:,:,:,:), kyyN(:,:,:,:), kyzN(:,:,:,:), kzzN(:,:,:,:)
  real(8), pointer :: rhoN(:,:,:,:), pressN(:,:,:,:), epsN(:,:,:,:)
  real(8), pointer :: veluxN(:,:,:,:), veluyN(:,:,:,:), veluzN(:,:,:,:), BvecuxN(:,:,:,:), BvecuyN(:,:,:,:), BvecuzN(:,:,:,:)

  real(8), pointer ::  va(:,:,:)  , vaxd(:,:,:) , vayd(:,:,:) , vazd(:,:,:)                           ! (phi,A)
  real(8), pointer :: fxd(:,:,:)  ,  fyd(:,:,:) ,  fzd(:,:,:) , fxyd(:,:,:) , fxzd(:,:,:), fyzd(:,:,:) ! Faraday tensor
  call CCTK_INFO("Cocal matrices and points declared")
!
  gxx1=0.0d0; gxy1=0.0d0; gxz1=0.0d0; gyy1=0.0d0; gyz1=0.0d0; gzz1=0.0d0
  kxx1=0.0d0; kxy1=0.0d0; kxz1=0.0d0; kyy1=0.0d0; kyz1=0.0d0; kzz1=0.0d0
  kxxca=0.0d0; kxyca=0.0d0; kxzca=0.0d0; kyyca=0.0d0; kyzca=0.0d0; kzzca=0.0d0

!   alphN = 0.0d0; betaxdN = 0.0d0; betaydN = 0.0d0; betazdN = 0.0d0
!   gxxN = 0.0d0; gxyN = 0.0d0; gxzN = 0.0d0; gyyN = 0.0d0; gyzN = 0.0d0; gzzN = 0.0d0
!   kxxN = 0.0d0; kxyN = 0.0d0; kxzN = 0.0d0; kyyN = 0.0d0; kyzN = 0.0d0; kzzN = 0.0d0
!   rhoN = 0.0d0; pressN = 0.0d0; epsN = 0.0d0
!   veluxN = 0.0d0; veluyN = 0.0d0; veluzN = 0.0d0; BvecuxN = 0.0d0; BvecuyN = 0.0d0; BvecuzN = 0.0d0

  vaca=0.0d0;  vaxdca=0.0d0;  vaydca=0.0d0;  vazdca=0.0d0;                              ! (phi,A)
  fxdca=0.0d0; fydca=0.0d0; fzdca=0.0d0; fxydca=0.0d0; fxzdca=0.0d0; fyzdca=0.0d0;     ! Faraday Tensor
  fxd1=0.0d0;  fyd1=0.0d0;  fzd1=0.0d0;  fxyd1=0.0d0;  fxzd1=0.0d0; fyzd1=0.0d0;          !
  va1=0.0d0;  vaxd1=0.0d0; vayd1=0.0d0; vazd1=0.0d0;                                      !(phi,A) (Aphi, Avec[3]) & (Bvec[3])
  utfca=0.0d0;  uxfca=0.0d0;  uyfca=0.0d0;  uzfca=0.0d0;                                  ! Fluid 4 vel (vel[3])

  expnr=0.0d0; expnv=0.0d0

  nx = cctk_lsh(1)
  ny = cctk_lsh(2)
  nz = cctk_lsh(3)
  ns = coc2cac_Nstars

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

  if (ns > 1) then 
    call CCTK_INFO("Beginning multiple star reader... ")
    allocate(alphN(1:nx,1:ny,1:nz,1:ns))
    allocate(betaxdN(1:nx,1:ny,1:nz,1:ns))
    allocate(betaydN(1:nx,1:ny,1:nz,1:ns))
    allocate(betazdN(1:nx,1:ny,1:nz,1:ns))
    allocate(gxxN(1:nx,1:ny,1:nz,1:ns))
    allocate(gxyN(1:nx,1:ny,1:nz,1:ns))
    allocate(gxzN(1:nx,1:ny,1:nz,1:ns))
    allocate(gyyN(1:nx,1:ny,1:nz,1:ns))
    allocate(gyzN(1:nx,1:ny,1:nz,1:ns))
    allocate(gzzN(1:nx,1:ny,1:nz,1:ns))
    allocate(kxxN(1:nx,1:ny,1:nz,1:ns))
    allocate(kxyN(1:nx,1:ny,1:nz,1:ns))
    allocate(kxzN(1:nx,1:ny,1:nz,1:ns))
    allocate(kyyN(1:nx,1:ny,1:nz,1:ns))
    allocate(kyzN(1:nx,1:ny,1:nz,1:ns))
    allocate(kzzN(1:nx,1:ny,1:nz,1:ns))
    allocate(rhoN(1:nx,1:ny,1:nz,1:ns))
    allocate(pressN(1:nx,1:ny,1:nz,1:ns))
    allocate(epsN(1:nx,1:ny,1:nz,1:ns))
    allocate(veluxN(1:nx,1:ny,1:nz,1:ns))
    allocate(veluyN(1:nx,1:ny,1:nz,1:ns))
    allocate(veluzN(1:nx,1:ny,1:nz,1:ns))
    allocate(BvecuxN(1:nx,1:ny,1:nz,1:ns))
    allocate(BvecuyN(1:nx,1:ny,1:nz,1:ns))
    allocate(BvecuzN(1:nx,1:ny,1:nz,1:ns))
  end if
  call CCTK_INFO("Beginning loop over all stars:")
  do l = 1, ns
      call CCTK_FortranString(dir_path_len,coc2cac_dir_path_ID,dir_path)  !For reading initial data

      !call CCTK_FortranString(cocout_dir_len,out_dir,cocout_dir) !For copy to output directory

      call CCTK_INFO("Beginning cocal_rns_thorn")
      
      !write(outstr,'(5A)') coc2cac_rnstype
      call CCTK_INFO("Execuiting Reader: "); call CCTK_PrintString(coc2cac_rnstype);

      call CCTK_FortranString(coc2cac_readformatlen,coc2cac_readformat,coc2cac_readformatf)


      call CCTK_INFO("Reading format: " // coc2cac_readformatf)
      ! -- Read parameters
      call CCTK_INFO("Cocal: Reading parameters...")
      

      ! -- Read parameters
      call read_parameter_cactus(dir_path)
      call peos_initialize_cactus(dir_path)
      if (CCTK_EQUALS(coc2cac_rnstype, "BHT_WL")) then
         call read_bht_parameter_cactus(dir_path) ! Unique to BHT
         call calc_bht_excision_radius ! Unique to BHT
         call grid_r_bht('eBH') ! Dfferent naming
         rexc = rg(0)
         write(outstr,'(E10.4)') rexc
         call CCTK_INFO("Excision at r="//outstr)
      else
         call grid_r
      end if
      call grid_theta
      call trig_grav_theta
      call grid_phi
      call allocate_trig_grav_mphi
      call trig_grav_phi
      call grid_extended

      ! -- Allocate arrays
      call CCTK_INFO("Cocal: Allocating local arrays...")

      write(outstr,'(6i5)') nrg, ntg, npg, nrf, ntf, npf
      call CCTK_INFO("Cocal Grav grid, Fluid Grid Sizes: "//outstr)
      !
      !    write(6,'(6i5)') nrg, ntg, npg, nrf, ntf, npf
      !  rr3 = 0.7d0*(rgout - rgmid)
      !  dis_cm = dis
      if (.not. CCTK_EQUALS(coc2cac_rnstype, "BHT_WL")) then
         allocate (  emd(0:nrf,0:ntf,0:npf))
         allocate ( omef(0:nrf,0:ntf,0:npf))
      else
         allocate (  emd(0:nrg,0:ntg,0:npg))
         allocate ( omef(0:nrg,0:ntg,0:npg))
      end if
      allocate ( psif(0:nrf,0:ntf,0:npf))
      allocate (alphf(0:nrf,0:ntf,0:npf))
      allocate (bvxuf(0:nrf,0:ntf,0:npf))
      allocate (bvyuf(0:nrf,0:ntf,0:npf))
      allocate (bvzuf(0:nrf,0:ntf,0:npf))
      allocate (   rs(0:ntf,0:npf))
      call CCTK_INFO("Fluid Grid Allocated...")
      allocate (  psi(0:nrg,0:ntg,0:npg))
      allocate ( alph(0:nrg,0:ntg,0:npg))
      allocate ( bvxd(0:nrg,0:ntg,0:npg))
      allocate ( bvyd(0:nrg,0:ntg,0:npg))
      allocate ( bvzd(0:nrg,0:ntg,0:npg))
      allocate ( bvxu(0:nrg,0:ntg,0:npg))
      allocate ( bvyu(0:nrg,0:ntg,0:npg))
      allocate ( bvzu(0:nrg,0:ntg,0:npg))
      allocate ( hxxd(0:nrg,0:ntg,0:npg))
      allocate ( hxyd(0:nrg,0:ntg,0:npg))
      allocate ( hxzd(0:nrg,0:ntg,0:npg))
      allocate ( hyyd(0:nrg,0:ntg,0:npg))
      allocate ( hyzd(0:nrg,0:ntg,0:npg))
      allocate ( hzzd(0:nrg,0:ntg,0:npg))
      allocate ( hxxu(0:nrg,0:ntg,0:npg))
      allocate ( hxyu(0:nrg,0:ntg,0:npg))
      allocate ( hxzu(0:nrg,0:ntg,0:npg))
      allocate ( hyyu(0:nrg,0:ntg,0:npg))
      allocate ( hyzu(0:nrg,0:ntg,0:npg))
      allocate ( hzzu(0:nrg,0:ntg,0:npg))
      allocate (  kxxa(0:nrg,0:ntg,0:npg))
      allocate (  kxya(0:nrg,0:ntg,0:npg))
      allocate (  kxza(0:nrg,0:ntg,0:npg))
      allocate (  kyya(0:nrg,0:ntg,0:npg))
      allocate (  kyza(0:nrg,0:ntg,0:npg))
      allocate (  kzza(0:nrg,0:ntg,0:npg))
      call CCTK_INFO("Gravitational Grid Allocated...")

      if (CCTK_EQUALS(coc2cac_rnstype,"MRNS_WL")) then
         allocate (  utf(0:nrf,0:ntf,0:npf)) ! fluid 4vel
         allocate (  uxf(0:nrf,0:ntf,0:npf)) !
         allocate (  uyf(0:nrf,0:ntf,0:npf)) !
         allocate (  uzf(0:nrf,0:ntf,0:npf)) !
         utf=0.0d0;  uxf=0.0d0;   uyf=0.0d0;   uzf=0.0d0;     !fluid 4 vel init
         if (coc2cac_readpot == 1) then
            allocate (   va(0:nrg,0:ntg,0:npg)) !(phi,A)
            allocate ( vaxd(0:nrg,0:ntg,0:npg))
            allocate ( vayd(0:nrg,0:ntg,0:npg))
            allocate ( vazd(0:nrg,0:ntg,0:npg))
            va=0.0d0;    vaxd=0.0d0;  vayd=0.0d0;   vazd=0.0d0;
         else if (coc2cac_readpot == 0) then
            allocate (  fxd(0:nrg,0:ntg,0:npg)) !Faraday tensor
            allocate (  fyd(0:nrg,0:ntg,0:npg))
            allocate (  fzd(0:nrg,0:ntg,0:npg))
            allocate ( fxyd(0:nrg,0:ntg,0:npg))
            allocate ( fxzd(0:nrg,0:ntg,0:npg))
            allocate ( fyzd(0:nrg,0:ntg,0:npg))
            fxd=0.0d0;   fyd=0.0d0;   fzd=0.0d0;    fxyd=0.0d0;   fxzd=0.0d0;  fyzd=0.0d0;
         end if
      end if
      call CCTK_INFO("Done Allocating...")
      emd=0.0d0;  rs  =0.0d0;  omef=0.0d0

      
      psif=0.0d0; alphf=0.0d0; bvxuf=0.0d0; bvyuf=0.0d0; bvzuf=0.0d0 !

      psi=0.0d0;  alph=0.0d0;  bvxd=0.0d0;  bvyd=0.0d0;  bvzd=0.0d0
      bvxu=0.0d0; bvyu=0.0d0;  bvzu=0.0d0
      kxxa=0.0d0;  kxya =0.0d0;  kxza =0.0d0;   kyya=0.0d0;   kyza=0.0d0;   kzza=0.0d0
      hxxd=0.0d0; hxyd=0.0d0;  hxzd=0.0d0;  hyyd=0.0d0;  hyzd=0.0d0;  hzzd=0.0d0;
      hxxu=0.0d0; hxyu=0.0d0;  hxzu=0.0d0;  hyyu=0.0d0;  hyzu=0.0d0;  hzzu=0.0d0;
      call CCTK_INFO("Done initializing Varriables...")

      call IO_input_CF_grav_export(trim(dir_path)//"/rnsgra_3D.las",coc2cac_readformatf,psi,alph,bvxd,bvyd,bvzd)

      if (CCTK_EQUALS(coc2cac_rnstype, "BHT_WL")) then
         call IO_input_matter_BHT_export(trim(dir_path)//"/rnsflu_3D.las",coc2cac_readformatf,emd,omef,ome,ber,radi)
      else
         call IO_input_CF_star_export(trim(dir_path)//"/rnsflu_3D.las",coc2cac_readformatf,emd,rs,omef,ome,ber,radi)
      end if

      call IO_input_grav_export_Kij(trim(dir_path)//"/rnsgra_Kij_3D.las",kxxa,kxya,kxza,kyya,kyza,kzza)

      ! LET CF h_ij's remain initialized to zero.
      ! Only CF change is that h_ij = 0, everything else remains the same.
      ! Only non-trivial change between CF and WL
      if (CCTK_EQUALS(coc2cac_rnstype, "RNS_WL") .or. CCTK_EQUALS(coc2cac_rnstype,"MRNS_WL") .or. (CCTK_EQUALS(coc2cac_rnstype, "BHT_WL"))) then
         call IO_input_WL_grav_export_hij(trim(dir_path)//"/rnsgra_hij_3D.las",coc2cac_readformatf,hxxd,hxyd,hxzd,hyyd,hyzd,hzzd)
      else if (CCTK_EQUALS(coc2cac_rnstype, "RNS_CF")) then 
         continue
      else 
         call CCTK_WARN(CCTK_WARN_ABORT, "Invalid cocal_rns_thorn: coc2cac_rnstype")
      end if


      if (CCTK_EQUALS(coc2cac_rnstype,"MRNS_WL")) then
         if (coc2cac_readpot == 1) then
            call IO_input_grav_export_Ai(trim(dir_path)//"/rnsEMF_3D.las",coc2cac_readformatf,va,vaxd,vayd,vazd)
         else if (coc2cac_readpot == 0) then
            call IO_input_grav_export_Faraday(trim(dir_path)//"/rnsEMF_faraday_3D.las",coc2cac_readformatf,fxd,fyd,fzd,fxyd,fxzd,fyzd)
         end if
         call IO_input_star4ve_export(trim(dir_path)//"/rns4ve_3D.las",coc2cac_readformatf,utf,uxf,uyf,uzf) ! fluid 4vel
      end if
      
         
      call invhij_WL_export(hxxd,hxyd,hxzd,hyyd,hyzd,hzzd,hxxu,hxyu,hxzu,hyyu,hyzu,hzzu) !trivial for CF 0=0

      call index_vec_down2up_export(hxxu,hxyu,hxzu,hyyu,hyzu,hzzu,bvxu,bvyu,bvzu,bvxd,bvyd,bvzd) !trivial for CF, u=d

      if (.not. CCTK_EQUALS(coc2cac_rnstype, "BHT_WL")) then
         call interpo_gr2fl_metric_CF_export(alph, psi, bvxu, bvyu, bvzu, alphf, psif, bvxuf, bvyuf, bvzuf, rs) !Main Interpolation
      end if
      

      !  call excurve_CF_gridpoint_export(alph,bvxd,bvyd,bvzd, & 
      !     &    axx, axy, axz, ayy, ayz, azz)
      

      write(6,'(2e20.12)') emd(0,0,0), omef(0,0,0)
      write(6,'(3e20.12)') ome, ber, radi
      
      !


      call CCTK_INFO("Cocal: Looping over local cartesian grid:")

      
      
      
         do k = 1, nz
            do j = 1, ny
               do i = 1, nx
               if (l == 1) then
                  xcac = x(i,j,k) - coc2cac_cen_x1
                  ycac = y(i,j,k) - coc2cac_cen_y1
                  zcac = z(i,j,k) - coc2cac_cen_z1
               else if (l == 2) then
                  xcac = x(i,j,k) - coc2cac_cen_x2
                  ycac = y(i,j,k) - coc2cac_cen_y2
                  zcac = z(i,j,k) - coc2cac_cen_z2
               else if (l == 3) then
                  xcac = x(i,j,k) - coc2cac_cen_x3
                  ycac = y(i,j,k) - coc2cac_cen_y3
                  zcac = z(i,j,k) - coc2cac_cen_z3
               end if
               !vcordx in carpetx for cell vertex
               !ccordx in carpetx for cell center
         !        write(6,*)' i, j, k, xcac, ycac, zcac', i, j, k, xcac, ycac, zcac

               xcoc = xcac/(radi)
               ycoc = ycac/(radi)
               zcoc = zcac/(radi)
               rcoc = sqrt(xcoc**2 + ycoc**2 + zcoc**2)

               if (CCTK_EQUALS(coc2cac_rnstype, "BHT_WL")) then
                  if ((rcoc < rexc) .and. (bool_hydro)) then
                     press(i,j,k) = 0.0
                     eps(i,j,k)   = 0.0
                     vel(i,j,k,1) = 0.0
                     vel(i,j,k,2) = 0.0
                     vel(i,j,k,3) = 0.0
                     Bvec(i,j,k,1) = 0.0
                     Bvec(i,j,k,2) = 0.0
                     Bvec(i,j,k,3) = 0.0
                  end if
               end if
         !        write(6,'(a23,3e20.12)') "Point given wrt COCAL:", xcoc,ycoc,zcoc

               call interpo_gr2cgr_4th(psi , psica , xcoc, ycoc, zcoc)
               call interpo_gr2cgr_4th(alph, alphca, xcoc, ycoc, zcoc)

               call interpo_gr2cgr_4th(bvxd, bvxdca, xcoc, ycoc, zcoc)
               call interpo_gr2cgr_4th(bvyd, bvydca, xcoc, ycoc, zcoc)
               call interpo_gr2cgr_4th(bvzd, bvzdca, xcoc, ycoc, zcoc)

               call interpo_gr2cgr_4th(bvxu, bvxuca, xcoc, ycoc, zcoc) !Repeats itself for CF
               call interpo_gr2cgr_4th(bvyu, bvyuca, xcoc, ycoc, zcoc)
               call interpo_gr2cgr_4th(bvzu, bvzuca, xcoc, ycoc, zcoc) 
         
               call interpo_gr2cgr_4th(hxxd, hxxdca, xcoc, ycoc, zcoc) !Trivial for CF, 0 = 0
               call interpo_gr2cgr_4th(hxyd, hxydca, xcoc, ycoc, zcoc)
               call interpo_gr2cgr_4th(hxzd, hxzdca, xcoc, ycoc, zcoc)
               call interpo_gr2cgr_4th(hyyd, hyydca, xcoc, ycoc, zcoc)
               call interpo_gr2cgr_4th(hyzd, hyzdca, xcoc, ycoc, zcoc)
               call interpo_gr2cgr_4th(hzzd, hzzdca, xcoc, ycoc, zcoc)
               
               
               call interpo_gr2cgr_4th(kxxa , kxxca , xcoc, ycoc, zcoc)
               call interpo_gr2cgr_4th(kxya , kxyca , xcoc, ycoc, zcoc)
               call interpo_gr2cgr_4th(kxza , kxzca , xcoc, ycoc, zcoc)
               call interpo_gr2cgr_4th(kyya , kyyca , xcoc, ycoc, zcoc)
               call interpo_gr2cgr_4th(kyza , kyzca , xcoc, ycoc, zcoc)
               call interpo_gr2cgr_4th(kzza , kzzca , xcoc, ycoc, zcoc)

               
               if (CCTK_EQUALS(coc2cac_rnstype, "BHT_WL")) then
                  call interpo_gr2cgr_4th(emd  , emdca   , xcoc, ycoc, zcoc)
                  call interpo_gr2cgr_4th(omef , omefca  , xcoc, ycoc, zcoc)
               else
                  call interpo_fl2cgr_4th_export(emd  , emdca   , xcoc, ycoc, zcoc, rs)
                  call interpo_fl2cgr_4th_export(omef , omefca  , xcoc, ycoc, zcoc, rs)

                  call interpo_fl2cgr_4th_export(bvxuf, bvxufca , xcoc, ycoc, zcoc, rs)
                  call interpo_fl2cgr_4th_export(bvyuf, bvyufca , xcoc, ycoc, zcoc, rs)
                  call interpo_fl2cgr_4th_export(bvzuf, bvzufca , xcoc, ycoc, zcoc, rs)

                  call interpo_fl2cgr_4th_export(psif , psifca  , xcoc, ycoc, zcoc, rs)
                  call interpo_fl2cgr_4th_export(alphf, alphfca , xcoc, ycoc, zcoc, rs)
               end if

               if (CCTK_EQUALS(coc2cac_rnstype,"MRNS_WL")) then
                  call interpo_fl2cgr_4th_export(utf, utfca , xcoc, ycoc, zcoc, rs)!4vel
                  call interpo_fl2cgr_4th_export(uxf, uxfca , xcoc, ycoc, zcoc, rs)!
                  call interpo_fl2cgr_4th_export(uyf, uyfca , xcoc, ycoc, zcoc, rs)!
                  call interpo_fl2cgr_4th_export(uzf, uzfca , xcoc, ycoc, zcoc, rs)!
                  if (coc2cac_readpot == 1) then
                     call interpo_gr2cgr_4th(  va,   vaca, xcoc, ycoc, zcoc)!potentials
                     call interpo_gr2cgr_4th(vaxd, vaxdca, xcoc, ycoc, zcoc)
                     call interpo_gr2cgr_4th(vayd, vaydca, xcoc, ycoc, zcoc)
                     call interpo_gr2cgr_4th(vazd, vazdca, xcoc, ycoc, zcoc)
                  else if (coc2cac_readpot == 0) then
                     call interpo_gr2cgr_4th(fxd,   fxdca, xcoc, ycoc, zcoc) !Faraday/fields
                     call interpo_gr2cgr_4th(fyd,   fydca, xcoc, ycoc, zcoc)
                     call interpo_gr2cgr_4th(fzd,   fzdca, xcoc, ycoc, zcoc)
                     call interpo_gr2cgr_4th(fxyd, fxydca, xcoc, ycoc, zcoc)
                     call interpo_gr2cgr_4th(fxzd, fxzdca, xcoc, ycoc, zcoc)
                     call interpo_gr2cgr_4th(fyzd, fyzdca, xcoc, ycoc, zcoc)
                  end if
               end if
               
               if ( CCTK_EQUALS(coc2cac_rnstype, "BHT_WL")) then
                  bxcor = bvxuca + omefca*(-ycoc)
                  bycor = bvyuca + omefca*(xcoc)
                  bzcor = bvzuca
               else
                  bxcor = bvxufca + omefca*(-ycoc)
                  bycor = bvyufca + omefca*(xcoc)
                  bzcor = bvzufca
                  psif4ca = psifca**4
               end if
               
               psi4ca = psica**4
               

               if (dabs(emdca) > 1.0d-14) then
                  if (CCTK_EQUALS(coc2cac_rnstype, "RNS_WL") .or. CCTK_EQUALS(coc2cac_rnstype,"RNS_CF")) then
                     vxu = bxcor/alphfca 
                     vyu = bycor/alphfca
                     vzu = bzcor/alphfca
                  else if  (CCTK_EQUALS(coc2cac_rnstype, "BHT_WL")) then
                     vxu = bxcor/alphca
                     vyu = bycor/alphca
                     vzu = bzcor/alphca
                  else if (CCTK_EQUALS(coc2cac_rnstype, "MRNS_WL")) then
                     vxu = (uxfca/utfca + bvxufca)/alphfca     !
                     vyu = (uyfca/utfca + bvyufca)/alphfca     !
                     vzu = (uzfca/utfca + bvzufca)/alphfca     !
                  end if
               else
                  emdca=0.0d0
                  vxu=0.0d0; vyu=0.0d0; vzu=0.0d0
               end if


               
               gxx1 = psi4ca*(1.0d0+hxxdca)
               gxy1 = psi4ca*(0.0d0+hxydca)
               gxz1 = psi4ca*(0.0d0+hxzdca)
               gyy1 = psi4ca*(1.0d0+hyydca)
               gyz1 = psi4ca*(0.0d0+hyzdca)
               gzz1 = psi4ca*(1.0d0+hzzdca)

               kxx1 = psi4ca*kxxca/(radi)
               kxy1 = psi4ca*kxyca/(radi)
               kxz1 = psi4ca*kxzca/(radi)
               kyy1 = psi4ca*kyyca/(radi)
               kyz1 = psi4ca*kyzca/(radi)
               kzz1 = psi4ca*kzzca/(radi)
               if (emdca .ne. emdca) then
                  emdca = 0.0d0 !NaN fixing condition, pray atmospheric reset fixes things afterwards.
               end if

               call peos_q2hprho(emdca, hca, preca, rhoca, eneca)
               !

               if (ns ==1) then
                  if (bool_Bvec) then
                     if (CCTK_EQUALS(coc2cac_rnstype, "MRNS_WL")) then
                        if (coc2cac_readpot == 1) then
                           Aphi(i,j,k)   = vaca
                           Avec(i,j,k,1) = vaxdca
                           Avec(i,j,k,2) = vaydca
                           Avec(i,j,k,3) = vazdca
                        else if (coc2cac_readpot == 0) then
                           fxd1  = fxdca/radi
                           fyd1  = fydca/radi
                           fzd1  = fzdca/radi
                           fxyd1 = fxydca/radi
                           fxzd1 = fxzdca/radi
                           fyzd1 = fyzdca/radi
                           
                           Bvec(i,j,k,1) = fxyd1 ! Probably not right
                           Bvec(i,j,k,2) = fxzd1
                           Bvec(i,j,k,3) = fyzd1
                           
                        end if
                     else 
                        Bvec(i,j,k,1) = 0.0d0
                        Bvec(i,j,k,2) = 0.0d0
                        Bvec(i,j,k,3) = 0.0d0
                     end if 
                  end if
            
                    

                  !SET ETK VARIABLES HERE
                  if (bool_lapse) then
                     if (alphca == alphca) then
                        alp(i,j,k) = alphca
                     else 
                        alp(i,j,k) = 1.0d0
                     end if
                  end if
                 
                  if (bool_shift) then
                     if (bvxdca == bvxdca) then
                        betax(i,j,k) = bvxdca
                     else 
                        betax(i,j,k) = 0.0d0
                     end if
                     if (bvydca == bvydca) then
                        betay(i,j,k) = bvydca
                     else 
                        betay(i,j,k) = 0.0d0
                     end if
                     if (bvzdca == bvzdca) then
                        betaz(i,j,k) = bvzdca
                     else 
                        betaz(i,j,k) = 0.0d0
                     end if
                  end if
                  if (bool_hydro) then
                     rho(i,j,k) = rhoca
                     press(i,j,k) = preca
                     eps(i,j,k) = eneca/rhoca - 1.0d0
                     vel(i,j,k,1) = vxu
                     vel(i,j,k,2) = vyu
                     vel(i,j,k,3) = vzu
                  end if


                  gxx(i,j,k) = gxx1
                  gxy(i,j,k) = gxy1
                  gxz(i,j,k) = gxz1
                  gyy(i,j,k) = gyy1
                  gyz(i,j,k) = gyz1
                  gzz(i,j,k) = gzz1
                  if (gxx1 == gxx1) then
                     continue
                  else 
                     call CCTK_WARN(CCTK_WARN_ABORT, "NaN in gxx")
                  end if 
                  
                  kxx(i,j,k) = kxx1
                  kxy(i,j,k) = kxy1
                  kxz(i,j,k) = kxz1
                  kyy(i,j,k) = kyy1
                  kyz(i,j,k) = kyz1
                  kzz(i,j,k) = kzz1
               else if (ns>1) then
                                       !ETK variables set here
                  alphn(i,j,k,l) = alphca
                  betaxdN(i,j,k,l) = bvxdca
                  betaydN(i,j,k,l) = bvydca
                  betazdN(i,j,k,l) = bvzdca

                  rhoN(i,j,k,l) = rhoca
                  pressN(i,j,k,l) = preca
                  epsN(i,j,k,l) = eneca/rhoca - 1.0d0
                  veluxN(i,j,k,l) = vxu
                  veluyN(i,j,k,l) = vyu
                  veluzN(i,j,k,l) = vzu

                  if (CCTK_EQUALS(coc2cac_rnstype, "MRNS_WL")) then

                     BvecuxN(i,j,k,l) = fxyd1 ! Probably not right
                     BvecuyN(i,j,k,l) = fxzd1
                     BvecuzN(i,j,k,l) = fyzd1
                  end if


                  gxxN(i,j,k,l) = gxx1
                  gxyN(i,j,k,l) = gxy1
                  gxzN(i,j,k,l) = gxz1
                  gyyN(i,j,k,l) = gyy1
                  gyzN(i,j,k,l) = gyz1
                  gzzN(i,j,k,l) = gzz1
                  if (gxx1 == gxx1) then
                     continue
                  else 
                     call CCTK_WARN(CCTK_WARN_ABORT, "NaN in gxx")
                  end if 
                  
                  kxxN(i,j,k,l) = kxx1
                  kxyN(i,j,k,l) = kxy1
                  kxzN(i,j,k,l) = kxz1
                  kyyN(i,j,k,l) = kyy1
                  kyzN(i,j,k,l) = kyz1
                  kzzN(i,j,k,l) = kzz1

               end if
               
               end do
            end do
         end do
         
            
      call CCTK_INFO("3D Loop 100% Done")
     

      call CCTK_INFO("Deallocating....")
      deallocate(  emd);  deallocate( omef);  deallocate( psif);  deallocate(alphf);    
      deallocate(bvxuf);  deallocate(bvyuf);  deallocate(bvzuf);  deallocate(   rs);    
      deallocate(  psi);  deallocate( alph);  deallocate( bvxd);  deallocate( bvyd);  
      deallocate( bvzd);  deallocate( bvxu);  deallocate( bvyu);  deallocate( bvzu);
      deallocate( hxxd);  deallocate( hxyd);  deallocate( hxzd);  deallocate( hyyd);  
      deallocate( hyzd);  deallocate( hzzd);  deallocate( hxxu);  deallocate( hxyu);  
      deallocate( hxzu);  deallocate( hyyu);  deallocate( hyzu);  deallocate( hzzu);  
      deallocate(  kxxa);  deallocate(  kxya);  deallocate(  kxza);  deallocate(  kyya);  
      deallocate(  kyza);  deallocate(  kzza);  

      
      
      if (CCTK_EQUALS(coc2cac_rnstype, "MRNS_WL")) then  
         deallocate(  utf);  deallocate(  uxf);  deallocate(  uyf);  deallocate(  uzf);!
         if (coc2cac_readpot == 1) then
            deallocate(   va);  deallocate( vaxd);  deallocate( vayd);  deallocate( vazd);
         else if (coc2cac_readpot == 0 ) then
            deallocate(  fxd);  deallocate(  fyd);  deallocate(  fzd);
            deallocate( fxyd);  deallocate( fxzd);  deallocate( fyzd);
         end if
      end if
      !
      call CCTK_INFO("Cocal: Done Deallocating...")
   end do
   if (ns > 1) then
      call CCTK_INFO("Superimposing multiple solutions onto a single one, weighting by lapse function...")
      if (bool_lapse) then
         alp = SUM(alphN,4) - 1.0d0*(ns - 1.0d0) ! Whoops, forgot to subtract ns-1
      end if
      if (bool_shift) then
         betax = SUM(betaxdN,4)
         betay = SUM(betaydN,4)
         betaz = SUM(betazdN,4)
      end if

      if (bool_Bvec) then
         if (CCTK_EQUALS(coc2cac_rnstype, "MRNS_WL")) then
            Bvec(:,:,:,1) = SUM(BvecuxN, 4) ! Probably not right
            Bvec(:,:,:,2) = SUM(BvecuyN, 4)
            Bvec(:,:,:,3) = SUM(BvecuzN, 4)
         else 
            Bvec(:,:,:,1) = 0.0d0
            Bvec(:,:,:,2) = 0.0d0
            Bvec(:,:,:,3) = 0.0d0
         end if
      end if

      if (bool_hydro) then
         !rho = SUM(rhoN, 4) ! Regular superposition
         rho = (alp**expnr)*SUM( (alphN**(-expnr))*rhoN, 4) ! Weighted superposition
         press = SUM(pressN, 4)
         eps = SUM(epsN, 4)
         !vel(:,:,:,1) = SUM(veluxN, 4) !Regular superposition
         !vel(:,:,:,2) = SUM(veluyN, 4) !Regular superposition
         !vel(:,:,:,3) = SUM(veluzN, 4) !Regular superposition
         vel(:,:,:,1) = (alp**(-expnv))*SUM( (alphN**expnv)*veluxN, 4) !Weighted superposition
         vel(:,:,:,2) = (alp**(-expnv))*SUM( (alphN**expnv)*veluyN, 4) !Weighted superposition
         vel(:,:,:,3) = (alp**(-expnv))*SUM( (alphN**expnv)*veluzN, 4) !Weighted superposition
      end if


      gxx = SUM(gxxN, 4) - 1.0d0*(ns - 1.0d0)
      gxy = SUM(gxyN, 4)
      gxz = SUM(gxzN, 4)
      gyy = SUM(gyyN, 4) - 1.0d0*(ns - 1.0d0)
      gyz = SUM(gyzN, 4)
      gzz = SUM(gzzN, 4) - 1.0d0*(ns - 1.0d0)
      if (gxx1 == gxx1) then
         continue
      else 
         call CCTK_WARN(CCTK_WARN_ABORT, "NaN in gxx")
      end if 
      
      kxx = SUM(kxxN,4)
      kxy = SUM(kxyN,4)
      kxz = SUM(kxzN,4)
      kyy = SUM(kyyN,4)
      kyz = SUM(kyzN,4)
      kzz = SUM(kzzN,4)

      deallocate(alphN)
      deallocate(betaxdN)
      deallocate(betaydN)
      deallocate(betazdN)
      deallocate(gxxN)
      deallocate(gxyN)
      deallocate(gxzN)
      deallocate(gyyN)
      deallocate(gyzN)
      deallocate(gzzN)
      deallocate(kxxN)
      deallocate(kxyN)
      deallocate(kxzN)
      deallocate(kyyN)
      deallocate(kyzN)
      deallocate(kzzN)
      deallocate(rhoN)
      deallocate(pressN)
      deallocate(epsN)
      deallocate(veluxN)
      deallocate(veluyN)
      deallocate(veluzN)
      deallocate(BvecuxN)
      deallocate(BvecuyN)
      deallocate(BvecuzN)
   end if


! if (coc2cac_copyID == 1) then
!   call CCTK_INFO("Copying Initial Data to Output Directory...")
!   call execute_command_line("cp -r " // trim(dir_path) // " " // trim(cocout_dir) // "/")
! end if

call CCTK_INFO("cocal_rns_thorn: finishing...")
END subroutine coc2cac_rns