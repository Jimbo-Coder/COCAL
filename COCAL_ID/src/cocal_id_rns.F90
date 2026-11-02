! Rotating neutron stars (RNS_CF, RNS_WL) and magnetized ones (MRNS_WL): one
! patch, with the fluid on COCAL's surface-fitted grid.

module COCAL_ID_rns
  use COCAL_ID_state
  implicit none
  private
  public :: rns_em, rns_em0, RNS_EM_POTENTIAL, RNS_EM_FARADAY, COCAL_ID_rns_em

  ! MRNS electromagnetic fields follow the metric fields, from index rns_em0:
  ! the potential A_t, A_x, A_y, A_z, or the Faraday tensor F_xy, F_xz, F_yz.
  integer, parameter :: RNS_EM_POTENTIAL = 1, RNS_EM_FARADAY = 2
  integer, save :: rns_em = 0, rns_em0 = 0

contains

  subroutine COCAL_ID_rns_em(q, s)
    real(8), intent(in) :: q(:)
    real(8), intent(inout) :: s(NSTATE)

    select case (rns_em)
    case (RNS_EM_POTENTIAL)
      s(S_APHI) = q(rns_em0)
      s(S_AVEC:S_AVEC + 2) = q(rns_em0 + 1:rns_em0 + 3)
    case (RNS_EM_FARADAY)
      s(S_BVEC:S_BVEC + 2) = [q(rns_em0 + 2), -q(rns_em0 + 1), q(rns_em0)]/cocal_radi
    end select
  end subroutine COCAL_ID_rns_em

end module COCAL_ID_rns

! Collective over the node, once COCAL's grid modules hold the patch's grids:
! the metric of a one-patch export from rnsgra_3D.las,
! rnsgra_Kij_3D.las and, for waveless data, rnsgra_hij_3D.las, with nextra
! further fields for the caller to read. COCAL exports the covariant shift of
! the conformal metric; the metric holds beta^i.
subroutine COCAL_ID_read_one_patch_metric(dir, fmt, wl, nextra)
  use COCAL_ID_grid_parameter, only : nrg, ntg, npg
  use COCAL_ID_coordinate_grav_r, only : rg
  use COCAL_ID_coordinate_grav_theta, only : thg
  use COCAL_ID_coordinate_grav_phi, only : phig
  use COCAL_ID_interface_IO_input_CF_grav_export
  use COCAL_ID_interface_IO_input_WL_grav_export_hij
  use COCAL_ID_interface_IO_input_grav_export_Kij
  use COCAL_ID_interface_invhij_WL_export
  use COCAL_ID_interface_index_vec_down2up_export
  use COCAL_ID_node_shared
  use COCAL_ID_state
  use COCAL_ID_domain
  implicit none
  character(*), intent(in) :: dir, fmt
  logical, intent(in) :: wl
  integer, intent(in) :: nextra
  real(8), pointer :: psi(:,:,:), alph(:,:,:), bvxu(:,:,:), bvyu(:,:,:), bvzu(:,:,:)
  real(8), pointer :: kxx(:,:,:), kxy(:,:,:), kxz(:,:,:), kyy(:,:,:), kyz(:,:,:), kzz(:,:,:)
  real(8), pointer :: hxxd(:,:,:), hxyd(:,:,:), hxzd(:,:,:), hyyd(:,:,:), hyzd(:,:,:), hzzd(:,:,:)
  real(8), pointer :: bvxd(:,:,:), bvyd(:,:,:), bvzd(:,:,:)
  real(8), pointer :: hxxu(:,:,:), hxyu(:,:,:), hxzu(:,:,:), hyyu(:,:,:), hyzu(:,:,:), hzzu(:,:,:)

  cocal_wl = wl
  call COCAL_ID_domain_init(cocal_metric, 1, merge(NMETRIC_WL, NMETRIC_CF, wl) + nextra, .true.)
  call COCAL_ID_patch_set_grids(cocal_metric%patch(1), rg(0:nrg), thg(0:ntg), phig(0:npg))
  cocal_metric%patch(1)%positive_indices = [M_PSI, M_ALP]
  if (.not. COCAL_ID_node_leader()) return

  ! Each name below is a view of one field of raw, which COCAL's reading
  ! routines fill in place and the spline fit then consumes.
  call COCAL_ID_domain_alloc_raw(cocal_metric, 1)
  associate(raw => cocal_metric%patch(1)%raw)
    psi(0:,0:,0:)  => raw(M_PSI,:,:,:);   alph(0:,0:,0:) => raw(M_ALP,:,:,:)
    bvxu(0:,0:,0:) => raw(M_BETA,:,:,:);  bvyu(0:,0:,0:) => raw(M_BETA+1,:,:,:); bvzu(0:,0:,0:) => raw(M_BETA+2,:,:,:)
    kxx(0:,0:,0:)  => raw(M_K,:,:,:);     kxy(0:,0:,0:)  => raw(M_K+1,:,:,:);    kxz(0:,0:,0:)  => raw(M_K+2,:,:,:)
    kyy(0:,0:,0:)  => raw(M_K+3,:,:,:);   kyz(0:,0:,0:)  => raw(M_K+4,:,:,:);    kzz(0:,0:,0:)  => raw(M_K+5,:,:,:)
    if (wl) then
      hxxd(0:,0:,0:) => raw(M_H,:,:,:);   hxyd(0:,0:,0:) => raw(M_H+1,:,:,:);    hxzd(0:,0:,0:) => raw(M_H+2,:,:,:)
      hyyd(0:,0:,0:) => raw(M_H+3,:,:,:); hyzd(0:,0:,0:) => raw(M_H+4,:,:,:);    hzzd(0:,0:,0:) => raw(M_H+5,:,:,:)
    end if
  end associate
  call COCAL_ID_IO_input_grav_export_Kij(dir//"/rnsgra_Kij_3D.las", kxx, kxy, kxz, kyy, kyz, kzz)
  if (.not. wl) then
    ! Conformally flat: beta_i = beta^i.
    call COCAL_ID_IO_input_CF_grav_export(dir//"/rnsgra_3D.las", fmt, psi, alph, bvxu, bvyu, bvzu)
    return
  end if
  ! beta^i = gt^ij beta_j, with gt^ij and beta_i needed only here.
  allocate(bvxd, bvyd, bvzd, hxxu, hxyu, hxzu, hyyu, hyzu, hzzu, mold=psi)
  call COCAL_ID_IO_input_CF_grav_export(dir//"/rnsgra_3D.las", fmt, psi, alph, bvxd, bvyd, bvzd)
  call COCAL_ID_IO_input_WL_grav_export_hij(dir//"/rnsgra_hij_3D.las", fmt, hxxd, hxyd, hxzd, hyyd, hyzd, hzzd)
  call COCAL_ID_invhij_WL_export(hxxd, hxyd, hxzd, hyyd, hyzd, hzzd, hxxu, hxyu, hxzu, hyyu, hyzu, hzzu)
  call COCAL_ID_index_vec_down2up_export(hxxu, hxyu, hxzu, hyyu, hyzu, hzzu, bvxu, bvyu, bvzu, bvxd, bvyd, bvzd)
  deallocate(bvxd, bvyd, bvzd, hxxu, hxyu, hxzu, hyyu, hyzu, hzzu)
end subroutine COCAL_ID_read_one_patch_metric

! Collective over the node.
subroutine COCAL_ID_read_rns(dir)
  ! COCAL's parameters, grids and the routines that build them.
  use COCAL_ID_grid_parameter, eps_cocal => eps
  use COCAL_ID_coordinate_grav_r
  use COCAL_ID_coordinate_grav_theta
  use COCAL_ID_coordinate_grav_phi
  use COCAL_ID_trigonometry_grav_theta
  use COCAL_ID_trigonometry_grav_phi
  use COCAL_ID_interface_IO_input_CF_star_export
  use COCAL_ID_interface_IO_input_star4ve_export
  use COCAL_ID_interface_IO_input_grav_export_Ai
  use COCAL_ID_interface_IO_input_grav_export_Faraday
  use COCAL_ID_node_shared
  use COCAL_ID_state
  use COCAL_ID_domain
  use COCAL_ID_rns
  implicit none
  DECLARE_CCTK_PARAMETERS
  DECLARE_CCTK_FUNCTIONS
  character(*), intent(in) :: dir
  ! COCAL's routines take the directory as character*400.
  character(400) :: cocal_dir
  character(*), parameter :: fmt = "20.12"
  real(8), pointer :: emd(:,:,:), omef(:,:,:), rs(:,:), utf(:,:,:), uxf(:,:,:), uyf(:,:,:), uzf(:,:,:)
  real(8), pointer :: va(:,:,:), vaxd(:,:,:), vayd(:,:,:), vazd(:,:,:)
  real(8), pointer :: fxd(:,:,:), fyd(:,:,:), fzd(:,:,:), fxyd(:,:,:), fxzd(:,:,:), fyzd(:,:,:)
  real(8) :: ome, ber, radi
  logical :: wl, mrns, leader, ok
  integer :: nextra

  leader = COCAL_ID_node_leader()
  cocal_dir = dir
  wl = .not. CCTK_EQUALS(COCAL_ID_rnstype, "RNS_CF")
  mrns = CCTK_EQUALS(COCAL_ID_rnstype, "MRNS_WL")
  call read_parameter_cactus(cocal_dir)
  if (cocal_tabulated_eos) then
    call COCAL_ID_teos_initialize_cactus(cocal_dir)
  else
    call COCAL_ID_peos_initialize_cactus(cocal_dir)
  end if
  call grid_r
  call grid_theta
  call trig_grav_theta
  call grid_phi
  call allocate_trig_grav_mphi
  call trig_grav_phi
  if (ntf /= ntg .or. npf /= npg) &
    call CCTK_WARN(CCTK_WARN_ABORT, "COCAL_ID RNS: fluid and gravitational angular grids must match.")

  nextra = 0
  if (mrns) then
    rns_em = merge(RNS_EM_POTENTIAL, RNS_EM_FARADAY, COCAL_ID_read_magnetic_potential == 1)
    rns_em0 = NMETRIC_WL + 1
    nextra = merge(4, 3, rns_em == RNS_EM_POTENTIAL)
  end if
  call COCAL_ID_read_one_patch_metric(dir, fmt, wl, nextra)

  call COCAL_ID_domain_init(cocal_matter, 1, merge(MAT_U + 3, MAT_OMEGA, mrns), .false.)
  associate(p => cocal_matter%patch(1))
    call COCAL_ID_patch_set_grids(p, rg(0:nrf), thg(0:ntg), phig(0:npg))
    p%surface_fitted = .true.
  end associate

  if (leader) then
    associate(raw => cocal_metric%patch(1)%raw)
      if (rns_em == RNS_EM_POTENTIAL) then
        va(0:,0:,0:) => raw(rns_em0,:,:,:);     vaxd(0:,0:,0:) => raw(rns_em0+1,:,:,:)
        vayd(0:,0:,0:) => raw(rns_em0+2,:,:,:); vazd(0:,0:,0:) => raw(rns_em0+3,:,:,:)
        call COCAL_ID_IO_input_grav_export_Ai(dir//"/rnsEMF_3D.las", fmt, va, vaxd, vayd, vazd)
      else if (rns_em == RNS_EM_FARADAY) then
        fxyd(0:,0:,0:) => raw(rns_em0,:,:,:); fxzd(0:,0:,0:) => raw(rns_em0+1,:,:,:); fyzd(0:,0:,0:) => raw(rns_em0+2,:,:,:)
        ! F_ti, needed only to read the file.
        allocate(fxd, fyd, fzd, mold=fxyd)
        call COCAL_ID_IO_input_grav_export_Faraday(dir//"/rnsEMF_faraday_3D.las", fmt, fxd, fyd, fzd, fxyd, fxzd, fyzd)
        deallocate(fxd, fyd, fzd)
      end if
    end associate

    call COCAL_ID_domain_alloc_raw(cocal_matter, 1)
    associate(raw => cocal_matter%patch(1)%raw)
      emd(0:,0:,0:) => raw(MAT_EMD,:,:,:); omef(0:,0:,0:) => raw(MAT_OMEGA,:,:,:)
      rs(0:,0:) => cocal_matter%patch(1)%raw_surface
      call COCAL_ID_IO_input_CF_star_export(dir//"/rnsflu_3D.las", fmt, emd, rs, omef, ome, ber, radi)
      if (mrns) then
        utf(0:,0:,0:) => raw(MAT_U,:,:,:);   uxf(0:,0:,0:) => raw(MAT_U+1,:,:,:)
        uyf(0:,0:,0:) => raw(MAT_U+2,:,:,:); uzf(0:,0:,0:) => raw(MAT_U+3,:,:,:)
        call COCAL_ID_IO_input_star4ve_export(dir//"/rns4ve_3D.las", fmt, utf, uxf, uyf, uzf)
      end if
    end associate
  end if
  call COCAL_ID_node_bcast(radi)
  if (.not. radi > 0.0d0) call CCTK_WARN(CCTK_WARN_ABORT, "COCAL_ID RNS: radi must be positive.")
  cocal_radi = radi

  call COCAL_ID_domain_prepare(cocal_metric, ok, 0.0d0)
  if (ok) call COCAL_ID_domain_prepare(cocal_matter, ok, 0.0d0)
  if (.not. ok) call CCTK_WARN(CCTK_WARN_ABORT, "COCAL_ID RNS: the exported grids must have increasing radii and uniform angles.")
end subroutine COCAL_ID_read_rns
