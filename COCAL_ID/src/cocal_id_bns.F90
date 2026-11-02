! Binary neutron stars: three conformally flat patches, one per star and the
! outer patch, with each star's fluid on its own surface-fitted grid.

module COCAL_ID_bns
  use COCAL_ID_state
  use COCAL_ID_domain
  implicit none
  private
  public :: bns_star, bns_stars, bns_id_type, bns_dis, bns_xunit, bns_ecc_velx
  public :: BNS_EMD, BNS_DVEP, BNS_W, COCAL_ID_bns_matter

  ! Matter fields: Emden function, gradient of the velocity potential
  ! (IR, SP) and spin vector (SP).
  integer, parameter :: BNS_EMD = 1, BNS_DVEP = 2, BNS_W = 5

  ! Each star's fluid and orbital constants.
  type :: bns_star
    real(8) :: ome = 0, ber = 0, radi = 0, xcm = 0, zcm = 0, confpow = 0
    real(8) :: orbital_cos = 1, orbital_sin = 0
  end type bns_star

  type(bns_star), save :: bns_stars(2)
  character(2), save :: bns_id_type = ''
  ! Distance of each star patch centre from the binary's centre.
  real(8), save :: bns_dis = 0
  ! Eccentricity-reduction velocity (param.ccl).
  logical, save :: bns_xunit = .false.
  real(8), save :: bns_ecc_velx = 0

contains

  ! Fluid at COCAL point (xc,yc,zc) with blended metric fields q there.
  ! Matter belongs to one star, chosen by the midplane.
  subroutine COCAL_ID_bns_matter(xc, yc, zc, q, s)
    real(8), intent(in) :: xc, yc, zc, q(:)
    real(8), intent(inout) :: s(NSTATE)
    real(8) :: f(cocal_matter%nvars), xs, ys, orbital_x, beta(3), dvep(3), w(3)
    real(8) :: alph, psi4, psicp, lam, wdvep, w2, huta, v(3)
    integer :: star
    logical :: ok

    if (xc <= 64*epsilon(1.0d0)*max(1.0d0, abs(bns_dis))) then
      star = 1
      xs = xc + bns_dis
      ys = yc
      orbital_x = xc - bns_stars(1)%xcm
    else
      star = 2
      xs = bns_dis - xc
      ys = -yc
      orbital_x = xc + bns_stars(2)%xcm
    end if
    call COCAL_ID_domain_patch_eval(cocal_matter, star, xs, ys, zc, f, ok)
    if (.not. ok) return

    associate(p => bns_stars(star))
      ! beta^i + Omega x r about the orbital axis (-orbital_sin, 0, orbital_cos).
      beta = q(M_BETA:M_BETA + 2) + p%ome*[-p%orbital_cos*yc, p%orbital_cos*orbital_x + p%orbital_sin*(zc - p%zcm), &
                                           -p%orbital_sin*yc]
      if (bns_xunit) then
        beta(1) = beta(1) + merge(-1.0d0, 1.0d0, star == 1)*bns_ecc_velx*bns_dis
      else
        beta = beta + bns_ecc_velx*[orbital_x, yc, zc]
      end if

      alph = q(M_ALP)
      psi4 = q(M_PSI)**4
      select case (bns_id_type)
      case ("CO")
        v = beta/alph
      case ("IR")
        dvep = f(BNS_DVEP:BNS_DVEP + 2)
        huta = (p%ber + dot_product(beta, dvep))/alph
        v = dvep/psi4/huta
      case ("SP")
        dvep = f(BNS_DVEP:BNS_DVEP + 2)
        psicp = q(M_PSI)**p%confpow
        w = psicp*f(BNS_W:BNS_W + 2)
        lam = p%ber + dot_product(beta, dvep)
        wdvep = dot_product(w, dvep)
        w2 = psi4*dot_product(w, w)
        huta = (lam + sqrt(lam*lam + 4.0d0*alph*alph*(wdvep + w2)))/(2.0d0*alph)
        v = (dvep/psi4 + w)/huta
      end select
    end associate
    call COCAL_ID_set_fluid(f(BNS_EMD), v, s)
  end subroutine COCAL_ID_bns_matter

end module COCAL_ID_bns

! Collective over the node.
subroutine COCAL_ID_read_bns(dir)
  use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
  use COCAL_ID_grid_parameter_binary_excision
  use COCAL_ID_grid_parameter, eps_cocal => eps
  use COCAL_ID_coordinate_grav_r
  use COCAL_ID_coordinate_grav_theta
  use COCAL_ID_coordinate_grav_phi
  use COCAL_ID_def_binary_parameter, only : dis
  use COCAL_ID_trigonometry_grav_phi, only : allocate_trig_grav_mphi
  use COCAL_ID_interface_IO_input_CF_grav_export
  use COCAL_ID_interface_IO_input_CF_flir_export
  use COCAL_ID_interface_IO_input_CF_flco_export
  use COCAL_ID_interface_IO_input_CF_flsp_export
  use COCAL_ID_interface_IO_input_CF_surf_export
  use COCAL_ID_interface_IO_input_grav_export_Kij
  use COCAL_ID_interface_IO_input_gradvep_export
  use COCAL_ID_node_shared
  use COCAL_ID_state
  use COCAL_ID_domain
  use COCAL_ID_bns
  implicit none
  DECLARE_CCTK_PARAMETERS
  DECLARE_CCTK_FUNCTIONS
  character(*), intent(in) :: dir
  ! COCAL's routines take the directory as character*400.
  character(400) :: cocal_dir
  ! Patch 2 is rotated by pi about z: x/y vectors and xz/yz tensors change sign.
  integer, parameter :: metric_flips(4) = [M_BETA, M_BETA + 1, M_K + 2, M_K + 4]
  character(*), parameter :: fmt = "20.12"
  character(400) :: msg
  character(2) :: id_type_p2
  logical :: leader, ok
  integer :: impt, nmatter, igrid
  real(8) :: radi_tol, star_extent, rs_max(2), dis_cm

  leader = COCAL_ID_node_leader()
  cocal_dir = dir
  call read_id_type(dir//"/rnspar_mpt1.dat", bns_id_type)
  call read_id_type(dir//"/rnspar_mpt2.dat", id_type_p2)
  if (id_type_p2 /= bns_id_type) &
    call CCTK_WARN(CCTK_WARN_ABORT, "COCAL_ID BNS: both star patches must use the same CO/IR/SP data format.")
  select case (bns_id_type)
  case ("CO")
    nmatter = 1
  case ("IR")
    nmatter = 4
  case ("SP")
    nmatter = 7
  case default
    call CCTK_WARN(CCTK_WARN_ABORT, "COCAL_ID BNS: unsupported fluid data type.")
  end select
  bns_xunit = COCAL_ID_bns_xunit == 1
  bns_ecc_velx = COCAL_ID_ecc_cor_velx
  igrid = 4     ! Legacy grids (rg(nrf) = 1).

  ! COCAL's parameters for each patch.
  call COCAL_ID_allocate_grid_parameter_mpt
  call COCAL_ID_allocate_grid_parameter_binary_excision_mpt
  radi_tol = 1.0d-12
  do impt = 1, 3
    call COCAL_ID_read_parameter_mpt_cactus(impt, cocal_dir)
    if (impt <= 2) radi_tol = max(radi_tol, abs(eps_cocal))
    indata_type = '3D'
    if (impt <= 2) call COCAL_ID_read_surf_parameter_mpt_cactus(impt, cocal_dir)
    call COCAL_ID_copy_grid_parameter_to_mpt(impt)
    call COCAL_ID_read_parameter_binary_excision_mpt_cactus(impt, cocal_dir)
    call COCAL_ID_copy_grid_parameter_binary_excision_to_mpt(impt)
    if (cocal_tabulated_eos) then
      if (impt <= 2) call COCAL_ID_teos_initialize_mpt_cactus(impt, cocal_dir)
      call COCAL_ID_copy_def_teos_parameter_to_mpt(impt)
    else
      if (impt <= 2) call COCAL_ID_peos_initialize_mpt_cactus(impt, cocal_dir)
      call COCAL_ID_copy_def_peos_parameter_to_mpt(impt)
    end if
  end do

  call COCAL_ID_set_allocate_size_mpt
  call allocate_trig_grav_mphi
  call COCAL_ID_domain_init(cocal_metric, 3, NMETRIC_CF, .true.)
  call COCAL_ID_domain_init(cocal_matter, 2, nmatter, .false.)
  do impt = 1, 3
    call COCAL_ID_copy_grid_parameter_from_mpt(impt)
    call COCAL_ID_copy_grid_parameter_binary_excision_from_mpt(impt)
    if (cocal_tabulated_eos) then
      call COCAL_ID_copy_def_teos_parameter_from_mpt(impt)
    else
      call COCAL_ID_copy_def_peos_parameter_from_mpt(impt)
    end if
    call COCAL_ID_coordinate_patch_kit_grav_grid_coc2cac_mpt(igrid, impt)
    call calc_parameter_binary_excision
    call COCAL_ID_copy_grid_parameter_to_mpt(impt)
    call COCAL_ID_copy_grid_parameter_binary_excision_to_mpt(impt)
    if (impt == 1) then
      cocal_metric%rr3 = 0.7d0*(rgout - rgmid)
      dis_cm = dis
    else if (impt == 2 .and. abs(dis - dis_cm) > 1.0d-8*max(1.0d0, abs(dis_cm))) then
      write(msg, '(a,es20.12,a,es20.12)') &
        'COCAL_ID BNS: patches disagree on the distance from their shared midpoint: mpt1=', dis_cm, '; mpt2=', dis
      call CCTK_WARN(CCTK_WARN_ABORT, trim(msg))
    end if
    call read_patch(impt)
  end do
  bns_dis = dis_cm

  ! The other ranks of the node take the constants from the leader.
  do impt = 1, 2
    associate(p => bns_stars(impt))
      call COCAL_ID_node_bcast(p%ome)
      call COCAL_ID_node_bcast(p%ber)
      call COCAL_ID_node_bcast(p%radi)
      call COCAL_ID_node_bcast(p%xcm)
      call COCAL_ID_node_bcast(p%zcm)
      call COCAL_ID_node_bcast(p%confpow)
      call COCAL_ID_node_bcast(rs_max(impt))
      if (.not. all(ieee_is_finite([p%radi, p%ome, p%ber, p%xcm, p%zcm]))) &
        call CCTK_WARN(CCTK_WARN_ABORT, "COCAL_ID BNS: non-finite fluid constants.")
      if (p%radi <= 0.0d0) call CCTK_WARN(CCTK_WARN_ABORT, "COCAL_ID BNS: radi must be positive.")
      call set_orbital_direction(p)
    end associate
  end do
  ! radi is a shared length unit, not the individual stellar radius. Allow iteration/print noise.
  if (abs(bns_stars(1)%radi - bns_stars(2)%radi) > radi_tol*max(bns_stars(1)%radi, bns_stars(2)%radi)) &
    call CCTK_WARN(CCTK_WARN_ABORT, "COCAL_ID BNS: unequal radi normalizations are not supported by this coordinate map.")
  cocal_radi = bns_stars(1)%radi

  star_extent = max(rs_max(1)*cocal_matter%patch(1)%rg(cocal_matter%patch(1)%nr), &
                    rs_max(2)*cocal_matter%patch(2)%rg(cocal_matter%patch(2)%nr))
  call COCAL_ID_domain_prepare(cocal_metric, ok, star_extent)
  if (ok) call COCAL_ID_domain_prepare(cocal_matter, ok, 0.0d0)
  if (.not. ok) call CCTK_WARN(CCTK_WARN_ABORT, "COCAL_ID BNS: the exported grids must have increasing radii and uniform angles.")
  if (COCAL_ID_verbose == 1) then
    do impt = 1, 2
      associate(p => bns_stars(impt))
        write(msg, '(a,i1,a,5es20.12)') 'COCAL_ID BNS star ', impt, ' ome, ber, radi, xcm, zcm: ', &
          p%ome, p%ber, p%radi, p%xcm, p%zcm
      end associate
      call CCTK_INFO(trim(msg))
    end do
  end if

contains

  ! Patch impt's grids from COCAL's grid modules, then, on the node leader,
  ! its fields: the metric, and for a star its fluid and surface.
  subroutine read_patch(impt)
    integer, intent(in) :: impt
    real(8), pointer :: psi(:,:,:), alph(:,:,:), bvxd(:,:,:), bvyd(:,:,:), bvzd(:,:,:)
    real(8), pointer :: axx(:,:,:), axy(:,:,:), axz(:,:,:), ayy(:,:,:), ayz(:,:,:), azz(:,:,:)
    real(8), pointer :: emd(:,:,:), vep(:,:,:), vepxf(:,:,:), vepyf(:,:,:), vepzf(:,:,:)
    real(8), pointer :: wxspf(:,:,:), wyspf(:,:,:), wzspf(:,:,:), rs(:,:)
    character(1) :: patch
    character(400) :: base

    write(patch, '(i1)') impt
    associate(p => cocal_metric%patch(impt))
      call COCAL_ID_patch_set_grids(p, rg(0:nrg), thg(0:ntg), phig(0:npg))
      p%dis = dis
      if (impt <= 2) p%ex_radius = p%rg(ex_nrg)
      p%positive_indices = [M_PSI, M_ALP]
      if (impt == 2) p%flip_indices = metric_flips
    end associate
    if (impt <= 2) then
      associate(p => cocal_matter%patch(impt))
        call COCAL_ID_patch_set_grids(p, rg(0:nrf), thg(0:ntg), phig(0:npg))
        p%surface_fitted = .true.
        if (impt == 2 .and. nmatter > 1) p%flip_indices = [BNS_DVEP, BNS_DVEP + 1]
        if (impt == 2 .and. nmatter > 4) p%flip_indices = [BNS_DVEP, BNS_DVEP + 1, BNS_W, BNS_W + 1]
      end associate
      if (ntf /= ntg .or. npf /= npg) &
        call CCTK_WARN(CCTK_WARN_ABORT, "COCAL_ID BNS: fluid and gravitational angular grids must match.")
    end if
    if (.not. leader) return

    ! Each name below is a view of one field of raw, which COCAL's reading
    ! routines fill in place and the spline fit then consumes.
    call COCAL_ID_domain_alloc_raw(cocal_metric, impt)
    associate(raw => cocal_metric%patch(impt)%raw)
      psi(0:,0:,0:)  => raw(M_PSI,:,:,:);    alph(0:,0:,0:) => raw(M_ALP,:,:,:)
      bvxd(0:,0:,0:) => raw(M_BETA,:,:,:);   bvyd(0:,0:,0:) => raw(M_BETA+1,:,:,:); bvzd(0:,0:,0:) => raw(M_BETA+2,:,:,:)
      axx(0:,0:,0:)  => raw(M_K,:,:,:);      axy(0:,0:,0:)  => raw(M_K+1,:,:,:);    axz(0:,0:,0:)  => raw(M_K+2,:,:,:)
      ayy(0:,0:,0:)  => raw(M_K+3,:,:,:);    ayz(0:,0:,0:)  => raw(M_K+4,:,:,:);    azz(0:,0:,0:)  => raw(M_K+5,:,:,:)
    end associate
    base = dir//"/bnsgra_3D_mpt"//patch//".las"
    call COCAL_ID_IO_input_CF_grav_export(trim(base), fmt, psi, alph, bvxd, bvyd, bvzd)
    base = dir//"/bnsgra_Kij_3D_mpt"//patch//".las"
    call COCAL_ID_IO_input_grav_export_Kij(trim(base), axx, axy, axz, ayy, ayz, azz)
    if (impt == 3) return

    call COCAL_ID_domain_alloc_raw(cocal_matter, impt)
    associate(raw => cocal_matter%patch(impt)%raw, p => bns_stars(impt))
      emd(0:,0:,0:) => raw(BNS_EMD,:,:,:)
      rs(0:,0:) => cocal_matter%patch(impt)%raw_surface
      base = dir//"/bnsflu_3D_mpt"//patch//".las"
      select case (bns_id_type)
      case ("CO")
        call COCAL_ID_IO_input_CF_flco_export(trim(base), fmt, emd, p%ome, p%ber, p%radi, p%xcm, p%zcm)
      case ("IR", "SP")
        ! The potential itself is not needed, only its gradient.
        allocate(vep, mold=emd)
        if (bns_id_type == "IR") then
          call COCAL_ID_IO_input_CF_flir_export(trim(base), fmt, emd, vep, p%ome, p%ber, p%radi, p%xcm, p%zcm)
        else
          wxspf(0:,0:,0:) => raw(BNS_W,:,:,:); wyspf(0:,0:,0:) => raw(BNS_W+1,:,:,:); wzspf(0:,0:,0:) => raw(BNS_W+2,:,:,:)
          call COCAL_ID_IO_input_CF_flsp_export(trim(base), fmt, emd, vep, wxspf, wyspf, wzspf, &
                                                p%ome, p%ber, p%radi, p%confpow, p%xcm, p%zcm)
        end if
        deallocate(vep)
        vepxf(0:,0:,0:) => raw(BNS_DVEP,:,:,:); vepyf(0:,0:,0:) => raw(BNS_DVEP+1,:,:,:); vepzf(0:,0:,0:) => raw(BNS_DVEP+2,:,:,:)
        base = dir//"/bnsdvep_3D_mpt"//patch//".las"
        call COCAL_ID_IO_input_gradvep_export(trim(base), fmt, vepxf, vepyf, vepzf)
      end select
      base = dir//"/bnssur_3D_mpt"//patch//".las"
      call COCAL_ID_IO_input_CF_surf_export(trim(base), fmt, rs)
      rs_max(impt) = maxval(rs)
    end associate
  end subroutine read_patch

  subroutine set_orbital_direction(p)
    type(bns_star), intent(inout) :: p
    real(8) :: norm

    p%orbital_cos = 1.0d0
    p%orbital_sin = 0.0d0
    if (p%zcm == 0.0d0) return
    norm = norm2([dis_cm + p%xcm, p%zcm])
    if (.not. ieee_is_finite(norm) .or. norm <= 0.0d0) &
      call CCTK_WARN(CCTK_WARN_ABORT, "COCAL_ID BNS: invalid orbital-axis normalization.")
    p%orbital_cos = (dis_cm + p%xcm)/norm
    p%orbital_sin = p%zcm/norm
  end subroutine set_orbital_direction

  subroutine read_id_type(filename, id_type)
    character(*), intent(in) :: filename
    character(2), intent(out) :: id_type
    integer :: u, nrg_, nrf_, nrf_deform, nrgin
    character(2) :: eq_point

    open(newunit=u, file=filename, status='old', action='read')
    read(u, '(4i5)') nrg_
    read(u, '(4i5)') nrf_
    read(u, '(2i5,2(3x,a2))') nrf_deform, nrgin, id_type, eq_point
    close(u)
  end subroutine read_id_type

end subroutine COCAL_ID_read_bns
