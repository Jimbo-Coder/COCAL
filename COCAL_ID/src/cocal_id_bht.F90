! Black hole with a torus: the one-patch waveless metric of RNS about an
! excised black hole, with the torus on the gravitational grid and the hole
! filled by cocal_id_holes.f90.

! Collective over the node.
subroutine COCAL_ID_read_bht(dir)
  ! COCAL's parameters, grids and the routines that build them.
  use COCAL_ID_grid_parameter, eps_cocal => eps
  use COCAL_ID_coordinate_grav_r
  use COCAL_ID_coordinate_grav_theta
  use COCAL_ID_coordinate_grav_phi
  use COCAL_ID_trigonometry_grav_theta
  use COCAL_ID_trigonometry_grav_phi
  use COCAL_ID_def_bh_parameter, only : mass_bh, spin_bh
  use COCAL_ID_interface_IO_input_matter_BHT_export
  use COCAL_ID_node_shared
  use COCAL_ID_state
  use COCAL_ID_domain
  use COCAL_ID_holes
  implicit none
  DECLARE_CCTK_PARAMETERS
  DECLARE_CCTK_FUNCTIONS
  character(*), intent(in) :: dir
  ! COCAL's routines take the directory as character*400.
  character(400) :: cocal_dir
  character(*), parameter :: fmt = "23.15"
  real(8), pointer :: emdg(:,:,:), omeg(:,:,:)
  real(8) :: ome, ber, radi, r_eh
  character(400) :: msg
  logical :: leader, ok

  leader = COCAL_ID_node_leader()
  cocal_dir = dir
  call read_parameter_cactus(cocal_dir)
  call COCAL_ID_read_bht_parameter_cactus(cocal_dir)
  if (mass_bh <= 0.0d0) call CCTK_WARN(CCTK_WARN_ABORT, "COCAL_ID BHT: mass_bh must be positive.")
  if (abs(spin_bh) > mass_bh) call CCTK_WARN(CCTK_WARN_ABORT, "COCAL_ID BHT: abs(spin_bh) exceeds mass_bh.")
  call COCAL_ID_calc_bht_excision_radius
  if (cocal_tabulated_eos) then
    call COCAL_ID_teos_initialize_cactus(cocal_dir)
  else
    call COCAL_ID_peos_initialize_cactus(cocal_dir)
  end if
  call grid_r_bht('eBH')
  call grid_theta
  call trig_grav_theta
  call grid_phi
  call allocate_trig_grav_mphi
  call trig_grav_phi

  call COCAL_ID_read_one_patch_metric(dir, fmt, .true., 0)
  call COCAL_ID_domain_init(cocal_matter, 1, MAT_OMEGA, .false.)
  associate(p => cocal_matter%patch(1))
    call COCAL_ID_patch_set_grids(p, rg(0:nrg), thg(0:ntg), phig(0:npg))
    p%support_index = MAT_EMD
    p%support_floor = EMD_VACUUM
  end associate
  if (leader) then
    call COCAL_ID_domain_alloc_raw(cocal_matter, 1)
    associate(raw => cocal_matter%patch(1)%raw)
      emdg(0:,0:,0:) => raw(MAT_EMD,:,:,:); omeg(0:,0:,0:) => raw(MAT_OMEGA,:,:,:)
      call COCAL_ID_IO_input_matter_BHT_export(dir//"/rnsflu_3D.las", fmt, emdg, omeg, ome, ber, radi)
    end associate
  end if
  call COCAL_ID_node_bcast(radi)
  if (.not. radi > 0.0d0) call CCTK_WARN(CCTK_WARN_ABORT, "COCAL_ID BHT: radi must be positive.")
  cocal_radi = radi

  call COCAL_ID_domain_prepare(cocal_metric, ok, 0.0d0)
  if (ok) call COCAL_ID_domain_prepare(cocal_matter, ok, 0.0d0)
  if (.not. ok) call CCTK_WARN(CCTK_WARN_ABORT, "COCAL_ID BHT: the exported grids must have increasing radii and uniform angles.")

  ! The hole is filled inside 0.99 of the Kerr horizon radius.
  cocal_fill_holes = .true.
  r_eh = mass_bh + sqrt(max(mass_bh**2 - spin_bh**2, 0.0d0))
  call COCAL_ID_holes_add(1, [0.0d0, 0.0d0, 0.0d0], [0.0d0, 0.0d0, 0.0d0], .false., 0.99d0*r_eh)
  call COCAL_ID_holes_extrapolate_centre(1, ok)
  if (.not. ok) call CCTK_WARN(CCTK_WARN_ABORT, "COCAL_ID BHT: cannot sample the fill centre.")
  associate(h => cocal_holes(1))
    if (h%rh <= rg(0)) call CCTK_WARN(CCTK_WARN_ABORT, "COCAL_ID BHT: the fill radius must exceed the excision radius.")
    if (h%dr <= 0.0d0) call CCTK_WARN(CCTK_WARN_ABORT, "COCAL_ID BHT: the fill sample spacing must be positive.")
    if (COCAL_ID_verbose == 1) then
      write(msg, '(a,1p,4e20.12)') "COCAL_ID BHT radii (COCAL units) excision, fill, horizon, ring: ", &
        rg(0), h%rh, r_eh, abs(spin_bh)
      call CCTK_INFO(trim(msg))
    end if
  end associate


end subroutine COCAL_ID_read_bht
