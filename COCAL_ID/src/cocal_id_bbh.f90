! Binary black holes (BBH) and binary black holes with a disk (BBHT): three
! conformally flat patches, read from their own export format. Coordinates
! are read from bbhgrids, never regenerated from rnspar (which can describe
! another mesh). Field order: psi, alpha, beta^x, beta^y, beta^z,
! Axx,Axy,Axz,Ayy,Ayz,Azz, as cocal_id_state.f90 orders the metric. The saved
! curvature is K_ij / psi**4 in the patch Cartesian basis. psi is fitted as
! r*psi on the hole patches; each hole is filled from its own patch by
! cocal_id_holes.f90.
!
! BBHT is Mai's BBHT_BHT_FLUID_V1: a disk on the outer patch's grid, in units
! G=c=M_sun=1, with rest-mass density rho, Omega = dphi/dt and a single-Gamma
! polytrope, and no magnetic fields or potentials.
!
! The readers return ok=.false. with a message in bbh_error; the Cactus
! wrapper aborts with it.
module COCAL_ID_bbh
  use COCAL_ID_phys_constant, only : pi
  use COCAL_ID_state
  use COCAL_ID_domain
  use COCAL_ID_holes
  use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
  implicit none
  private
  public :: bbh_error, COCAL_ID_read_bbh, COCAL_ID_bbht_matter

  ! Patch 2 is rotated by pi about z: x/y vectors and xz/yz tensors change sign.
  integer, parameter :: BBH_ROTATED(4) = [M_BETA, M_BETA + 1, M_K + 2, M_K + 4]
  ! Disk fields on the outer patch's grid.
  integer, parameter :: BBHT_RHO = 1, BBHT_OMEGA = 2
  character(1024), save :: bbh_error = ''
  ! Disk EOS and inner edge (the zero-density surface, inside a cell).
  real(8), save :: bbht_k = 0, bbht_gamma = 0, bbht_rin = 0

contains

  subroutine fail(message)
    character(*), intent(in) :: message
    bbh_error = message
  end subroutine fail

  logical function read_ok(ios, path)
    integer, intent(in) :: ios
    character(*), intent(in) :: path
    read_ok = ios == 0
    if (.not. read_ok) write(bbh_error, '(a,a,a,i0,a)') 'COCAL_ID BBH: cannot read ', trim(path), &
      ' (iostat=', ios, '): invalid or incomplete export'
  end function read_ok

  ! Collective over the node.
  subroutine COCAL_ID_read_bbh(dir_path, ok, fill_interior, with_disk)
    use COCAL_ID_node_shared
    character(*), intent(in) :: dir_path
    logical, intent(out) :: ok
    logical, intent(in) :: fill_interior, with_disk
    character(1024) :: path
    character(1) :: npatch
    integer :: m, u, ios, nr, nt, np, ir, it, ip, ex_nr, ex_dis
    real(8) :: status, constants(3, 3)
    real(8), allocatable :: rg(:), thg(:), phig(:)
    logical :: disk_exists, leader
    ok = .false.
    if (with_disk .and. cocal_tabulated_eos) then
      call fail('COCAL_ID BBHT: the disk has a single-Gamma polytrope, not a tabulated EOS'); return
    end if
    inquire(file=trim(dir_path)//'/disk_state.las', exist=disk_exists)
    if (disk_exists .neqv. with_disk) then
      if (disk_exists) call fail('COCAL_ID BBH: these are BBHT data (disk_state.las); use CocalBBHT')
      if (.not. disk_exists) call fail('COCAL_ID BBHT: no disk_state.las')
      return
    end if
    ! Every rank reads the (small) grids and excision indices; the fields are
    ! read once per node, below.
    call COCAL_ID_domain_init(cocal_metric, 3, NMETRIC_CF, .true.)
    do m = 1, 3
      write(npatch, '(i1)') m
      path = trim(dir_path)//'/bbhgrids_3D_mpt'//npatch//'.las'
      open(newunit=u, file=trim(path), status='old', action='read', iostat=ios)
      if (.not. read_ok(ios, path)) return
      read(u, *, iostat=ios) nr, nt, np
      if (.not. read_ok(ios, path)) return
      if (nr < 3 .or. nt < 4 .or. np < 4 .or. mod(np, 4) /= 0) then
        call fail('COCAL_ID BBH: need four radial nodes and full angular grids'); return
      end if
      allocate(rg(0:nr), thg(0:nt), phig(0:np))
      do ir = 0, nr
        read(u, *, iostat=ios) rg(ir)
        if (.not. read_ok(ios, path)) return
      end do
      do it = 0, nt
        read(u, *, iostat=ios) thg(it)
        if (.not. read_ok(ios, path)) return
      end do
      do ip = 0, np
        read(u, *, iostat=ios) phig(ip)
        if (.not. read_ok(ios, path)) return
      end do
      close(u)
      associate(p => cocal_metric%patch(m))
        call COCAL_ID_patch_set_grids(p, rg, thg, phig)
        deallocate(rg, thg, phig)
        if (.not. all(ieee_is_finite(p%rg)) .or. .not. all(ieee_is_finite(p%thg)) .or. &
            .not. all(ieee_is_finite(p%phig))) then
          call fail('COCAL_ID BBH: nonfinite grid'); return
        end if
        if (any(p%rg(1:) <= p%rg(:nr-1)) .or. any(p%thg(1:) <= p%thg(:nt-1)) .or. &
            any(p%phig(1:) <= p%phig(:np-1))) then
          call fail('COCAL_ID BBH: nonmonotone grid'); return
        end if
        if (p%rg(0) < 0 .or. abs(p%thg(0)) > 1d-10 .or. abs(p%thg(nt)-pi) > 1d-10 .or. &
            abs(p%phig(0)) > 1d-10 .or. abs(p%phig(np)-2*pi) > 1d-10) then
          call fail('COCAL_ID BBH: unsupported angular domain'); return
        end if
        path = trim(dir_path)//'/bin_ex_par_mpt'//npatch//'.dat'
        open(newunit=u, file=trim(path), status='old', action='read', iostat=ios)
        if (.not. read_ok(ios, path)) return
        read(u, *, iostat=ios) ex_nr, ex_dis
        if (.not. read_ok(ios, path)) return
        close(u)
        if (ex_nr < 0 .or. ex_nr > nr .or. ex_dis < 0 .or. ex_nr+ex_dis > nr) then
          call fail('COCAL_ID BBH: invalid companion-excision indices'); return
        end if
        if (m <= 2 .and. ex_nr > 0) p%ex_radius = p%rg(ex_nr)
        p%positive_indices = [M_PSI, M_ALP]
        if (m <= 2) p%r_scaled_indices = [M_PSI]
        if (m == 2) p%flip_indices = BBH_ROTATED
      end associate
    end do
    leader = COCAL_ID_node_leader()
    ok = .true.
    constants = 0
    if (leader) then
      do m = 1, 3
        call COCAL_ID_domain_alloc_raw(cocal_metric, m)
        call read_fields(dir_path, m, constants(:, m), ok)
        if (.not. ok) exit
      end do
    end if
    ! The other ranks of the node take the leader's status and constants.
    status = merge(1d0, 0d0, ok)
    call COCAL_ID_node_bcast(status)
    do m = 1, 3
      do ir = 1, 3
        call COCAL_ID_node_bcast(constants(ir, m))
      end do
    end do
    ok = status == 1
    if (.not. ok) then
      if (.not. leader) call fail('COCAL_ID BBH: the node leader could not read the export')
      return
    end if
    ! constants(:,m) = dis, rgin, sepa of patch m.
    ok = .false.
    if (.not. all(ieee_is_finite(constants))) then
      call fail('COCAL_ID BBH: nonfinite constants'); return
    end if
    do m = 1, 2
      if (constants(2, m) <= 0 .or. constants(1, m) <= 0 .or. &
          abs(constants(2, m)-cocal_metric%patch(m)%rg(0)) > 1d-10) then
        call fail('COCAL_ID BBH: inconsistent patch geometry'); return
      end if
      cocal_metric%patch(m)%dis = constants(1, m)
    end do
    if (abs(constants(3, 1)-constants(3, 2)) > 1d-10 .or. &
        abs(constants(1, 1)+constants(1, 2)-constants(3, 1)) > 1d-10) then
      call fail('COCAL_ID BBH: inconsistent centers'); return
    end if
    if (constants(2, 1)+constants(2, 2) >= constants(3, 1)) then
      call fail('COCAL_ID BBH: overlapping BH excisions'); return
    end if
    ! Same interior/outer patch split used by COCAL_ID's BNS reader.
    associate(p1 => cocal_metric%patch(1), p2 => cocal_metric%patch(2))
      cocal_metric%rr3 = 0.7d0*min(p1%rg(p1%nr)-p1%dis, p2%rg(p2%nr)-p2%dis)
    end associate
    call COCAL_ID_domain_prepare(cocal_metric, ok, 0.0d0)
    if (.not. ok) then
      call fail('COCAL_ID BBH: the exported grids must have increasing radii and uniform angles'); return
    end if
    ! Each hole is the excision of its own patch, filled from that patch.
    cocal_fill_holes = fill_interior
    call COCAL_ID_holes_add(1, [-constants(1, 1), 0d0, 0d0], [0d0, 0d0, 0d0], .false., constants(2, 1))
    call COCAL_ID_holes_add(2, [constants(1, 2), 0d0, 0d0], [0d0, 0d0, 0d0], .true., constants(2, 2))
    if (cocal_fill_holes) then
      do m = 1, 2
        call COCAL_ID_holes_extrapolate_centre(m, ok)
        if (.not. ok) then
          call fail('COCAL_ID BBH: cannot sample automatic fill center'); return
        end if
      end do
    end if
    if (with_disk) call read_disk_domain(dir_path, ok)
  end subroutine COCAL_ID_read_bbh

  ! Node leader only: the metric and curvature of patch m, and its dis, rgin, sepa.
  subroutine read_fields(dir_path, m, geometry, ok)
    character(*), intent(in) :: dir_path
    integer, intent(in) :: m
    real(8), intent(out) :: geometry(3)
    logical, intent(out) :: ok
    character(1024) :: path, line
    character(1) :: npatch
    character(32) :: fmt
    integer :: u, ios, nr, nt, np, ir, it, ip, width
    real(8) :: constants(6), dist2, sepa
    ok = .false.
    write(npatch, '(i1)') m
    associate(p => cocal_metric%patch(m))
      path = trim(dir_path)//'/bbhgra_3D_mpt'//npatch//'.las'
      open(newunit=u, file=trim(path), status='old', action='read', iostat=ios)
      if (.not. read_ok(ios, path)) return
      read(u, *, iostat=ios) nr, nt, np
      if (.not. read_ok(ios, path)) return
      if (nr /= p%nr .or. nt /= p%nt .or. np /= p%np) then
        call fail('COCAL_ID BBH: field/grid shape mismatch'); return
      end if
      read(u, '(a)', iostat=ios) line
      if (.not. read_ok(ios, path)) return
      width = 20
      if (len_trim(line) > 100) width = 23
      write(fmt, '("(5e",i0,".0)")') width
      backspace(u)
      do ip = 0, np
        do it = 0, nt
          do ir = 0, nr
            read(u, fmt, iostat=ios) p%raw(M_PSI:M_BETA+2, ir, it, ip)
            if (.not. read_ok(ios, path)) return
          end do
        end do
      end do
      read(u, '(a)', iostat=ios) line
      if (.not. read_ok(ios, path)) return
      ! Historical exports have ome,0,0,dis,rgin,sepa; newer exports have
      ! ome,dis,rgin,sepa. Preserve both original fixed-width formats.
      if (len_trim(line) > 5*width) then
        write(fmt, '("(6e",i0,".0)")') width
        read(line, fmt, iostat=ios) constants
        geometry = constants(4:6)
      else
        write(fmt, '("(4e",i0,".0)")') width
        read(line, fmt, iostat=ios) constants(1:4)
        geometry = constants(2:4)
      end if
      if (.not. read_ok(ios, path)) return
      close(u)

      path = trim(dir_path)//'/bnsgra_Kij_3D_mpt'//npatch//'.las'
      open(newunit=u, file=trim(path), status='old', action='read', iostat=ios)
      if (.not. read_ok(ios, path)) return
      read(u, *, iostat=ios) nr, nt, np
      if (.not. read_ok(ios, path)) return
      if (nr /= p%nr .or. nt /= p%nt .or. np /= p%np) then
        call fail('COCAL_ID BBH: curvature/grid shape mismatch'); return
      end if
      do ip = 0, np
        do it = 0, nt
          do ir = 0, nr
            read(u, '(6e23.15)', iostat=ios) p%raw(M_K:M_K+5, ir, it, ip)
            if (.not. read_ok(ios, path)) return
          end do
        end do
      end do
      close(u)
      ! The published BBH files include undefined curvature values inside
      ! each companion excision. Those nodes are never interpolation data.
      sepa = geometry(3)
      do ip = 0, np
        do it = 0, nt
          do ir = 0, nr
            dist2 = p%rg(ir)**2+sepa**2-2*sepa*p%rg(ir)*sin(p%thg(it))*cos(p%phig(ip))
            if (m <= 2 .and. dist2 < p%ex_radius**2*(1-1d-12)) cycle
            if (.not. all(ieee_is_finite(p%raw(:, ir, it, ip)))) then
              call fail('COCAL_ID BBH: nonfinite exterior field'); return
            end if
          end do
        end do
      end do
    end associate
    ok = .true.
  end subroutine read_fields

  ! Collective over the node: the disk, as a matter domain on the outer
  ! patch's grid.
  subroutine read_disk_domain(dir_path, ok)
    use COCAL_ID_node_shared
    character(*), intent(in) :: dir_path
    logical, intent(out) :: ok
    real(8) :: status
    call COCAL_ID_domain_init(cocal_matter, 1, 2, .false.)
    associate(p => cocal_matter%patch(1), p3 => cocal_metric%patch(3))
      call COCAL_ID_patch_set_grids(p, p3%rg, p3%thg, p3%phig)
      p%support_index = BBHT_RHO
    end associate
    ok = .true.
    if (COCAL_ID_node_leader()) then
      call COCAL_ID_domain_alloc_raw(cocal_matter, 1)
      call read_disk(dir_path, ok)
    end if
    ! The other ranks of the node take the leader's status and constants.
    status = merge(1d0, 0d0, ok)
    call COCAL_ID_node_bcast(status)
    call COCAL_ID_node_bcast(bbht_k)
    call COCAL_ID_node_bcast(bbht_gamma)
    call COCAL_ID_node_bcast(bbht_rin)
    ok = status == 1
    if (.not. ok) then
      if (.not. COCAL_ID_node_leader()) call fail('COCAL_ID BBHT: the node leader could not read the disk')
      return
    end if
    call COCAL_ID_domain_prepare(cocal_matter, ok, 0.0d0)
  end subroutine read_disk_domain

  ! Node leader only: disk fluid, EOS and edges.
  subroutine read_disk(dir_path, ok)
    character(*), intent(in) :: dir_path
    logical, intent(out) :: ok
    character(1024) :: path
    character(17) :: version
    character(8) :: disk_start
    character(256) :: disk_restart_file
    integer :: u, ios, nr, nt, np, phase, i
    real(8), pointer :: rho(:, :, :), omega(:, :, :)
    real(8) :: logk_rot, ber, disk_rin, disk_rout, disk_rho_target, index_DR_disk
    real(8) :: disk_logk, disk_relax, disk_tolerance, rho_ref, pre_ref, rho_break, gamma_piece
    real(8), allocatable :: r(:), theta(:), phi(:)
    logical :: disk_seed_vacuum
    namelist /bbht_disk/ disk_rin, disk_rout, disk_rho_target, index_DR_disk, &
      disk_logk, disk_relax, disk_tolerance, disk_start, disk_seed_vacuum, disk_restart_file
    ok = .false.
    path = trim(dir_path)//'/disk_state.las'
    open(newunit=u, file=trim(path), status='old', action='read', form='unformatted', iostat=ios)
    if (.not. read_ok(ios, path)) return
    read(u, iostat=ios) version, nr, nt, np
    if (.not. read_ok(ios, path)) return
    if (version /= 'BBHT_BHT_FLUID_V1') then
      call fail('COCAL_ID BBHT: unsupported disk format'); return
    end if
    associate(p => cocal_metric%patch(3))
      if (nr /= p%nr .or. nt /= p%nt .or. np /= p%np) then
        call fail('COCAL_ID BBHT: disk/patch-3 shape mismatch'); return
      end if
      allocate(r(0:nr), theta(0:nt), phi(0:np))
      read(u, iostat=ios) r, theta, phi
      if (.not. read_ok(ios, path)) return
      if (.not. all(ieee_is_finite(r)) .or. .not. all(ieee_is_finite(theta)) .or. &
          .not. all(ieee_is_finite(phi))) then
        call fail('COCAL_ID BBHT: nonfinite disk grid'); return
      end if
      if (any(abs(r-p%rg) > 64*epsilon(1d0)*max(1d0, abs(p%rg))) .or. &
          any(abs(theta-p%thg) > 64*epsilon(1d0)) .or. &
          any(abs(phi-p%phig) > 64*epsilon(1d0))) then
        call fail('COCAL_ID BBHT: disk grid differs from patch 3'); return
      end if
      deallocate(r, theta, phi)
    end associate
    ! K is fitted by the disk solver: peos_parameter.dat's K is not authoritative.
    read(u, iostat=ios) bbht_k, bbht_gamma, logk_rot, ber
    if (.not. read_ok(ios, path)) return
    if (.not. all(ieee_is_finite([bbht_k, bbht_gamma, logk_rot, ber]))) then
      call fail('COCAL_ID BBHT: nonfinite EOS/rotation constants'); return
    end if
    if (bbht_k <= 0 .or. bbht_gamma <= 1) then
      call fail('COCAL_ID BBHT: invalid polytropic EOS'); return
    end if
    associate(raw => cocal_matter%patch(1)%raw)
      rho(0:, 0:, 0:) => raw(BBHT_RHO, :, :, :); omega(0:, 0:, 0:) => raw(BBHT_OMEGA, :, :, :)
    end associate
    read(u, iostat=ios) rho, omega
    if (.not. read_ok(ios, path)) return
    close(u)
    if (.not. all(ieee_is_finite(rho)) .or. .not. all(ieee_is_finite(omega))) then
      call fail('COCAL_ID BBHT: nonfinite fluid data'); return
    end if
    if (any(rho < 0)) then
      call fail('COCAL_ID BBHT: negative density'); return
    end if
    ! V1 saves one Gamma. Explicitly reject a genuinely piecewise-Gamma EOS.
    path = trim(dir_path)//'/peos_parameter.dat'
    open(newunit=u, file=trim(path), status='old', action='read', iostat=ios)
    if (.not. read_ok(ios, path)) return
    read(u, '(8x,i5,es13.5)', iostat=ios) phase, rho_ref
    if (.not. read_ok(ios, path)) return
    if (phase < 0 .or. phase > 100) then
      call fail('COCAL_ID BBHT: invalid EOS phase count'); return
    end if
    read(u, '(2es13.5)', iostat=ios) rho_ref, pre_ref
    if (.not. read_ok(ios, path)) return
    do i = phase, 0, -1
      read(u, '(2es13.5)', iostat=ios) rho_break, gamma_piece
      if (.not. read_ok(ios, path)) return
      if (.not. ieee_is_finite(gamma_piece)) then
        call fail('COCAL_ID BBHT: nonfinite EOS Gamma'); return
      end if
      ! peos_parameter.dat prints Gamma as es13.5 (six significant digits), so
      ! allow its rounding; a full-precision disk Gamma must still match.
      if (abs(gamma_piece-bbht_gamma) > 1d-5*abs(bbht_gamma)) then
        call fail('COCAL_ID BBHT: V1 reader requires a single Gamma matching disk_state.las'); return
      end if
    end do
    close(u)
    ! The exact zero-density inner surface is not generally a radial grid node.
    ! Require saved edges; accept the other solver controls for namelist compatibility.
    disk_rin = -1; disk_rout = -1; disk_rho_target = 0; index_DR_disk = 0
    disk_logk = 0; disk_relax = 0; disk_tolerance = 0
    disk_start = ''; disk_seed_vacuum = .false.; disk_restart_file = ''
    path = trim(dir_path)//'/diskpar.dat'
    open(newunit=u, file=trim(path), status='old', action='read', iostat=ios)
    if (.not. read_ok(ios, path)) return
    read(u, nml=bbht_disk, iostat=ios)
    if (.not. read_ok(ios, path)) return
    close(u)
    bbht_rin = disk_rin
    if (.not. all(ieee_is_finite([disk_rin, disk_rout]))) then
      call fail('COCAL_ID BBHT: nonfinite disk edges'); return
    end if
    if (disk_rin < cocal_metric%patch(3)%rg(0) .or. disk_rout <= disk_rin .or. &
        disk_rout > cocal_metric%patch(3)%rg(nr)) then
      call fail('COCAL_ID BBHT: unsupported disk edges'); return
    end if
    ok = .true.
  end subroutine read_disk

  ! Disk fluid at (x,y,z) with metric fields q there; saved Omega is dphi/dt
  ! about the global z axis.
  subroutine COCAL_ID_bbht_matter(x, y, z, q, s)
    real(8), intent(in) :: x, y, z, q(:)
    real(8), intent(inout) :: s(NSTATE)
    real(8) :: f(2), emd
    logical :: ok
    if (x*x+y*y+z*z <= bbht_rin**2) return
    call COCAL_ID_domain_eval(cocal_matter, x, y, z, f, ok)
    if (.not. ok .or. .not. f(BBHT_RHO) > 0) return
    emd = bbht_k*f(BBHT_RHO)**(bbht_gamma-1)
    s(S_RHO) = f(BBHT_RHO)
    s(S_PRESS) = f(BBHT_RHO)*emd
    s(S_EPS) = emd/(bbht_gamma-1)
    s(S_VEL:S_VEL+2) = (f(BBHT_OMEGA)*[-y, x, 0d0]+q(M_BETA:M_BETA+2))/q(M_ALP)
  end subroutine COCAL_ID_bbht_matter

end module COCAL_ID_bbh
