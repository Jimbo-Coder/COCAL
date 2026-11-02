! COCAL fields on COCAL's spherical coordinate patches.
!
! A domain holds one patch (RNS, a star's fluid) or the three patches of a
! binary (BNS): one centred on each star, rotated by pi about z for the
! second, and the outer patch about the common centre. Every
! field is the C2 spline of cocal_id_spline.f90, and the three metric patches
! are blended across their overlaps with C2 weights (below). The node leader
! reads the exported fields into raw and fits the node-shared coefficients.
!
! Per patch, a reader can mark fields that
!   - change sign between the stored and the global Cartesian basis (flip),
!   - must be positive at a valid point (positive),
! A surface-fitted patch has radial coordinate r/rs(theta,phi), as on
! COCAL's fluid grids, with rs exported beside the fields.

module COCAL_ID_domain
  use COCAL_ID_spline
  use COCAL_ID_node_shared
  use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
  implicit none
  private


  public :: cocal_patch, cocal_domain
  public :: COCAL_ID_domain_init, COCAL_ID_patch_set_grids, COCAL_ID_domain_alloc_raw, COCAL_ID_domain_prepare
  public :: COCAL_ID_domain_eval, COCAL_ID_domain_patch_eval, COCAL_ID_domain_free

  type :: cocal_patch
    integer :: nr = 0, nt = 0, np = 0
    real(8) :: dis = 0.0d0       ! distance of the patch centre from the domain centre
    real(8) :: ex_radius = 0.0d0 ! companion excision radius (star patches)
    integer, allocatable :: flip_indices(:), positive_indices(:)
    logical :: surface_fitted = .false.
    real(8), allocatable :: rg(:), thg(:), phig(:) ! rg(0:nr), thg(0:nt), phig(0:np)

    ! Exported fields raw(1:nvars, 0:nr, 0:nt, 0:np) and surface
    ! raw_surface(0:nt, 0:np), node leader only until the fit. Readers point
    ! a 3D view at each field, e.g.
    !   psi(0:,0:,0:) => p%raw(PSI,:,:,:)
    ! and COCAL's reading routines fill it in place.
    real(8), allocatable :: raw(:,:,:,:), raw_surface(:,:)

    ! Node-shared spline coefficients (1:nvars, 0:nr, 0:nt, 0:np-1), surface
    ! coefficients (1, 0, 0:nt, 0:np-1) and, for matter, the marks of rough
    ! angular data (cocal_id_spline.f90) (0:nr, 0:nt, 0:np-1).
    real(8), pointer :: sf(:,:,:,:) => null(), sm(:,:,:,:) => null(), rough(:,:,:) => null()
    real(8), pointer :: surface(:,:,:,:) => null()
  end type cocal_patch

  type :: cocal_domain
    integer :: npatches = 0
    integer :: nvars = 0
    ! A metric keeps the C2 spline everywhere (the constraints take its second
    ! derivatives), with COCAL's centre value replaced by the smooth
    ! continuation of its neighbours. Matter fields need only their values and
    ! fall back to the local stencil where their data are rough; a metric's
    ! rough spots are the seams of COCAL's binary patches, where that does not
    ! help (cocal_id_spline.f90).
    logical :: metric = .false.
    ! Radius beyond which only the outer patch is used, clamped into the
    ! outer blend shell; the hard switch radius when the patches do not overlap.
    real(8) :: rr3 = 0.0d0
    ! Separation of the star patch centres and their midpoint on the x axis.
    real(8) :: sepa = 0.0d0, xmid = 0.0d0
    real(8) :: slab_halfwidth = 0.0d0
    real(8) :: shell_lo = 0.0d0, shell_hi = 0.0d0
    type(cocal_patch), allocatable :: patch(:)
  end type cocal_domain

contains

  subroutine COCAL_ID_domain_init(dom, npatches, nvars, metric)
    type(cocal_domain), intent(inout) :: dom
    integer, intent(in) :: npatches, nvars
    logical, intent(in) :: metric

    call COCAL_ID_domain_free(dom)
    dom%npatches = npatches
    dom%nvars = nvars
    dom%metric = metric
    allocate(dom%patch(npatches))
  end subroutine COCAL_ID_domain_init

  subroutine COCAL_ID_patch_set_grids(p, rg, thg, phig)
    type(cocal_patch), intent(inout) :: p
    real(8), intent(in) :: rg(0:), thg(0:), phig(0:)

    p%nr = ubound(rg, 1)
    p%nt = ubound(thg, 1)
    p%np = ubound(phig, 1)
    p%rg = rg
    p%thg = thg
    p%phig = phig
  end subroutine COCAL_ID_patch_set_grids

  ! Node leader only, once patch m's grids are set; the leader then fills
  ! raw(:, 0:nr, 0:nt, 0:np), and raw_surface if surface fitted, both zeroed
  ! here.
  subroutine COCAL_ID_domain_alloc_raw(dom, m)
    type(cocal_domain), intent(inout) :: dom
    integer, intent(in) :: m

    associate(p => dom%patch(m))
      allocate(p%raw(dom%nvars, 0:p%nr, 0:p%nt, 0:p%np))
      p%raw = 0.0d0
      if (p%surface_fitted) then
        allocate(p%raw_surface(0:p%nt, 0:p%np))
        p%raw_surface = 0.0d0
      end if
    end associate
  end subroutine COCAL_ID_domain_alloc_raw

  ! Collective over the node: blend zones on every rank; companion fill and
  ! spline fit by the node leader into node-shared coefficients.
  ! ok=.false. when a patch's grids are not increasing radii and the uniform
  ! angular grids the spline assumes.
  ! star_extent: largest stellar radius about a star patch centre, kept clear
  ! of the midplane blend (0 for holes).
  subroutine COCAL_ID_domain_prepare(dom, ok, star_extent)
    type(cocal_domain), intent(inout) :: dom
    logical, intent(out) :: ok
    real(8), intent(in) :: star_extent
    logical :: leader
    integer :: m, iv
    real(8) :: scale(dom%nvars)

    leader = COCAL_ID_node_leader()
    if (dom%npatches == 3) call set_blend(dom, star_extent)
    ok = .true.
    do m = 1, dom%npatches
      associate(p => dom%patch(m))
        ok = ok .and. all(p%rg(1:) > p%rg(:p%nr - 1)) .and. COCAL_ID_spline_grids_ok(p%thg, p%phig)
        call COCAL_ID_node_alloc_4d(p%sf, dom%nvars, p%nr, p%nt, p%np)
        call COCAL_ID_node_alloc_4d(p%sm, dom%nvars, p%nr, p%nt, p%np)
        if (.not. dom%metric) call COCAL_ID_node_alloc_3d(p%rough, p%nr, p%nt, p%np - 1)
        if (p%surface_fitted) call COCAL_ID_node_alloc_4d(p%surface, 1, 0, p%nt, p%np)
      end associate
    end do
    if (.not. ok) return

    if (leader) then
      do m = 1, dom%npatches
        associate(p => dom%patch(m))
          p%sf = p%raw(:, 0:p%nr, :, 0:p%np - 1)
        end associate
      end do
      if (dom%npatches == 3) call fill_companion(dom)
      do m = 1, dom%npatches
        associate(p => dom%patch(m))
          if (associated(p%rough)) scale = [(maxval(abs(p%sf(iv, :, :, :))), iv = 1, dom%nvars)]
          call COCAL_ID_spline_fit(p%rg, p%nt, p%np, p%sf, p%sm, dom%metric)
          if (associated(p%rough)) call COCAL_ID_spline_rough(p%nt, p%np, p%sf, scale, p%rough)
          deallocate(p%raw)
          if (associated(p%surface)) then
            p%surface(1, 0, :, :) = p%raw_surface(:, 0:p%np - 1)
            call COCAL_ID_spline_angular(p%surface, 0, p%nt, p%np)
            deallocate(p%raw_surface)
          end if
        end associate
      end do
    end if
    call COCAL_ID_node_barrier()
  end subroutine COCAL_ID_domain_prepare

  ! Fields at domain coordinates (x,y,z), in global Cartesian components.
  subroutine COCAL_ID_domain_eval(dom, x, y, z, vals, ok)
    type(cocal_domain), intent(in) :: dom
    real(8), intent(in) :: x, y, z
    real(8), intent(out) :: vals(dom%nvars)
    logical, intent(out) :: ok
    real(8) :: r, wmid, wout, q(dom%nvars)
    logical :: okp

    if (dom%npatches == 1) then
      call patch_eval(dom, 1, x, y, z, vals, ok, .false.)
      return
    end if

    r = sqrt(x*x + y*y + z*z)
    wmid = blend_weight(x - dom%xmid, -dom%slab_halfwidth, dom%slab_halfwidth)
    wout = blend_weight(r, dom%shell_lo, dom%shell_hi)
    vals = 0.0d0
    ok = .true.
    if (wout < 1.0d0 .and. wmid < 1.0d0) then
      call patch_eval(dom, 1, x + dom%patch(1)%dis, y, z, q, okp, .false.)
      vals = vals + (1.0d0 - wout)*(1.0d0 - wmid)*q
      ok = ok .and. okp
    end if
    if (wout < 1.0d0 .and. wmid > 0.0d0) then
      call patch_eval(dom, 2, dom%patch(2)%dis - x, -y, z, q, okp, .false.)
      vals = vals + (1.0d0 - wout)*wmid*q
      ok = ok .and. okp
    end if
    if (wout > 0.0d0) then
      call patch_eval(dom, 3, x, y, z, q, okp, .false.)
      vals = vals + wout*q
      ok = ok .and. okp
    end if
  end subroutine COCAL_ID_domain_eval

  ! Patch m alone, at patch coordinates (xc,yc,zc), in global components.
  subroutine COCAL_ID_domain_patch_eval(dom, m, xc, yc, zc, vals, ok)
    type(cocal_domain), intent(in) :: dom
    integer, intent(in) :: m
    real(8), intent(in) :: xc, yc, zc
    real(8), intent(out) :: vals(dom%nvars)
    logical, intent(out) :: ok

    call patch_eval(dom, m, xc, yc, zc, vals, ok, .false.)
  end subroutine COCAL_ID_domain_patch_eval

  ! Node-shared blocks are released by COCAL_ID_node_free_all.
  subroutine COCAL_ID_domain_free(dom)
    type(cocal_domain), intent(inout) :: dom

    if (allocated(dom%patch)) deallocate(dom%patch)
    dom = cocal_domain()
  end subroutine COCAL_ID_domain_free

  ! Patch m at patch coordinates (xc,yc,zc), in global Cartesian components:
  ! the spline, or (lagrange, before the fit) the Lagrange stencil on raw.
  subroutine patch_eval(dom, m, xc, yc, zc, vals, ok, lagrange)
    type(cocal_domain), intent(in) :: dom
    integer, intent(in) :: m
    real(8), intent(in) :: xc, yc, zc
    real(8), intent(out) :: vals(dom%nvars)
    logical, intent(out) :: ok
    logical, intent(in) :: lagrange
    real(8) :: s

    associate(p => dom%patch(m))
      if (lagrange) then
        call lagrange_eval(p, dom%nvars, dom%sepa, xc, yc, zc, vals, ok)
      else
        s = 1.0d0
        if (associated(p%surface)) s = 1.0d0/COCAL_ID_spline_eval_angular(p%nt, p%np, p%surface, xc, yc, zc)
        ! A metric's rough is null, which leaves the optional argument absent.
        call COCAL_ID_spline_eval(p%rg, p%nt, p%np, p%sf, p%sm, s*xc, s*yc, s*zc, vals, ok, rough=p%rough)
      end if
      if (.not. ok) return
      if (allocated(p%flip_indices)) vals(p%flip_indices) = -vals(p%flip_indices)
      ok = all(ieee_is_finite(vals))
      if (ok .and. allocated(p%positive_indices)) ok = all(vals(p%positive_indices) > 0.0d0)
    end associate
  end subroutine patch_eval

  ! Star-patch nodes inside the companion excision take the companion patch's
  ! value at the same point, as COCAL's interpolation_fillup does, so the
  ! fit sees no stale data.
  subroutine fill_companion(dom)
    type(cocal_domain), intent(inout) :: dom
    real(8) :: xc, yc, zc, v(dom%nvars)
    integer :: m, ir, it, ip
    logical :: ok

    do m = 1, 2
      associate(p => dom%patch(m))
        do ip = 0, p%np - 1
          do it = 0, p%nt
            do ir = 0, p%nr
              xc = p%rg(ir)*sin(p%thg(it))*cos(p%phig(ip))
              yc = p%rg(ir)*sin(p%thg(it))*sin(p%phig(ip))
              zc = p%rg(ir)*cos(p%thg(it))
              if ((xc - dom%sepa)**2 + yc*yc + zc*zc >= p%ex_radius**2) cycle
              ! The companion's frame is this one rotated by pi about z and
              ! centred sepa along x.
              call patch_eval(dom, 3 - m, dom%sepa - xc, -yc, zc, v, ok, .true.)
              if (.not. ok) cycle
              if (allocated(p%flip_indices)) v(p%flip_indices) = -v(p%flip_indices)
              p%sf(:, ir, it, ip) = v
            end do
          end do
        end do
      end associate
    end do
  end subroutine fill_companion

  ! COCAL's 4th-order Lagrange stencil on the exported fields raw of patch p,
  ! on uniform angular grids, continued over the poles and in phi. Near a
  ! companion excision the radial stencil is shifted so that it reads no
  ! excised node while still bracketing the point.
  subroutine lagrange_eval(p, nvars, sepa, xc, yc, zc, vals, ok)
    type(cocal_patch), intent(in) :: p
    integer, intent(in) :: nvars
    real(8), intent(in) :: sepa, xc, yc, zc
    real(8), intent(out) :: vals(nvars)
    logical, intent(out) :: ok
    real(8), parameter :: pi = 3.14159265358979323846264338327950d0
    real(8) :: rc, thc, phic, tol, r4(4), th4(4), phi4(4), wr4(4), wth4(4), wphi4(4), w, f(nvars)
    integer :: lower, ir0, it0, ip0, ii, jj, kk, it, ip, trial, starts(4), its(4), ips(4)

    vals = 0.0d0
    ok = .false.
    rc = sqrt(xc*xc + yc*yc + zc*zc)
    if (p%ex_radius > 0.0d0 .and. (xc - sepa)**2 + yc*yc + zc*zc < p%ex_radius**2) return
    thc = atan2(sqrt(xc*xc + yc*yc), zc)
    phic = modulo(atan2(yc, xc), 2.0d0*pi)
    lower = lower_index(p%rg(0:p%nr), rc)
    it0 = lower_index(p%thg(0:p%nt), thc) - 1
    ip0 = lower_index(p%phig(0:p%np), phic) - 1
    do jj = 1, 4
      it = it0 + jj - 1
      if (it < 0) then
        th4(jj) = -p%thg(-it)
      else if (it > p%nt) then
        th4(jj) = 2.0d0*pi - p%thg(2*p%nt - it)
      else
        th4(jj) = p%thg(it)
      end if
      its(jj) = min(abs(it), 2*p%nt - it)
    end do
    do kk = 1, 4
      ip = ip0 + kk - 1
      ips(kk) = modulo(ip, p%np)
      phi4(kk) = p%phig(ips(kk)) + floor(real(ip, 8)/p%np)*2.0d0*pi
    end do
    ! Across a pole the stencil continues at phi+pi.
    do jj = 1, 4
      if (it0 + jj - 1 < 0 .or. it0 + jj - 1 > p%nt) its(jj) = -its(jj) - 1
    end do

    tol = 64.0d0*epsilon(rc)*max(1.0d0, rc)
    starts = [max(0, min(lower - 1, p%nr - 3)), max(0, min(lower - 2, p%nr - 3)), &
              max(0, min(lower - 3, p%nr - 3)), max(0, min(lower, p%nr - 3))]
    do trial = 1, 4
      ir0 = starts(trial)
      ok = rc >= p%rg(ir0) - tol .and. rc <= p%rg(ir0 + 3) + tol
      do kk = 1, 4
        do jj = 1, 4
          call stencil_node(jj, kk, it, ip)
          do ii = 0, 3
            if (p%ex_radius > 0.0d0 .and. p%rg(ir0 + ii)**2 + sepa**2 - 2.0d0*sepa*p%rg(ir0 + ii)* &
                sin(p%thg(it))*cos(p%phig(ip)) < p%ex_radius**2*(1.0d0 - 1.0d-12)) ok = .false.
          end do
        end do
      end do
      if (ok) exit
    end do
    if (.not. ok) return

    r4 = p%rg(ir0:ir0 + 3)
    call COCAL_ID_lagint_4th_weights(r4, rc, wr4)
    call COCAL_ID_lagint_4th_weights(th4, thc, wth4)
    call COCAL_ID_lagint_4th_weights(phi4, phic, wphi4)
    do kk = 1, 4
      do jj = 1, 4
        call stencil_node(jj, kk, it, ip)
        do ii = 1, 4
          w = wr4(ii)*wth4(jj)*wphi4(kk)
          f = p%raw(:, ir0 + ii - 1, it, ip)
          vals = vals + w*f
        end do
      end do
    end do

  contains
    subroutine stencil_node(jj, kk, it, ip)
      integer, intent(in) :: jj, kk
      integer, intent(out) :: it, ip
      it = its(jj)
      ip = ips(kk)
      if (it < 0) then
        it = -it - 1
        ip = modulo(ip + p%np/2, p%np)
      end if
    end subroutine stencil_node
  end subroutine lagrange_eval

  pure integer function lower_index(grid, x) result(lo)
    real(8), intent(in) :: grid(0:), x
    integer :: hi, mid
    lo = 0
    hi = ubound(grid, 1)
    do while (hi - lo > 1)
      mid = (hi + lo)/2
      if (x < grid(mid)) then
        hi = mid
      else
        lo = mid
      end if
    end do
  end function lower_index

  ! Blending patches across their overlaps.
  !
  ! COCAL solves each patch on its own spherical grid and couples the patches
  ! only through surface data (sourceterm_exsurf/insurf/outsurf in PCOCAL), so
  ! where two patches overlap both approximate the same field and differ only
  ! by COCAL's discretization error D. A discrete switch from one patch to the
  ! other leaves a slope jump that Cartesian differencing turns into a
  ! constraint error ~|grad D|/h, which grows without bound under refinement.
  ! Blending q = (1-w) q_a + w q_b across a zone of width 2L instead adds only
  ! 2 grad(w).grad(D) + D lap(w) ~ |grad D|/L + |D|/L**2, independent of h.
  ! Both terms fall with L, so each zone spans the whole overlap in which both
  ! patches are valid, less two cells at each edge to keep the 4-point Lagrange
  ! stencils clear of the excised region. The quintic weight is C2, as the
  ! second derivatives in the constraints require. There is no free parameter.
  !
  ! The zones follow from the patch grids alone: the midplane slab both star
  ! patches cover, kept clear of the stars, and the outer shell between the
  ! star patches' edges and the outer patch.
  subroutine set_blend(dom, star_extent)
    type(cocal_domain), intent(inout) :: dom
    real(8), intent(in) :: star_extent
    real(8) :: shell_edge
    integer :: m

    dom%sepa = dom%patch(1)%dis + dom%patch(2)%dis
    dom%xmid = (dom%patch(2)%dis - dom%patch(1)%dis)/2.0d0
    dom%slab_halfwidth = huge(1.0d0)
    dom%shell_hi = huge(1.0d0)
    do m = 1, 2
      associate(p => dom%patch(m))
        ! COCAL puts the companion at distance sepa with an excision ball of
        ! radius ex_radius, so the patch covers the midplane slab out to
        ! sepa/2 - ex_radius on the companion side.
        dom%slab_halfwidth = min(dom%slab_halfwidth, dom%sepa/2.0d0 - p%ex_radius &
          - 2.0d0*blend_cell(p%rg, dom%sepa - p%ex_radius))
        shell_edge = p%rg(p%nr) - 2.0d0*blend_cell(p%rg, p%rg(p%nr))
        dom%shell_hi = min(dom%shell_hi, shell_edge - p%dis)
      end associate
    end do
    associate(p1 => dom%patch(1), p3 => dom%patch(3))
      if (star_extent > 0.0d0) dom%slab_halfwidth = min(dom%slab_halfwidth, &
        p1%dis - star_extent - 2.0d0*blend_cell(p1%rg, star_extent))
      dom%shell_lo = max(p3%rg(0), max(p1%dis, dom%patch(2)%dis) + star_extent) + &
        2.0d0*blend_cell(p3%rg, p3%rg(0))
    end associate
    dom%slab_halfwidth = max(0.0d0, dom%slab_halfwidth)

    if (dom%shell_hi > dom%shell_lo) then
      dom%rr3 = min(max(dom%rr3, dom%shell_lo), dom%shell_hi)
    else
      ! No usable overlap: hard switch at rr3.
      dom%shell_lo = dom%rr3
      dom%shell_hi = dom%rr3
    end if
  end subroutine set_blend

  ! Weight of the second patch: 0 for u <= lo, 1 for u >= hi; with no
  ! overlap (lo = hi) a step, 0 at u = hi.
  pure real(8) function blend_weight(u, lo, hi)
    real(8), intent(in) :: u, lo, hi
    real(8) :: t
    if (hi <= lo) then
      blend_weight = merge(1.0d0, 0.0d0, u > hi)
      return
    end if
    t = min(1.0d0, max(0.0d0, (u - lo)/(hi - lo)))
    blend_weight = t**3*(10.0d0 - 15.0d0*t + 6.0d0*t*t)
  end function blend_weight

  ! Larger of the two radial cells of rg(0:) adjacent to radius r.
  pure real(8) function blend_cell(rg, r)
    real(8), intent(in) :: rg(0:), r
    integer :: i, n
    n = ubound(rg, 1)
    i = 1
    do while (i < n .and. rg(i) < r)
      i = i + 1
    end do
    blend_cell = rg(i) - rg(i - 1)
    if (i < n) blend_cell = max(blend_cell, rg(i + 1) - rg(i))
  end function blend_cell

end module COCAL_ID_domain
