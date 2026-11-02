! Black-hole interiors, shared by BBH, BBHT and BHT.
!
! COCAL excises each hole; inside a fill radius rh the metric is continued
! from the exterior data of the hole's own patch, in fill variables that keep
! it valid: log psi**4, log alpha, the covariant shift beta_i, the physical
! tracefree curvature A_ij = K_ij - gamma_ij K/3, h_ij for waveless data, and
! the trace K. Each is continued as follows:
!   - centre value: four samples r_k = rh + (k-1)*dr along 14 directions (6
!     axes, 8 cube diagonals), each extrapolated to r=0 as a cubic in r**2
!     (even in r), averaged over the directions; K is taken to zero there;
!   - boundary data: value and first two radial derivatives at rh along the
!     point's own ray, from one-sided differences of the exterior on six
!     samples rh + (k-1)*0.02*rh: fourth-order accurate, with roundoff
!     amplified only ~2500 times;
!   - a quintic Hermite in s=r/rh matching the boundary data (C2 at rh) and
!     flat through second order at the centre.
! The centre spacing dr follows the exterior grid spacing near rh, bounded to
! [0.2, 0.25]*rh. Holes hold no matter.

module COCAL_ID_holes
  use COCAL_ID_state
  use COCAL_ID_domain
  use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
  implicit none
  private

  ! Centre samples per direction, and boundary samples along a ray.
  integer, parameter :: nsample = 4, nray = 6
  ! Boundary sample spacing, relative to rh.
  real(8), parameter :: ray_step = 0.02d0

  public :: cocal_hole, cocal_holes, cocal_fill_holes
  public :: COCAL_ID_holes_add, COCAL_ID_holes_extrapolate_centre
  public :: COCAL_ID_holes_eval, COCAL_ID_holes_free

  type :: cocal_hole
    integer :: patch = 1               ! metric patch whose exterior fills the hole
    real(8) :: centre(3) = 0.0d0      ! in domain coordinates
    real(8) :: patch_centre(3) = 0.0d0 ! in that patch's coordinates
    logical :: rotated = .false.       ! patch frame rotated by pi about z
    real(8) :: rh = 0.0d0, dr = 0.0d0
    real(8), allocatable :: q0(:)     ! fill variables at the centre
  end type cocal_hole

  type(cocal_hole), allocatable, save :: cocal_holes(:)
  ! Fill the interiors; otherwise points inside a hole are rejected.
  logical, save :: cocal_fill_holes = .false.

contains

  ! Hole about patch_centre of metric patch m, at centre in domain
  ! coordinates, with fill radius rh.
  subroutine COCAL_ID_holes_add(m, centre, patch_centre, rotated, rh)
    integer, intent(in) :: m
    real(8), intent(in) :: centre(3), patch_centre(3), rh
    logical, intent(in) :: rotated
    type(cocal_hole), allocatable :: holes(:)
    type(cocal_hole) :: h
    integer :: ir

    h%patch = m
    h%centre = centre
    h%patch_centre = patch_centre
    h%rotated = rotated
    h%rh = rh
    associate(rg => cocal_metric%patch(m)%rg)
      ir = 0
      do while (ir < ubound(rg, 1) - 1 .and. rg(ir + 1) <= rh)
        ir = ir + 1
      end do
      h%dr = fill_spacing(rh, rg(ir + 1) - rg(ir))
    end associate
    allocate(h%q0(nfill()))
    h%q0 = 0.0d0
    if (allocated(cocal_holes)) then
      holes = [cocal_holes, h]
    else
      holes = [h]
    end if
    call move_alloc(holes, cocal_holes)
  end subroutine COCAL_ID_holes_add

  ! Centre values of hole ih extrapolated from the exterior.
  subroutine COCAL_ID_holes_extrapolate_centre(ih, ok)
    integer, intent(in) :: ih
    logical, intent(out) :: ok

    associate(h => cocal_holes(ih))
      call fill_centre(ih, nfill(), h%rh, h%dr, h%q0, ok)
      h%q0(nfill()) = 0.0d0
    end associate
  end subroutine COCAL_ID_holes_extrapolate_centre

  ! Metric fields q at domain coordinates (x,y,z) inside a hole, if any.
  subroutine COCAL_ID_holes_eval(x, y, z, q, inside, ok)
    real(8), intent(in) :: x, y, z
    real(8), intent(out) :: q(:)
    logical, intent(out) :: inside, ok
    real(8) :: rel(3), f(nfill())
    integer :: ih

    inside = .false.
    ok = .false.
    if (.not. allocated(cocal_holes)) return
    do ih = 1, size(cocal_holes)
      associate(h => cocal_holes(ih))
        rel = [x, y, z] - h%centre
        if (h%rotated) rel(1:2) = -rel(1:2)
        if (sum(rel**2) >= h%rh**2) cycle
        inside = .true.
        if (.not. cocal_fill_holes) return
        call fill_ray(ih, nfill(), h%rh, h%q0, rel(1), rel(2), rel(3), f, ok)
        if (ok) call from_fill(f, q, ok)
        return
      end associate
    end do
  end subroutine COCAL_ID_holes_eval

  subroutine COCAL_ID_holes_free()
    if (allocated(cocal_holes)) deallocate(cocal_holes)
    cocal_fill_holes = .false.
  end subroutine COCAL_ID_holes_free

  ! The metric fields, then the trace K.
  pure integer function nfill()
    nfill = cocal_metric%nvars + 1
  end function nfill

  ! Fill variables at hole-centred (x,y,z): the hole's own patch, in the
  ! global Cartesian basis.
  subroutine sample(ih, x, y, z, f, ok)
    integer, intent(in) :: ih
    real(8), intent(in) :: x, y, z
    real(8), intent(out) :: f(:)
    logical, intent(out) :: ok
    real(8) :: q(cocal_metric%nvars), c(3)

    f = 0.0d0
    c = cocal_holes(ih)%patch_centre
    call COCAL_ID_domain_patch_eval(cocal_metric, cocal_holes(ih)%patch, x + c(1), y + c(2), z + c(3), q, ok)
    if (ok) call to_fill(q, f, ok)
  end subroutine sample

  subroutine to_fill(q, f, ok)
    real(8), intent(in) :: q(:)
    real(8), intent(out) :: f(:)
    logical, intent(out) :: ok
    real(8) :: gt(6), gi(6), det, psi4, trk

    call conformal_metric(q, gt)
    call sym_inverse(gt, gi, det)
    ok = det > 0.0d0
    if (.not. ok) return
    psi4 = q(M_PSI)**4
    trk = sym_trace(gi, q(M_K:M_K + 5))
    f = q(:size(f) - 1)
    f(M_PSI) = log(psi4)
    f(M_ALP) = log(q(M_ALP))
    f(M_BETA:M_BETA + 2) = psi4*sym_mul(gt, q(M_BETA:M_BETA + 2))
    f(M_K:M_K + 5) = psi4*(q(M_K:M_K + 5) - gt*trk/3.0d0)
    f(size(f)) = trk
  end subroutine to_fill

  ! Inverse of to_fill, with A_ij made tracefree again through A_yy.
  subroutine from_fill(f, q, ok)
    real(8), intent(in) :: f(:)
    real(8), intent(out) :: q(:)
    logical, intent(out) :: ok
    real(8) :: gt(6), gi(6), det, psi4, a(6)

    q = f(:size(q))
    psi4 = exp(f(M_PSI))
    q(M_PSI) = psi4**0.25d0
    q(M_ALP) = exp(f(M_ALP))
    call conformal_metric(q, gt)
    call sym_inverse(gt, gi, det)
    ok = det > 0.0d0
    if (.not. ok) return
    q(M_BETA:M_BETA + 2) = sym_mul(gi, f(M_BETA:M_BETA + 2))/psi4
    a = f(M_K:M_K + 5)/psi4
    a(4) = -(gi(1)*a(1) + gi(6)*a(6) + 2.0d0*(gi(2)*a(2) + gi(3)*a(3) + gi(5)*a(5)))/gi(4)
    q(M_K:M_K + 5) = a + gt*f(size(f))/3.0d0
    ok = all(ieee_is_finite(q))
  end subroutine from_fill

  subroutine conformal_metric(q, gt)
    real(8), intent(in) :: q(:)
    real(8), intent(out) :: gt(6)

    gt = [1.0d0, 0.0d0, 0.0d0, 1.0d0, 0.0d0, 1.0d0]
    if (cocal_wl) gt = gt + q(M_H:M_H + 5)
  end subroutine conformal_metric

  ! Symmetric 3x3 tensors as (xx, xy, xz, yy, yz, zz).
  pure subroutine sym_inverse(g, gi, det)
    real(8), intent(in) :: g(6)
    real(8), intent(out) :: gi(6), det

    gi(1) = g(4)*g(6) - g(5)*g(5)
    gi(2) = g(3)*g(5) - g(2)*g(6)
    gi(3) = g(2)*g(5) - g(3)*g(4)
    gi(4) = g(1)*g(6) - g(3)*g(3)
    gi(5) = g(2)*g(3) - g(1)*g(5)
    gi(6) = g(1)*g(4) - g(2)*g(2)
    det = g(1)*gi(1) + g(2)*gi(2) + g(3)*gi(3)
    if (det > 0.0d0) gi = gi/det
  end subroutine sym_inverse

  pure function sym_mul(g, v) result(w)
    real(8), intent(in) :: g(6), v(3)
    real(8) :: w(3)
    w = [g(1)*v(1) + g(2)*v(2) + g(3)*v(3), g(2)*v(1) + g(4)*v(2) + g(5)*v(3), &
         g(3)*v(1) + g(5)*v(2) + g(6)*v(3)]
  end function sym_mul

  pure real(8) function sym_trace(gi, k)
    real(8), intent(in) :: gi(6), k(6)
    sym_trace = gi(1)*k(1) + gi(4)*k(4) + gi(6)*k(6) + 2.0d0*(gi(2)*k(2) + gi(3)*k(3) + gi(5)*k(5))
  end function sym_trace

  pure real(8) function fill_spacing(rh, dr_grid)
    real(8), intent(in) :: rh, dr_grid
    fill_spacing = min(0.25d0*rh, max(2.0d0*dr_grid, 0.20d0*rh))
  end function fill_spacing

  subroutine fill_centre(hole, nq, rh, dr, center, ok)
    integer, intent(in) :: hole, nq
    real(8), intent(in) :: rh, dr
    real(8), intent(out) :: center(nq)
    logical, intent(out) :: ok
    real(8) :: dirs(3, 14), q(nq, nsample), r(nsample), w(nsample)
    integer :: d, k
    external :: COCAL_ID_lagint_4th_weights

    dirs(:, 1) = [ 1.0d0,  0.0d0,  0.0d0]
    dirs(:, 2) = [-1.0d0,  0.0d0,  0.0d0]
    dirs(:, 3) = [ 0.0d0,  1.0d0,  0.0d0]
    dirs(:, 4) = [ 0.0d0, -1.0d0,  0.0d0]
    dirs(:, 5) = [ 0.0d0,  0.0d0,  1.0d0]
    dirs(:, 6) = [ 0.0d0,  0.0d0, -1.0d0]
    do d = 0, 7
      dirs(:, d + 7) = [real(2*ibits(d, 0, 1) - 1, 8), real(2*ibits(d, 1, 1) - 1, 8), &
                        real(2*ibits(d, 2, 1) - 1, 8)]/sqrt(3.0d0)
    end do

    r = [(rh + dr*(k - 1), k = 1, nsample)]
    call COCAL_ID_lagint_4th_weights(r*r, 0.0d0, w)
    center = 0.0d0
    do d = 1, 14
      do k = 1, nsample
        call sample(hole, r(k)*dirs(1, d), r(k)*dirs(2, d), r(k)*dirs(3, d), q(:, k), ok)
        if (.not. ok) return
      end do
      center = center + matmul(q, w)/14.0d0
    end do
  end subroutine fill_centre

  ! Fill variables at hole-centred (x,y,z) with r = |(x,y,z)| <= rh.
  subroutine fill_ray(hole, nq, rh, center, x, y, z, q, ok)
    integer, intent(in) :: hole, nq
    real(8), intent(in) :: rh, center(nq), x, y, z
    real(8), intent(out) :: q(nq)
    logical, intent(out) :: ok
    real(8) :: u(3), r, s, qs(nq, nray), rk(nray), c(nray, 0:2)
    real(8) :: h0, h1, hd1, hdd1
    integer :: k

    q = 0.0d0
    r = sqrt(x*x + y*y + z*z)
    u = [1.0d0, 0.0d0, 0.0d0]
    if (r > 0.0d0) u = [x, y, z]/r
    rk = [(rh*(1.0d0 + ray_step*(k - 1)), k = 1, nray)]
    do k = 1, nray
      call sample(hole, rk(k)*u(1), rk(k)*u(2), rk(k)*u(3), qs(:, k), ok)
      if (.not. ok) return
    end do
    call derivative_weights(rk, rh, c)

    s = min(1.0d0, r/rh)
    h1 = s**3*(10.0d0 - 15.0d0*s + 6.0d0*s*s)
    h0 = 1.0d0 - h1
    hd1 = rh*s**3*(-4.0d0 + 7.0d0*s - 3.0d0*s*s)
    hdd1 = 0.5d0*rh*rh*s**3*(1.0d0 - 2.0d0*s + s*s)
    q = h0*center + h1*qs(:, 1) + hd1*matmul(qs, c(:, 1)) + hdd1*matmul(qs, c(:, 2))
  end subroutine fill_ray

  ! Weights c(:,m) of the m-th derivative at x0 of the polynomial through xs
  ! (Fornberg, Math. Comp. 51, 699 (1988)).
  pure subroutine derivative_weights(xs, x0, c)
    real(8), intent(in) :: xs(:), x0
    real(8), intent(out) :: c(size(xs), 0:2)
    real(8) :: c1, c2, c3, c4, c5
    integer :: i, j, k, n, mn

    n = size(xs)
    c = 0.0d0
    c(1, 0) = 1.0d0
    c1 = 1.0d0
    c4 = xs(1) - x0
    do i = 2, n
      mn = min(i - 1, 2)
      c2 = 1.0d0
      c5 = c4
      c4 = xs(i) - x0
      do j = 1, i - 1
        c3 = xs(i) - xs(j)
        c2 = c2*c3
        if (j == i - 1) then
          do k = mn, 1, -1
            c(i, k) = c1*(k*c(i - 1, k - 1) - c5*c(i - 1, k))/c2
          end do
          c(i, 0) = -c1*c5*c(i - 1, 0)/c2
        end if
        do k = mn, 1, -1
          c(j, k) = (c4*c(j, k) - k*c(j, k - 1))/c3
        end do
        c(j, 0) = c4*c(j, 0)/c3
      end do
      c1 = c2
    end do
  end subroutine derivative_weights

end module COCAL_ID_holes
