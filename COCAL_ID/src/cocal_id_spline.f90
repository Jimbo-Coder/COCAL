! C2 interpolant of COCAL fields on one spherical patch. The 4-point Lagrange
! stencil is only C0: its slope jumps at every COCAL node, and Cartesian
! second derivatives (the constraints) then grow like 1/h once h is below the
! COCAL spacing. A tensor-product cubic spline is C2 with the same
! fourth-order accuracy: a clamped spline in r (end slopes from the
! four-point one-sided derivative; through the centre along whole diameters
! when the grid starts at r=0) and periodic uniform B-splines in theta
! (continued over the poles to phi+pi) and phi. Coefficients take one extra
! copy of the patch fields.
!
! Where a field's angular data is not smooth (COCAL's exported Cartesian
! gradients keep a finite m=2 harmonic in phi up to the grid pole, while the
! pole row is single valued), the spline rings over several cells and the
! local 4-point Lagrange stencil does not. The spline's third derivative
! jumps at each node by the 4th difference of its coefficients; for smooth
! data that is small next to the 2nd differences, (k h)**2 of them for a
! harmonic of wavenumber k, and comparable to them at a kink or a bad node.
! A node is marked rough, per field, where in theta or phi it exceeds
! ROUGH_RATIO times the 2nd differences (plus ROUGH_FLOOR of the field's
! magnitude) and is isolated, ROUGH_ISOLATION times the 4th differences three
! nodes away (an under-resolved smooth field is large there too). Within two
! nodes of a mark the field blends to the Lagrange stencil, with the C2
! weight of the marks taken as B-spline coefficients (between radial shells,
! a smoothstep). Smooth fields are marked nowhere and keep the spline. The
! domain (cocal_id_domain.f90) applies this to matter fields only. Exports
! whose gradients COCAL takes at fourth order (flgrad_4th_gridpoint in
! calc_gradvep) are smooth there, so they keep the spline everywhere.
module COCAL_ID_spline
  implicit none
  private
  public :: COCAL_ID_spline_fit, COCAL_ID_spline_eval, COCAL_ID_spline_grids_ok
  public :: COCAL_ID_spline_angular, COCAL_ID_spline_eval_angular, COCAL_ID_spline_rough

  real(8), parameter :: ROUGH_RATIO = 2.0d0, ROUGH_FLOOR = 1.0d-3, ROUGH_ISOLATION = 6.0d0
  integer, parameter :: ROUGH_MIN_NODES = 24

contains

  ! Second derivatives of the clamped cubic spline through y(:,0:n) at rg.
  subroutine COCAL_ID_spline_radial_m(rg, y)
    real(8), intent(in) :: rg(0:)
    real(8), intent(inout) :: y(:, 0:)
    integer :: n, i
    real(8), allocatable :: h(:), d(:), cc(:), r(:, :)
    real(8) :: w(4), f
    n = ubound(rg, 1)
    allocate (h(0:n - 1), d(0:n), cc(0:n), r(size(y, 1), 0:n))
    h = rg(1:n) - rg(0:n - 1)
    call COCAL_ID_spline_dlag(rg(0:3), rg(0), w)
    r(:, 0) = 6*((y(:, 1) - y(:, 0))/h(0) - matmul(y(:, 0:3), w))
    call COCAL_ID_spline_dlag(rg(n - 3:n), rg(n), w)
    r(:, n) = 6*(matmul(y(:, n - 3:n), w) - (y(:, n) - y(:, n - 1))/h(n - 1))
    do i = 1, n - 1
      r(:, i) = 6*((y(:, i + 1) - y(:, i))/h(i) - (y(:, i) - y(:, i - 1))/h(i - 1))
    end do
    ! Thomas: rows i have sub h(i-1), diag 2(h(i-1)+h(i)), super h(i).
    d(0) = 2*h(0)
    cc(0) = h(0)
    do i = 1, n
      if (i < n) then
        f = h(i - 1)/d(i - 1)
        d(i) = 2*(h(i - 1) + h(i)) - f*cc(i - 1)
        cc(i) = h(i)
      else
        f = h(n - 1)/d(n - 1)
        d(n) = 2*h(n - 1) - f*cc(n - 1)
      end if
      r(:, i) = r(:, i) - f*r(:, i - 1)
    end do
    y(:, n) = r(:, n)/d(n)
    do i = n - 1, 0, -1
      y(:, i) = (r(:, i) - cc(i)*y(:, i + 1))/d(i)
    end do
  end subroutine COCAL_ID_spline_radial_m
  ! Derivative weights of the 4-point Lagrange polynomial through x at v.
  subroutine COCAL_ID_spline_dlag(x, v, w)
    real(8), intent(in) :: x(4), v
    real(8), intent(out) :: w(4)
    integer :: a, b, c
    real(8) :: term
    w = 0
    do a = 1, 4
      do b = 1, 4
        if (b == a) cycle
        term = 1/(x(a) - x(b))
        do c = 1, 4
          if (c /= a .and. c /= b) term = term*(v - x(c))/(x(a) - x(c))
        end do
        w(a) = w(a) + term
      end do
    end do
  end subroutine COCAL_ID_spline_dlag
  ! Replace a(:,:,0:nt,0:np-1) by periodic uniform cubic B-spline coefficients:
  ! along each great circle (theta at phi, then 2pi-theta at phi+pi), then phi.
  subroutine COCAL_ID_spline_angular(a, nr, nt, np)
    integer, intent(in) :: nr, nt, np
    real(8), intent(inout) :: a(:, 0:, 0:, 0:)
    real(8), allocatable :: circ(:, :, :), ring(:, :, :)
    integer :: j, k, kk
    allocate (circ(size(a, 1), 0:nr, 0:2*nt - 1), ring(size(a, 1), 0:nr, 0:np - 1))
    do k = 0, np/2 - 1
      kk = k + np/2
      do j = 0, nt
        circ(:, :, j) = a(:, :, j, k)
      end do
      do j = nt + 1, 2*nt - 1
        circ(:, :, j) = a(:, :, 2*nt - j, kk)
      end do
      call COCAL_ID_spline_periodic(circ, 2*nt)
      do j = 0, nt
        a(:, :, j, k) = circ(:, :, j)
        a(:, :, j, kk) = circ(:, :, modulo(2*nt - j, 2*nt))
      end do
    end do
    do j = 0, nt
      ring = a(:, :, j, :)
      call COCAL_ID_spline_periodic(ring, np)
      a(:, :, j, :) = ring
    end do
  end subroutine COCAL_ID_spline_angular
  ! Solve (c(i-1)+4c(i)+c(i+1))/6 = f(i), periodic in the last index (length n),
  ! for every leading index (cyclic Thomas with a Sherman-Morrison correction).
  subroutine COCAL_ID_spline_periodic(f, n)
    integer, intent(in) :: n
    real(8), intent(inout) :: f(:, :, 0:)
    real(8) :: d(0:n - 1), zz(0:n - 1), fac
    real(8), allocatable :: x(:, :, :), corr(:, :)
    integer :: i
    allocate (x(size(f, 1), size(f, 2), 0:n - 1), corr(size(f, 1), size(f, 2)))
    ! Modified diagonal: gamma=-4 at both ends, so d0=8 and d(n-1)=4+1/4.
    d(0) = 8
    do i = 1, n - 1
      d(i) = 4 - 1/d(i - 1)
      if (i == n - 1) d(i) = 4.25d0 - 1/d(i - 1)
    end do
    x(:, :, 0) = 6*f(:, :, 0)
    zz(0) = -4
    do i = 1, n - 1
      x(:, :, i) = 6*f(:, :, i) - x(:, :, i - 1)/d(i - 1)
      zz(i) = merge(1d0, 0d0, i == n - 1) - zz(i - 1)/d(i - 1)
    end do
    x(:, :, n - 1) = x(:, :, n - 1)/d(n - 1)
    zz(n - 1) = zz(n - 1)/d(n - 1)
    do i = n - 2, 0, -1
      x(:, :, i) = (x(:, :, i) - x(:, :, i + 1))/d(i)
      zz(i) = (zz(i) - zz(i + 1))/d(i)
    end do
    fac = 1 + zz(0) - zz(n - 1)/4
    corr = (x(:, :, 0) - x(:, :, n - 1)/4)/fac
    do i = 0, n - 1
      f(:, :, i) = x(:, :, i) - corr*zz(i)
    end do
  end subroutine COCAL_ID_spline_periodic
  ! Uniform cubic B-spline weights for nodes i-1..i+2 at fraction t of [i,i+1].
  pure subroutine COCAL_ID_spline_bw(t, w)
    real(8), intent(in) :: t
    real(8), intent(out) :: w(4)
    w(1) = (1 - t)**3/6
    w(2) = (3*t**3 - 6*t*t + 4)/6
    w(3) = (-3*t**3 + 3*t*t + 3*t + 1)/6
    w(4) = t**3/6
  end subroutine COCAL_ID_spline_bw

  ! Fit: f(:,0:nr,0:nt,0:np-1) on rg(0:nr) and uniform theta/phi becomes the
  ! angular coefficients of the values; sm those of the radial second
  ! derivatives. Callers fill unusable nodes first.
  subroutine COCAL_ID_spline_fit(rg, nt, np, f, sm, correct_centre)
    real(8), intent(in) :: rg(0:)
    integer, intent(in) :: nt, np
    real(8), intent(inout) :: f(:, 0:, 0:, 0:)
    real(8), intent(out) :: sm(:, 0:, 0:, 0:)
    real(8), allocatable :: line(:, :), rr(:)
    integer :: it, ip, nr, ia, ipa
    logical, optional, intent(in) :: correct_centre
    logical :: repair_centre
    repair_centre = .true.
    if (present(correct_centre)) repair_centre = correct_centre
    nr = ubound(rg, 1)
    if (rg(0) == 0) then
      ! Grid centred on the origin: fit each diameter, ray (theta,phi) joined to
      ! (pi-theta,phi+pi) at r=0, so the interpolant is C2 through the centre
      ! rather than two independently clamped rays meeting there.
      if (repair_centre) call COCAL_ID_spline_centre_value(rg, nt, np, f)
      allocate (line(size(f, 1), 0:2*nr), rr(0:2*nr))
      rr(0:nr - 1) = -rg(nr:1:-1)
      rr(nr:2*nr) = rg(0:nr)
      do ip = 0, np/2 - 1
        do it = 0, nt
          ia = nt - it
          ipa = ip + np/2
          line(:, 0:nr - 1) = f(:, nr:1:-1, ia, ipa)
          line(:, nr:2*nr) = f(:, 0:nr, it, ip)
          call COCAL_ID_spline_radial_m(rr, line)
          sm(:, 0:nr, it, ip) = line(:, nr:2*nr)
          sm(:, 1:nr, ia, ipa) = line(:, nr - 1:0:-1)
          sm(:, 0, ia, ipa) = line(:, nr)
        end do
      end do
    else
      allocate (line(size(f, 1), 0:nr))
      do ip = 0, np - 1
        do it = 0, nt
          line = f(:, :, it, ip)
          call COCAL_ID_spline_radial_m(rg, line)
          sm(:, :, it, ip) = line
        end do
      end do
    end if
    call COCAL_ID_spline_angular(f, nr, nt, np)
    call COCAL_ID_spline_angular(sm, nr, nt, np)
  end subroutine COCAL_ID_spline_fit

  ! COCAL's exported value at r=0 is offset from the smooth continuation of
  ! its neighbours (~1e-5 in psi and alph for the RNS and BNS exports; every
  ! other node lies on the smooth curve to ~1e-7), which puts a spike of
  ! ~16 pi rho into the Hamiltonian constraint at each star centre. Replace it
  ! by the even extrapolation from nodes 1..4 (4-point Lagrange in r**2),
  ! averaged over all diameters so the centre keeps a single value.
  subroutine COCAL_ID_spline_centre_value(rg, nt, np, f)
    real(8), intent(in) :: rg(0:)
    integer, intent(in) :: nt, np
    real(8), intent(inout) :: f(:, 0:, 0:, 0:)
    real(8) :: w(4), centre(size(f, 1))
    integer :: it, ip, ia, ipa, i
    call COCAL_ID_lagint_4th_weights(rg(1:4)**2, 0.0d0, w)
    centre = 0
    do ip = 0, np/2 - 1
      do it = 0, nt
        ia = nt - it
        ipa = ip + np/2
        do i = 1, 4
          centre = centre + w(i)*0.5d0*(f(:, i, it, ip) + f(:, i, ia, ipa))
        end do
      end do
    end do
    centre = centre/((np/2)*(nt + 1))
    do ip = 0, np - 1
      do it = 0, nt
        f(:, 0, it, ip) = centre
      end do
    end do
  end subroutine COCAL_ID_spline_centre_value

  ! Node (j,k) of the angular grid for any j: continued over the poles at phi+pi.
  pure subroutine wrap(nt, np, j, k)
    integer, intent(in) :: nt, np
    integer, intent(inout) :: j, k
    if (j < 0 .or. j > nt) then
      j = merge(-j, 2*nt - j, j < 0)
      k = k + np/2
    end if
    k = modulo(k, np)
  end subroutine wrap

  ! Marks of rough angular data (see above), from the fitted coefficients
  ! c(1:nvars,0:nr,0:nt,0:np-1) and each field's magnitude scale(1:nvars):
  ! rough(ir,it,ip) is the bitmask of the fields marked at the node, held as
  ! an exact integer. Each shell is copied into e, padded by 5 nodes and
  ! continued over the poles; most nodes fail the first test and go no further.
  subroutine COCAL_ID_spline_rough(nt, np, c, scale, rough)
    integer, intent(in) :: nt, np
    real(8), intent(in) :: c(:, 0:, 0:, 0:), scale(:)
    real(8), intent(out) :: rough(0:, 0:, 0:)
    real(8), allocatable :: e(:, :), me(:, :), d4(:, :)
    integer(8), allocatable :: mark(:, :)
    real(8) :: q, d2
    integer :: ir, iv, it, ip, dir, di, dk, a, b

    rough = 0
    ! The test reads 11 nodes along a line: on coarser grids it wraps onto itself.
    if (nt < ROUGH_MIN_NODES .or. np < ROUGH_MIN_NODES) return
    allocate (e(-5:nt + 5, -5:np + 4), me(-2:nt + 2, -2:np + 1), d4(-3:nt + 3, -3:np + 2), mark(0:nt, 0:np - 1))
    do ir = 0, ubound(c, 2)
      mark = 0
      do iv = 1, size(c, 1)
        e(0:nt, 0:np - 1) = c(iv, ir, :, :)
        call pad(e, 5)
        do dir = 1, 2
          di = merge(1, 0, dir == 1)
          dk = 1 - di
          ! |4th differences| out to three nodes beyond the grid, for the
          ! isolation test.
          d4 = abs(e(-3 - 2*di:nt + 3 - 2*di, -3 - 2*dk:np + 2 - 2*dk) - 4*e(-3 - di:nt + 3 - di, -3 - dk:np + 2 - dk) &
                   + 6*e(-3:nt + 3, -3:np + 2) - 4*e(-3 + di:nt + 3 + di, -3 + dk:np + 2 + dk) &
                   + e(-3 + 2*di:nt + 3 + 2*di, -3 + 2*dk:np + 2 + 2*dk))
          do ip = 0, np - 1
            do it = 0, nt
              q = d4(it, ip)
              if (q <= ROUGH_RATIO*ROUGH_FLOOR*scale(iv)) cycle
              if (q <= ROUGH_ISOLATION*max(d4(it + 3*di, ip + 3*dk), d4(it - 3*di, ip - 3*dk))) cycle
              d2 = max(abs(v(-2) - 2*v(-1) + v(0)), abs(v(-1) - 2*v(0) + v(1)), abs(v(0) - 2*v(1) + v(2)))
              if (q > ROUGH_RATIO*(d2 + ROUGH_FLOOR*scale(iv))) mark(it, ip) = ibset(mark(it, ip), iv - 1)
            end do
          end do
        end do
      end do
      ! Marks spread to the nodes up to two away, so that the blend weight is
      ! 1 on and next to a marked node.
      me(0:nt, 0:np - 1) = real(mark, 8)
      call pad(me, 2)
      mark = 0
      do b = -2, 2
        do a = -2, 2
          mark = ior(mark, int(me(a:nt + a, b:np - 1 + b), 8))
        end do
      end do
      rough(ir, :, :) = real(mark, 8)
    end do

  contains
    ! Fill the n-node border of f(-n:nt+n, -n:np-1+n) from f(0:nt, 0:np-1):
    ! over the poles at phi+pi, then periodic in phi.
    subroutine pad(f, n)
      real(8), intent(inout) :: f(-n:, -n:)
      integer, intent(in) :: n
      integer :: j
      do j = 1, n
        f(-j, 0:np - 1) = cshift(f(j, 0:np - 1), np/2)
        f(nt + j, 0:np - 1) = cshift(f(nt - j, 0:np - 1), np/2)
      end do
      f(:, -n:-1) = f(:, np - n:np - 1)
      f(:, np:np - 1 + n) = f(:, 0:n - 1)
    end subroutine pad

    ! e at s nodes from (it,ip) along dir.
    real(8) function v(s)
      integer, intent(in) :: s
      v = e(it + s*di, ip + s*dk)
    end function v
  end subroutine COCAL_ID_spline_rough

  ! Angular stencil at the direction of (xc,yc,zc): the 6x6 coefficient
  ! window (j(jj),k(kk)) around the cell, continued over the poles, with the
  ! B-spline weights w (outer ring zero) and the Lagrange weights wl of the
  ! 4x4 nodes around the cell (a node value is (c(n-1)+4c(n)+c(n+1))/6 of the
  ! coefficients). it and ip are the cell's lower theta and phi nodes, and tt
  ! and tp the fractions of the cell.
  pure subroutine angular_stencil(nt, np, xc, yc, zc, j, k, w, wl, it, ip, tt, tp)
    integer, intent(in) :: nt, np
    real(8), intent(in) :: xc, yc, zc
    integer, intent(out) :: j(6, 6), k(6, 6), it, ip
    real(8), intent(out) :: w(6, 6), wl(6, 6), tt, tp
    real(8), parameter :: pi = 3.14159265358979323846264338327950d0
    real(8) :: u, wt(6), wp(6), lt(6), lp(6)
    integer :: jj, kk

    u = atan2(sqrt(xc*xc + yc*yc), zc)/(pi/nt)
    it = min(int(u), nt - 1)
    tt = u - it
    u = modulo(atan2(yc, xc), 2*pi)/(2*pi/np)
    ip = min(int(u), np - 1)
    tp = u - ip
    wt = 0
    wp = 0
    call COCAL_ID_spline_bw(tt, wt(2:5))
    call COCAL_ID_spline_bw(tp, wp(2:5))
    call node_lagrange(tt, lt)
    call node_lagrange(tp, lp)
    do kk = 1, 6
      do jj = 1, 6
        w(jj, kk) = wt(jj)*wp(kk)
        wl(jj, kk) = lt(jj)*lp(kk)
        j(jj, kk) = it - 3 + jj
        k(jj, kk) = ip - 3 + kk
        call wrap(nt, np, j(jj, kk), k(jj, kk))
      end do
    end do

  contains
    ! 4-point Lagrange weights at fraction t for nodes -1..2, as weights on
    ! coefficients -2..3.
    pure subroutine node_lagrange(t, c)
      real(8), intent(in) :: t
      real(8), intent(out) :: c(6)
      real(8) :: l(4)
      integer :: a
      l(1) = -t*(t - 1)*(t - 2)/6
      l(2) = (t + 1)*(t - 1)*(t - 2)/2
      l(3) = -(t + 1)*t*(t - 2)/2
      l(4) = (t + 1)*t*(t - 1)/6
      c = 0
      do a = 1, 4
        c(a:a + 2) = c(a:a + 2) + l(a)*[1, 4, 1]/6.0d0
      end do
    end subroutine node_lagrange
  end subroutine angular_stencil

  ! C2 value at patch coordinates (xc,yc,zc); ok=.false. outside rg's range.
  ! cell and frac, if present, return the lower (r, theta, phi) nodes of the
  ! cell holding the point and the point's fractions of it. rough, if
  ! present, holds the marks of COCAL_ID_spline_rough.
  subroutine COCAL_ID_spline_eval(rg, nt, np, sf, sm, xc, yc, zc, q, ok, cell, frac, rough)
    real(8), intent(in) :: rg(0:), sf(:, 0:, 0:, 0:), sm(:, 0:, 0:, 0:), xc, yc, zc
    integer, intent(in) :: nt, np
    real(8), intent(out) :: q(size(sf, 1))
    logical, intent(out) :: ok
    integer, intent(out), optional :: cell(3)
    real(8), intent(out), optional :: frac(3)
    real(8), intent(in), optional :: rough(0:, 0:, 0:)
    real(8) :: rc, h, a, b, c, d, w(6, 6), wl(6, 6), tt, tp, sb, s, ql(size(sf, 1)), f(size(sf, 1))
    integer :: nr, i, it, ip, jj, kk, iv, j(6, 6), k(6, 6), lo, hi, mid, n0, n1
    integer(8) :: m0(2:5, 2:5), m1(2:5, 2:5)
    logical :: local
    q = 0
    nr = ubound(rg, 1)
    rc = sqrt(xc*xc + yc*yc + zc*zc)
    ok = rc >= rg(0) - 64*epsilon(rc)*max(1d0, rc) .and. rc <= rg(nr) + 64*epsilon(rc)*max(1d0, rc)
    if (.not. ok) return
    rc = max(rg(0), min(rc, rg(nr)))
    lo = 0
    hi = nr
    do while (hi - lo > 1)
      mid = (lo + hi)/2
      if (rg(mid) <= rc) then
        lo = mid
      else
        hi = mid
      end if
    end do
    i = lo
    h = rg(i + 1) - rg(i)
    a = (rg(i + 1) - rc)/h
    b = 1 - a
    c = (a**3 - a)*h*h/6
    d = (b**3 - b)*h*h/6
    call angular_stencil(nt, np, xc, yc, zc, j, k, w, wl, it, ip, tt, tp)
    local = .false.
    if (present(rough)) then
      do kk = 2, 5
        do jj = 2, 5
          m0(jj, kk) = int(rough(i, j(jj, kk), k(jj, kk)), 8)
          m1(jj, kk) = int(rough(i + 1, j(jj, kk), k(jj, kk)), 8)
        end do
      end do
      local = any(m0 /= 0) .or. any(m1 /= 0)
    end if
    ! The B-spline needs the inner 4x4 coefficients, the Lagrange stencil all 6x6.
    n0 = merge(1, 2, local)
    n1 = merge(6, 5, local)
    ql = 0
    do kk = n0, n1
      do jj = n0, n1
        f = a*sf(:, i, j(jj, kk), k(jj, kk)) + b*sf(:, i + 1, j(jj, kk), k(jj, kk)) &
            + c*sm(:, i, j(jj, kk), k(jj, kk)) + d*sm(:, i + 1, j(jj, kk), k(jj, kk))
        q = q + w(jj, kk)*f
        if (local) ql = ql + wl(jj, kk)*f
      end do
    end do
    if (local) then
      ! Blend weight of field iv: its marks as B-spline coefficients in angle,
      ! a smoothstep between the shells.
      sb = b**3*(10 - 15*b + 6*b*b)
      do iv = 1, size(q)
        s = 0
        do kk = 2, 5
          do jj = 2, 5
            if (btest(m0(jj, kk), iv - 1)) s = s + (1 - sb)*w(jj, kk)
            if (btest(m1(jj, kk), iv - 1)) s = s + sb*w(jj, kk)
          end do
        end do
        q(iv) = q(iv) + s*(ql(iv) - q(iv))
      end do
    end if
    if (present(cell)) cell = [i, it, ip]
    if (present(frac)) frac = [b, tt, tp]
  end subroutine COCAL_ID_spline_eval

  ! Value of an angular-only field (a stellar surface) in the direction of
  ! (xc,yc,zc), from coefficients a(1,0,0:nt,0:np-1) of COCAL_ID_spline_angular.
  pure real(8) function COCAL_ID_spline_eval_angular(nt, np, a, xc, yc, zc) result(v)
    integer, intent(in) :: nt, np
    real(8), intent(in) :: a(:, 0:, 0:, 0:), xc, yc, zc
    real(8) :: w(6, 6), wl(6, 6), tt, tp
    integer :: j(6, 6), k(6, 6), it, ip, jj, kk
    call angular_stencil(nt, np, xc, yc, zc, j, k, w, wl, it, ip, tt, tp)
    v = 0
    do kk = 2, 5
      do jj = 2, 5
        v = v + w(jj, kk)*a(1, 0, j(jj, kk), k(jj, kk))
      end do
    end do
  end function COCAL_ID_spline_eval_angular

  ! True when thg(0:nt) and phig(0:np) are the uniform grids the fit assumes.
  pure logical function COCAL_ID_spline_grids_ok(thg, phig)
    real(8), intent(in) :: thg(0:), phig(0:)
    real(8), parameter :: pi = 3.14159265358979323846264338327950d0
    integer :: i, nt, np
    nt = ubound(thg, 1)
    np = ubound(phig, 1)
    COCAL_ID_spline_grids_ok = mod(np, 2) == 0 .and. &
                               all(abs(thg - [(i*pi/nt, i=0, nt)]) < 1d-9) .and. all(abs(phig - [(i*2*pi/np, i=0, np)]) < 1d-9)
  end function COCAL_ID_spline_grids_ok

end module COCAL_ID_spline
