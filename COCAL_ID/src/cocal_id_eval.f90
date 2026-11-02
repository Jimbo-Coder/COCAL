! The state at one Cactus point, for every ID: the metric at the point, then
! the ID's matter, composed with that metric.

module COCAL_ID_eval
  use COCAL_ID_state
  use COCAL_ID_domain
  use COCAL_ID_rns, only : COCAL_ID_rns_em
  use COCAL_ID_bns, only : COCAL_ID_bns_matter
  use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
  implicit none
  private
  public :: COCAL_ID_eval_point, COCAL_ID_speed2

contains

  ! State s at Cactus point (x,y,z); ok=.false. outside the COCAL data or
  ! for an invalid state (non-finite, nonpositive psi or alpha, v >= 1).
  subroutine COCAL_ID_eval_point(x, y, z, s, ok)
    real(8), intent(in) :: x, y, z
    real(8), intent(out) :: s(NSTATE)
    logical, intent(out) :: ok
    real(8) :: xc, yc, zc, q(cocal_metric%nvars)

    s = 0.0d0
    xc = x/cocal_radi
    yc = y/cocal_radi
    zc = z/cocal_radi
    call COCAL_ID_domain_eval(cocal_metric, xc, yc, zc, q, ok)
    if (.not. ok) return
    call COCAL_ID_metric_to_state(q, s)
    select case (cocal_kind)
    case (COCAL_RNS)
      call COCAL_ID_rotating_matter(xc, yc, zc, q, s)
      call COCAL_ID_rns_em(q, s)
    case (COCAL_BNS)
      call COCAL_ID_bns_matter(xc, yc, zc, q, s)
    end select
    ok = all(ieee_is_finite(s)) .and. COCAL_ID_speed2(s) < 1.0d0
  end subroutine COCAL_ID_eval_point

  ! g_ij v^i v^j.
  pure real(8) function COCAL_ID_speed2(s)
    real(8), intent(in) :: s(NSTATE)

    associate(g => s(S_G:S_G + 5), v => s(S_VEL:S_VEL + 2))
      COCAL_ID_speed2 = g(1)*v(1)**2 + g(4)*v(2)**2 + g(6)*v(3)**2 &
                        + 2.0d0*(g(2)*v(1)*v(2) + g(3)*v(1)*v(3) + g(5)*v(2)*v(3))
    end associate
  end function COCAL_ID_speed2

end module COCAL_ID_eval
