! The COCAL ID in memory and the values it gives at one Cactus point.
!
! Every reader stores the same metric fields on a metric domain and its
! matter fields on a matter domain (cocal_id_domain.f90). COCAL_ID_eval_point
! (cocal_id_eval.f90) turns them into one state vector at a point, which one
! routine writes to ADMBase and HydroBase (COCAL_ID_master.F90).

module COCAL_ID_state
  use COCAL_ID_domain
  use COCAL_ID_node_shared, only : COCAL_ID_node_free_all
  implicit none


  ! COCAL's metric fields, in this order for every ID: conformal factor psi,
  ! lapse alpha, shift beta^i, conformal curvature Kt_ij = K_ij/psi**4 in
  ! COCAL units and, for waveless data, h_ij = gt_ij - delta_ij. Symmetric
  ! tensors are ordered xx, xy, xz, yy, yz, zz. A reader may append fields of
  ! its own (MRNS electromagnetic fields).
  integer, parameter :: M_PSI = 1, M_ALP = 2, M_BETA = 3, M_K = 6, M_H = 12
  integer, parameter :: NMETRIC_CF = 11, NMETRIC_WL = 17

  ! Matter fields of a rotating star (RNS, MRNS): Emden function
  ! emd = P/rho, angular velocity Omega and, for MRNS, the 4-velocity u^mu.
  integer, parameter :: MAT_EMD = 1, MAT_OMEGA = 2, MAT_U = 3

  ! ADMBase and HydroBase values at a point, in Cactus units.
  integer, parameter :: S_ALP = 1, S_BETA = 2, S_G = 5, S_K = 11
  integer, parameter :: S_RHO = 17, S_PRESS = 18, S_EPS = 19, S_VEL = 20
  integer, parameter :: S_AVEC = 23, S_APHI = 26, S_BVEC = 27, NSTATE = 29

  integer, parameter :: COCAL_RNS = 1, COCAL_BNS = 2

  ! COCAL's vacuum threshold for the Emden function.
  real(8), parameter :: EMD_VACUUM = 1.0d-14

  ! The ID: its kind, metric and matter domains, the COCAL length unit radi
  ! in Cactus units, whether the metric has h_ij, and the EOS.
  integer, save :: cocal_kind = 0
  type(cocal_domain), target, save :: cocal_metric, cocal_matter
  real(8), save :: cocal_radi = 1.0d0
  logical, save :: cocal_wl = .false.
  logical, save :: cocal_tabulated_eos = .false.

contains

  ! ADM state from the metric fields q.
  subroutine COCAL_ID_metric_to_state(q, s)
    real(8), intent(in) :: q(:)
    real(8), intent(out) :: s(NSTATE)
    real(8), parameter :: delta(6) = [1.0d0, 0.0d0, 0.0d0, 1.0d0, 0.0d0, 1.0d0]
    real(8) :: psi4

    s = 0.0d0
    psi4 = q(M_PSI)**4
    s(S_ALP) = q(M_ALP)
    s(S_BETA:S_BETA + 2) = q(M_BETA:M_BETA + 2)
    s(S_G:S_G + 5) = psi4*delta
    if (cocal_wl) s(S_G:S_G + 5) = s(S_G:S_G + 5) + psi4*q(M_H:M_H + 5)
    s(S_K:S_K + 5) = psi4*q(M_K:M_K + 5)/cocal_radi
  end subroutine COCAL_ID_metric_to_state

  ! Fluid from the Emden function through COCAL's EOS, with Eulerian
  ! velocity v; vacuum at or below EMD_VACUUM.
  subroutine COCAL_ID_set_fluid(emd, v, s)
    real(8), intent(in) :: emd, v(3)
    real(8), intent(inout) :: s(NSTATE)
    real(8) :: q, h, press, rho, ene

    s(S_RHO:S_VEL + 2) = 0.0d0
    if (.not. emd > EMD_VACUUM) return
    q = emd
    if (cocal_tabulated_eos) then
      call COCAL_ID_teos_q2hprho(q, h, press, rho, ene)
    else
      call COCAL_ID_peos_q2hprho(q, h, press, rho, ene)
    end if
    s(S_RHO) = rho
    s(S_PRESS) = press
    if (rho > 0.0d0) s(S_EPS) = ene/rho - 1.0d0
    s(S_VEL:S_VEL + 2) = v
  end subroutine COCAL_ID_set_fluid

  ! Fluid of a rotating star at COCAL point (xc,yc,zc), with metric
  ! fields q there: v^i = (u^i/u^t + beta^i)/alpha, with u^i/u^t = Omega
  ! (-y, x, 0) unless the 4-velocity is exported (MRNS).
  subroutine COCAL_ID_rotating_matter(xc, yc, zc, q, s)
    real(8), intent(in) :: xc, yc, zc, q(:)
    real(8), intent(inout) :: s(NSTATE)
    real(8) :: f(cocal_matter%nvars), dxdt(3)
    logical :: ok

    call COCAL_ID_domain_eval(cocal_matter, xc, yc, zc, f, ok)
    if (.not. ok) return
    if (cocal_matter%nvars > MAT_U) then
      dxdt = f(MAT_U + 1:MAT_U + 3)/f(MAT_U)
    else
      dxdt = f(MAT_OMEGA)*[-yc, xc, 0.0d0]
    end if
    call COCAL_ID_set_fluid(f(MAT_EMD), (dxdt + q(M_BETA:M_BETA + 2))/q(M_ALP), s)
  end subroutine COCAL_ID_rotating_matter

  ! Collective over the node.
  subroutine COCAL_ID_state_free()
    call COCAL_ID_domain_free(cocal_metric)
    call COCAL_ID_domain_free(cocal_matter)
    call COCAL_ID_node_free_all()
    cocal_kind = 0
    cocal_radi = 1.0d0
    cocal_wl = .false.
    cocal_tabulated_eos = .false.
  end subroutine COCAL_ID_state_free

end module COCAL_ID_state
