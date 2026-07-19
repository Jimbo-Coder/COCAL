! COCAL_ID BHT fill/smoothing routines.
!
! Provenance:
!   - COCAL_ID_bht_read_IDpar_cenvals preserves the legacy read_bin_AH.par
!     center-anchor input used by the old AH filling workflow.
!   - The other routines below are COCAL_ID additions around the BHT data
!     already loaded by COCAL_ID_read_bht_data. They are not copied as named
!     upstream COCAL routines.

! COCAL_ID-only reset for the fill machinery. Mode-specific values are either
! read from read_bin_AH.par or computed by the smooth-fill path.
subroutine COCAL_ID_bht_reset_fill_state()

  use COCAL_ID_data_bht
  implicit none

    rexc = 0.0d0; r_ah = 0.0d0; r_eh = 0.0d0; r_ring = 0.0d0
    xx_bh = 0.0d0; yy_bh = 0.0d0; zz_bh = 0.0d0
    extra_points_dr = 0.0d0
    trk_xi = 0.0d0
    trk_npow = 0.0d0

    lap_cent = 0.0d0
    betax_cent = 0.0d0; betay_cent = 0.0d0; betaz_cent = 0.0d0

    gxx_cent =  0.0d0
    gxy_cent =  0.0d0
    gxz_cent =  0.0d0
    gyy_cent =  0.0d0
    gyz_cent =  0.0d0
    gzz_cent =  0.0d0

    kxx_cent =  0.0d0
    kxy_cent =  0.0d0
    kxz_cent =   0.0d0
    kyy_cent =   0.0d0
    kyz_cent =   0.0d0
    kzz_cent =   0.0d0

    rho_cent = 0.0d0
    vx_cent = 0.0d0; vy_cent = 0.0d0; vz_cent = 0.0d0
    smooth_center = 0.0d0
    smooth_center_attempted = .false.
    smooth_center_valid = .false.
end subroutine COCAL_ID_bht_reset_fill_state

subroutine COCAL_ID_bht_read_IDpar_cenvals(dir_path, loaded)

  use COCAL_ID_data_bht
  implicit none
    character(len=*), intent(in) :: dir_path
    logical, intent(out) :: loaded
    character(len=500) :: filenm
    logical :: exists
    integer :: unitno, extra_points
    real(8) :: psi_cent_dummy

    filenm = trim(dir_path)//"/read_bin_AH.par"
    inquire(file=trim(filenm), exist=exists)
    loaded = .false.
    if (.not. exists) return

    unitno = 99
    open(unitno, file=trim(filenm), status='old')
    read(unitno,*) xx_bh
    read(unitno,*) yy_bh
    read(unitno,*) zz_bh
    read(unitno,*) r_ah
    read(unitno,*) extra_points
    read(unitno,*) extra_points_dr
    read(unitno,*) lap_cent
    read(unitno,*) psi_cent_dummy
    read(unitno,*) betax_cent
    read(unitno,*) betay_cent
    read(unitno,*) betaz_cent
    read(unitno,*) gxx_cent
    read(unitno,*) gxy_cent
    read(unitno,*) gxz_cent
    read(unitno,*) gyy_cent
    read(unitno,*) gyz_cent
    read(unitno,*) gzz_cent
    read(unitno,*) kxx_cent
    read(unitno,*) kxy_cent
    read(unitno,*) kxz_cent
    read(unitno,*) kyz_cent
    read(unitno,*) kzz_cent
    read(unitno,*) rho_cent
    read(unitno,*) vx_cent
    read(unitno,*) vy_cent
    read(unitno,*) vz_cent
    read(unitno,*) trk_xi
    read(unitno,*) trk_npow
    close(unitno)
    loaded = .true.
end subroutine COCAL_ID_bht_read_IDpar_cenvals

! COCAL_ID-only automatic scale choice for smooth filling. This replaces the
! legacy trial-and-error choice of read_bin_AH.par spacing/trace damping values.
subroutine COCAL_ID_bht_set_smooth_scales()

  use COCAL_ID_data_bht
  use COCAL_ID_grid_parameter
  use COCAL_ID_coordinate_grav_r
  implicit none
    integer :: ir
    real(8) :: dr_local, gap, rah_safe

    rah_safe = max(r_ah, 1.0d-12)
    dr_local = 0.0d0
    do ir = 0, nrg-1
       if (rg(ir) <= r_ah .and. rg(ir+1) > r_ah) then
          dr_local = rg(ir+1) - rg(ir)
          exit
       end if
    end do
    if (dr_local <= 0.0d0) dr_local = 0.05d0*rah_safe

    extra_points_dr = max(2.0d0*dr_local, 0.20d0*rah_safe)
    extra_points_dr = min(extra_points_dr, 0.25d0*rah_safe)

    gap = max(0.0d0, r_ah - rexc)
    if (gap > 0.0d0) then
       trk_xi = max(0.04d0, min(0.16d0, 1.25d0*gap/rah_safe))
    else
       trk_xi = 0.08d0
    end if
    trk_npow = 2.0d0
end subroutine COCAL_ID_bht_set_smooth_scales

! Internal point-evaluation helper extracted from the exterior interpolation
! block of PCOCAL's coc2pri_bht_wl. This is not a named upstream routine.
subroutine COCAL_ID_bht_eval_raw(xcac, ycac, zcac, allow_inside, ok, bht_point_vals)

  use COCAL_ID_grid_parameter, eps_cocal => eps
  use COCAL_ID_interface_modules_cartesian
  use COCAL_ID_data_bht
  implicit none
    real(8), intent(in) :: xcac, ycac, zcac
    logical, intent(in) :: allow_inside
    logical, intent(out) :: ok
    real(8), intent(out) :: bht_point_vals(bht_nvars)
    integer :: gr_irgex4(4,4,4), gr_itgex4(4,4,4), gr_ipgex4(4,4,4)
    real(8) :: gr_wr(4), gr_wth(4), gr_wphi(4)
    real(8) :: xcoc, ycoc, zcoc, rcoc
    real(8) :: psica, alphca, psi4ca
    real(8) :: bvxdca, bvydca, bvzdca, bvxuca, bvyuca, bvzuca
    real(8) :: hxxdca, hxydca, hxzdca, hyydca, hyzdca, hzzdca
    real(8) :: kxxca, kxyca, kxzca, kyyca, kyzca, kzzca
    real(8) :: emdgca, omegca, hca, preca, rhoca, eneca
    real(8) :: beta_x, beta_y, beta_z, bxcor, bycor, bzcor

    bht_point_vals = 0.0d0
    ok = .false.

    xcoc = xcac/radi
    ycoc = ycac/radi
    zcoc = zcac/radi
    rcoc = dsqrt(xcoc*xcoc + ycoc*ycoc + zcoc*zcoc)
    if ((.not. allow_inside) .and. (rcoc < rexc)) return

    call COCAL_ID_gr2cgr_4th_setup(xcoc, ycoc, zcoc, gr_irgex4, gr_itgex4, gr_ipgex4, gr_wr, gr_wth, gr_wphi)
    call COCAL_ID_cgr_4th_apply(psi,  psica,  gr_irgex4, gr_itgex4, gr_ipgex4, gr_wr, gr_wth, gr_wphi)
    call COCAL_ID_cgr_4th_apply(alph, alphca, gr_irgex4, gr_itgex4, gr_ipgex4, gr_wr, gr_wth, gr_wphi)
    call COCAL_ID_cgr_4th_apply(bvxd, bvxdca, gr_irgex4, gr_itgex4, gr_ipgex4, gr_wr, gr_wth, gr_wphi)
    call COCAL_ID_cgr_4th_apply(bvyd, bvydca, gr_irgex4, gr_itgex4, gr_ipgex4, gr_wr, gr_wth, gr_wphi)
    call COCAL_ID_cgr_4th_apply(bvzd, bvzdca, gr_irgex4, gr_itgex4, gr_ipgex4, gr_wr, gr_wth, gr_wphi)
    call COCAL_ID_cgr_4th_apply(bvxu, bvxuca, gr_irgex4, gr_itgex4, gr_ipgex4, gr_wr, gr_wth, gr_wphi)
    call COCAL_ID_cgr_4th_apply(bvyu, bvyuca, gr_irgex4, gr_itgex4, gr_ipgex4, gr_wr, gr_wth, gr_wphi)
    call COCAL_ID_cgr_4th_apply(bvzu, bvzuca, gr_irgex4, gr_itgex4, gr_ipgex4, gr_wr, gr_wth, gr_wphi)

    call COCAL_ID_cgr_4th_apply(hxxd, hxxdca, gr_irgex4, gr_itgex4, gr_ipgex4, gr_wr, gr_wth, gr_wphi)
    call COCAL_ID_cgr_4th_apply(hxyd, hxydca, gr_irgex4, gr_itgex4, gr_ipgex4, gr_wr, gr_wth, gr_wphi)
    call COCAL_ID_cgr_4th_apply(hxzd, hxzdca, gr_irgex4, gr_itgex4, gr_ipgex4, gr_wr, gr_wth, gr_wphi)
    call COCAL_ID_cgr_4th_apply(hyyd, hyydca, gr_irgex4, gr_itgex4, gr_ipgex4, gr_wr, gr_wth, gr_wphi)
    call COCAL_ID_cgr_4th_apply(hyzd, hyzdca, gr_irgex4, gr_itgex4, gr_ipgex4, gr_wr, gr_wth, gr_wphi)
    call COCAL_ID_cgr_4th_apply(hzzd, hzzdca, gr_irgex4, gr_itgex4, gr_ipgex4, gr_wr, gr_wth, gr_wphi)

    call COCAL_ID_cgr_4th_apply(kxxa, kxxca, gr_irgex4, gr_itgex4, gr_ipgex4, gr_wr, gr_wth, gr_wphi)
    call COCAL_ID_cgr_4th_apply(kxya, kxyca, gr_irgex4, gr_itgex4, gr_ipgex4, gr_wr, gr_wth, gr_wphi)
    call COCAL_ID_cgr_4th_apply(kxza, kxzca, gr_irgex4, gr_itgex4, gr_ipgex4, gr_wr, gr_wth, gr_wphi)
    call COCAL_ID_cgr_4th_apply(kyya, kyyca, gr_irgex4, gr_itgex4, gr_ipgex4, gr_wr, gr_wth, gr_wphi)
    call COCAL_ID_cgr_4th_apply(kyza, kyzca, gr_irgex4, gr_itgex4, gr_ipgex4, gr_wr, gr_wth, gr_wphi)
    call COCAL_ID_cgr_4th_apply(kzza, kzzca, gr_irgex4, gr_itgex4, gr_ipgex4, gr_wr, gr_wth, gr_wphi)

    call COCAL_ID_cgr_4th_apply(emdg, emdgca, gr_irgex4, gr_itgex4, gr_ipgex4, gr_wr, gr_wth, gr_wphi)
    call COCAL_ID_cgr_4th_apply(omeg, omegca, gr_irgex4, gr_itgex4, gr_ipgex4, gr_wr, gr_wth, gr_wphi)

    psi4ca = psica**4
    bht_point_vals(BHT_GXX) = psi4ca*(1.0d0 + hxxdca)
    bht_point_vals(BHT_GXY) = psi4ca*hxydca
    bht_point_vals(BHT_GXZ) = psi4ca*hxzdca
    bht_point_vals(BHT_GYY) = psi4ca*(1.0d0 + hyydca)
    bht_point_vals(BHT_GYZ) = psi4ca*hyzdca
    bht_point_vals(BHT_GZZ) = psi4ca*(1.0d0 + hzzdca)

    bht_point_vals(BHT_KXX) = psi4ca*kxxca/radi
    bht_point_vals(BHT_KXY) = psi4ca*kxyca/radi
    bht_point_vals(BHT_KXZ) = psi4ca*kxzca/radi
    bht_point_vals(BHT_KYY) = psi4ca*kyyca/radi
    bht_point_vals(BHT_KYZ) = psi4ca*kyzca/radi
    bht_point_vals(BHT_KZZ) = psi4ca*kzzca/radi

    beta_x = psi4ca*bvxdca
    beta_y = psi4ca*bvydca
    beta_z = psi4ca*bvzdca
    bht_point_vals(BHT_ALP) = alphca
    bht_point_vals(BHT_BETAX) = beta_x
    bht_point_vals(BHT_BETAY) = beta_y
    bht_point_vals(BHT_BETAZ) = beta_z

    bxcor = bvxuca + omegca*(-ycoc)
    bycor = bvyuca + omegca*xcoc
    bzcor = bvzuca
    if (dabs(emdgca) > 1.0d-14) then
       bht_point_vals(BHT_VX) = bxcor/alphca
       bht_point_vals(BHT_VY) = bycor/alphca
       bht_point_vals(BHT_VZ) = bzcor/alphca
       if (bht_tabulated_eos) then
          call COCAL_ID_teos_q2hprho(emdgca, hca, preca, rhoca, eneca)
       else
          call COCAL_ID_peos_q2hprho(emdgca, hca, preca, rhoca, eneca)
       end if
       bht_point_vals(BHT_RHO) = rhoca
       bht_point_vals(BHT_PRESS) = preca
       if (rhoca > 0.0d0) bht_point_vals(BHT_EPS) = eneca/rhoca - 1.0d0
    end if

    ok = .true.
end subroutine COCAL_ID_bht_eval_raw

subroutine COCAL_ID_bht_metric_det(bht_point_vals, detg)

  use COCAL_ID_data_bht
  implicit none
    real(8), intent(in) :: bht_point_vals(bht_nvars)
    real(8), intent(out) :: detg

    detg = -bht_point_vals(BHT_GXZ)*bht_point_vals(BHT_GXZ)*bht_point_vals(BHT_GYY) &
         & + 2.0d0*bht_point_vals(BHT_GXY)*bht_point_vals(BHT_GXZ)*bht_point_vals(BHT_GYZ) &
         & - bht_point_vals(BHT_GXX)*bht_point_vals(BHT_GYZ)*bht_point_vals(BHT_GYZ) &
         & - bht_point_vals(BHT_GXY)*bht_point_vals(BHT_GXY)*bht_point_vals(BHT_GZZ) &
         & + bht_point_vals(BHT_GXX)*bht_point_vals(BHT_GYY)*bht_point_vals(BHT_GZZ)
end subroutine COCAL_ID_bht_metric_det

subroutine COCAL_ID_bht_metric_inverse(bht_point_vals, detg, gupxx, gupxy, gupxz, gupyy, gupyz, gupzz, ok)

  use COCAL_ID_data_bht
  implicit none
    real(8), intent(in) :: bht_point_vals(bht_nvars)
    real(8), intent(out) :: detg, gupxx, gupxy, gupxz, gupyy, gupyz, gupzz
    logical, intent(out) :: ok

    ok = .false.
    gupxx = 0.0d0; gupxy = 0.0d0; gupxz = 0.0d0
    gupyy = 0.0d0; gupyz = 0.0d0; gupzz = 0.0d0

    call COCAL_ID_bht_metric_det(bht_point_vals, detg)
    if (detg <= 0.0d0) return

    gupxx = (-bht_point_vals(BHT_GYZ)*bht_point_vals(BHT_GYZ) + bht_point_vals(BHT_GYY)*bht_point_vals(BHT_GZZ))/detg
    gupxy =  (bht_point_vals(BHT_GXZ)*bht_point_vals(BHT_GYZ) - bht_point_vals(BHT_GXY)*bht_point_vals(BHT_GZZ))/detg
    gupxz = (-bht_point_vals(BHT_GXZ)*bht_point_vals(BHT_GYY) + bht_point_vals(BHT_GXY)*bht_point_vals(BHT_GYZ))/detg
    gupyy = (-bht_point_vals(BHT_GXZ)*bht_point_vals(BHT_GXZ) + bht_point_vals(BHT_GXX)*bht_point_vals(BHT_GZZ))/detg
    gupyz =  (bht_point_vals(BHT_GXY)*bht_point_vals(BHT_GXZ) - bht_point_vals(BHT_GXX)*bht_point_vals(BHT_GYZ))/detg
    gupzz = (-bht_point_vals(BHT_GXY)*bht_point_vals(BHT_GXY) + bht_point_vals(BHT_GXX)*bht_point_vals(BHT_GYY))/detg
    ok = .true.
end subroutine COCAL_ID_bht_metric_inverse

subroutine COCAL_ID_bht_zero_hydro(bht_point_vals)

  use COCAL_ID_data_bht
  implicit none
    real(8), intent(inout) :: bht_point_vals(bht_nvars)

    bht_point_vals(BHT_RHO) = 0.0d0
    bht_point_vals(BHT_PRESS) = 0.0d0
    bht_point_vals(BHT_EPS) = 0.0d0
    bht_point_vals(BHT_VX) = 0.0d0
    bht_point_vals(BHT_VY) = 0.0d0
    bht_point_vals(BHT_VZ) = 0.0d0
end subroutine COCAL_ID_bht_zero_hydro

! COCAL_ID-only guard used by the fill path to reject NaNs, nonpositive lapse,
! or nonpositive 3-metric determinant before writing Cactus variables.
subroutine COCAL_ID_bht_validate_state(ok, bht_point_vals)

  use COCAL_ID_data_bht
  implicit none
    logical, intent(inout) :: ok
    real(8), intent(in) :: bht_point_vals(bht_nvars)
    real(8) :: detg
    integer :: iv

    if (.not. ok) return

    do iv = 1, bht_nvars
       if (bht_point_vals(iv) /= bht_point_vals(iv)) then
          ok = .false.
          return
       end if
    end do

    call COCAL_ID_bht_metric_det(bht_point_vals, detg)
    if (bht_point_vals(BHT_ALP) <= 0.0d0 .or. detg <= 0.0d0) ok = .false.
end subroutine COCAL_ID_bht_validate_state

! COCAL_ID-only K_ij trace split used by the smooth fill. This is not a PCOCAL
! routine; it keeps the interpolated trace handling explicit and local.
subroutine COCAL_ID_bht_split_curvature(bht_point_vals, avec, trk, ok)

  use COCAL_ID_data_bht
  implicit none
    real(8), intent(in) :: bht_point_vals(bht_nvars)
    real(8), intent(out) :: avec(6), trk
    logical, intent(out) :: ok
    real(8) :: detg, gupxx, gupxy, gupxz, gupyy, gupyz, gupzz

    ok = .false.
    avec = 0.0d0
    trk = 0.0d0

    call COCAL_ID_bht_metric_inverse(bht_point_vals, detg, gupxx, gupxy, gupxz, gupyy, gupyz, gupzz, ok)
    if (.not. ok) return

    trk = gupxx*bht_point_vals(BHT_KXX) + gupyy*bht_point_vals(BHT_KYY) + gupzz*bht_point_vals(BHT_KZZ) &
        & + 2.0d0*(gupxy*bht_point_vals(BHT_KXY) + gupxz*bht_point_vals(BHT_KXZ) + gupyz*bht_point_vals(BHT_KYZ))

    avec(1) = bht_point_vals(BHT_KXX) - bht_point_vals(BHT_GXX)*trk/3.0d0
    avec(2) = bht_point_vals(BHT_KXY) - bht_point_vals(BHT_GXY)*trk/3.0d0
    avec(3) = bht_point_vals(BHT_KXZ) - bht_point_vals(BHT_GXZ)*trk/3.0d0
    avec(4) = bht_point_vals(BHT_KYY) - bht_point_vals(BHT_GYY)*trk/3.0d0
    avec(5) = bht_point_vals(BHT_KYZ) - bht_point_vals(BHT_GYZ)*trk/3.0d0
    avec(6) = bht_point_vals(BHT_KZZ) - bht_point_vals(BHT_GZZ)*trk/3.0d0
    ok = .true.
end subroutine COCAL_ID_bht_split_curvature

! COCAL_ID-only helper for storing tracefree K_ij components in the fill vector.
subroutine COCAL_ID_bht_store_tracefree(bht_point_vals, avec)

  use COCAL_ID_data_bht
  implicit none
    real(8), intent(inout) :: bht_point_vals(bht_nvars)
    real(8), intent(in) :: avec(6)

    bht_point_vals(BHT_KXX) = avec(1)
    bht_point_vals(BHT_KXY) = avec(2)
    bht_point_vals(BHT_KXZ) = avec(3)
    bht_point_vals(BHT_KYY) = avec(4)
    bht_point_vals(BHT_KYZ) = avec(5)
    bht_point_vals(BHT_KZZ) = avec(6)
end subroutine COCAL_ID_bht_store_tracefree

! COCAL_ID-only projection used after smoothing so the interior tracefree part
! is consistent with the filled metric.
subroutine COCAL_ID_bht_project_tracefree(bht_point_vals, ok)

  use COCAL_ID_data_bht
  implicit none
    real(8), intent(inout) :: bht_point_vals(bht_nvars)
    logical, intent(out) :: ok
    real(8) :: detg, gupxx, gupxy, gupxz, gupyy, gupyz, gupzz

    ok = .false.
    call COCAL_ID_bht_metric_inverse(bht_point_vals, detg, gupxx, gupxy, gupxz, gupyy, gupyz, gupzz, ok)
    if (.not. ok) return
    if (dabs(gupyy) <= 1.0d-14) return

    bht_point_vals(BHT_KYY) = -(gupxx*bht_point_vals(BHT_KXX) + gupzz*bht_point_vals(BHT_KZZ) &
         & + 2.0d0*(gupxy*bht_point_vals(BHT_KXY) + gupxz*bht_point_vals(BHT_KXZ) + gupyz*bht_point_vals(BHT_KYZ)))/gupyy
    ok = .true.
end subroutine COCAL_ID_bht_project_tracefree

! COCAL_ID-only helper that restores the controlled trace profile after the
! tracefree projection.
subroutine COCAL_ID_bht_add_trace(bht_point_vals, trk)

  use COCAL_ID_data_bht
  implicit none
    real(8), intent(inout) :: bht_point_vals(bht_nvars)
    real(8), intent(in) :: trk

    bht_point_vals(BHT_KXX) = bht_point_vals(BHT_KXX) + bht_point_vals(BHT_GXX)*trk/3.0d0
    bht_point_vals(BHT_KXY) = bht_point_vals(BHT_KXY) + bht_point_vals(BHT_GXY)*trk/3.0d0
    bht_point_vals(BHT_KXZ) = bht_point_vals(BHT_KXZ) + bht_point_vals(BHT_GXZ)*trk/3.0d0
    bht_point_vals(BHT_KYY) = bht_point_vals(BHT_KYY) + bht_point_vals(BHT_GYY)*trk/3.0d0
    bht_point_vals(BHT_KYZ) = bht_point_vals(BHT_KYZ) + bht_point_vals(BHT_GYZ)*trk/3.0d0
    bht_point_vals(BHT_KZZ) = bht_point_vals(BHT_KZZ) + bht_point_vals(BHT_GZZ)*trk/3.0d0
end subroutine COCAL_ID_bht_add_trace

! COCAL_ID-only automatic center estimate. Samples exterior COCAL data along
! several directions and extrapolates to r=0 instead of reading tuned values.
subroutine COCAL_ID_bht_compute_smooth_center()

  use COCAL_ID_data_bht
  implicit none
    integer, parameter :: ndir = 14
    real(8) :: dirs(3,ndir), norm
    real(8) :: xfit(bht_fill_points), w(bht_fill_points)
    real(8) :: qs(bht_nvars,bht_fill_points), qdir(bht_nvars)
    real(8) :: sample_q(bht_nvars), sample_r
    real(8) :: avec(6), trk
    integer :: idir, ipt, iv, nvalid
    logical :: sample_ok, dir_ok, split_ok
    real(8), external :: COCAL_ID_lagint_4th_apply

    dirs(:,1)  = (/ 1.0d0, 0.0d0, 0.0d0/)
    dirs(:,2)  = (/-1.0d0, 0.0d0, 0.0d0/)
    dirs(:,3)  = (/ 0.0d0, 1.0d0, 0.0d0/)
    dirs(:,4)  = (/ 0.0d0,-1.0d0, 0.0d0/)
    dirs(:,5)  = (/ 0.0d0, 0.0d0, 1.0d0/)
    dirs(:,6)  = (/ 0.0d0, 0.0d0,-1.0d0/)
    dirs(:,7)  = (/ 1.0d0, 1.0d0, 1.0d0/)
    dirs(:,8)  = (/ 1.0d0, 1.0d0,-1.0d0/)
    dirs(:,9)  = (/ 1.0d0,-1.0d0, 1.0d0/)
    dirs(:,10) = (/ 1.0d0,-1.0d0,-1.0d0/)
    dirs(:,11) = (/-1.0d0, 1.0d0, 1.0d0/)
    dirs(:,12) = (/-1.0d0, 1.0d0,-1.0d0/)
    dirs(:,13) = (/-1.0d0,-1.0d0, 1.0d0/)
    dirs(:,14) = (/-1.0d0,-1.0d0,-1.0d0/)

    smooth_center = 0.0d0
    nvalid = 0
    do idir = 1, ndir
       norm = dsqrt(dirs(1,idir)**2 + dirs(2,idir)**2 + dirs(3,idir)**2)
       dirs(:,idir) = dirs(:,idir)/norm
       dir_ok = .true.
       do ipt = 1, bht_fill_points
          sample_r = r_ah + extra_points_dr*(dble(ipt)-1.0d0+1.0d-8)
          xfit(ipt) = sample_r*sample_r
          call COCAL_ID_bht_eval_raw(xx_bh + radi*sample_r*dirs(1,idir), &
             & yy_bh + radi*sample_r*dirs(2,idir), zz_bh + radi*sample_r*dirs(3,idir), &
             & .true., sample_ok, sample_q)
          if (.not. sample_ok) then
             dir_ok = .false.
             exit
          end if
          call COCAL_ID_bht_split_curvature(sample_q, avec, trk, split_ok)
          if (.not. split_ok) then
             dir_ok = .false.
             exit
          end if
          call COCAL_ID_bht_store_tracefree(sample_q, avec)
          qs(:,ipt) = sample_q
       end do
       if (.not. dir_ok) cycle
       call COCAL_ID_lagint_4th_weights(xfit, 0.0d0, w)
       do iv = 1, bht_nvars
          qdir(iv) = COCAL_ID_lagint_4th_apply(w, qs(iv,:))
       end do
       smooth_center = smooth_center + qdir
       nvalid = nvalid + 1
    end do

    if (nvalid > 0) then
       smooth_center = smooth_center/dble(nvalid)
       smooth_center_valid = .true.
    else
       smooth_center_valid = .false.
    end if

    call COCAL_ID_bht_zero_hydro(smooth_center)
    smooth_center_attempted = .true.
end subroutine COCAL_ID_bht_compute_smooth_center

! COCAL_ID-only BHT point fill. It uses raw PCOCAL-style exterior interpolation
! outside the fill radius and the automatic/legacy smooth fill inside it.
subroutine COCAL_ID_bht_fill_point(xcac, ycac, zcac, use_IDpar, ok, bht_point_vals)

  use COCAL_ID_data_bht
  implicit none
    real(8), intent(in) :: xcac, ycac, zcac
    logical, intent(in) :: use_IDpar
    logical, intent(out) :: ok
    real(8), intent(out) :: bht_point_vals(bht_nvars)
    real(8) :: rr, rr_coc, ux, uy, uz, xrel, yrel, zrel
    real(8) :: xfit(bht_fill_points), w(bht_fill_points)
    real(8) :: qs(bht_nvars,bht_fill_points)
    real(8) :: sample_q(bht_nvars)
    real(8) :: avec(6)
    real(8) :: sample_r, ss, trk, trk_first, trk_inner
    real(8) :: h00, h01, h11, dqdr
    integer :: ipt, iv
    logical :: sample_ok, split_ok, project_ok
    real(8), external :: COCAL_ID_lagint_4th_apply

    bht_point_vals = 0.0d0
    trk_first = 0.0d0
    ok = .false.

    xrel = xcac - xx_bh
    yrel = ycac - yy_bh
    zrel = zcac - zz_bh
    rr = dsqrt(xrel*xrel + yrel*yrel + zrel*zrel)

    rr_coc = rr/radi

    if (rr_coc >= r_ah) then
       call COCAL_ID_bht_eval_raw(xcac, ycac, zcac, .false., ok, bht_point_vals)
       call COCAL_ID_bht_validate_state(ok, bht_point_vals)
       return
    end if

    if (rr > 1.0d-14) then
       ux = xrel/rr
       uy = yrel/rr
       uz = zrel/rr
    else
       ux = 1.0d0
       uy = 0.0d0
       uz = 0.0d0
    end if

    if (use_IDpar) then
       xfit(1) = 0.0d0
       qs(:,1) = 0.0d0
       qs(BHT_ALP,1) = lap_cent
       qs(BHT_BETAX,1) = betax_cent
       qs(BHT_BETAY,1) = betay_cent
       qs(BHT_BETAZ,1) = betaz_cent
       qs(BHT_GXX,1) = gxx_cent
       qs(BHT_GXY,1) = gxy_cent
       qs(BHT_GXZ,1) = gxz_cent
       qs(BHT_GYY,1) = gyy_cent
       qs(BHT_GYZ,1) = gyz_cent
       qs(BHT_GZZ,1) = gzz_cent
       qs(BHT_KXX,1) = kxx_cent
       qs(BHT_KXY,1) = kxy_cent
       qs(BHT_KXZ,1) = kxz_cent
       qs(BHT_KYY,1) = kyy_cent
       qs(BHT_KYZ,1) = kyz_cent
       qs(BHT_KZZ,1) = kzz_cent
       qs(BHT_RHO,1) = rho_cent
       qs(BHT_VX,1) = vx_cent
       qs(BHT_VY,1) = vy_cent
       qs(BHT_VZ,1) = vz_cent
       do ipt = 2, bht_fill_points
          sample_r = r_ah + extra_points_dr*(dble(ipt)-2.0d0+1.0d-8)
          xfit(ipt) = sample_r
          call COCAL_ID_bht_eval_raw(xx_bh + radi*sample_r*ux, yy_bh + radi*sample_r*uy, &
             & zz_bh + radi*sample_r*uz, &
             & .true., sample_ok, sample_q)
          if (.not. sample_ok) return
          if (ipt == 2) then
             call COCAL_ID_bht_split_curvature(sample_q, avec, trk_first, split_ok)
          else
             call COCAL_ID_bht_split_curvature(sample_q, avec, trk, split_ok)
          end if
          if (.not. split_ok) return
          call COCAL_ID_bht_store_tracefree(sample_q, avec)
          qs(:,ipt) = sample_q
       end do
       call COCAL_ID_lagint_4th_weights(xfit, rr_coc, w)
       do iv = 1, bht_nvars
          bht_point_vals(iv) = COCAL_ID_lagint_4th_apply(w, qs(iv,:))
       end do
    else
       if (.not. smooth_center_attempted) call COCAL_ID_bht_compute_smooth_center()
       if (.not. smooth_center_valid) return
       ! Hermite matching only needs the exterior value and radial slope.
       do ipt = 1, 2
          sample_r = r_ah + extra_points_dr*(dble(ipt)-1.0d0+1.0d-8)
          call COCAL_ID_bht_eval_raw(xx_bh + radi*sample_r*ux, yy_bh + radi*sample_r*uy, &
             & zz_bh + radi*sample_r*uz, &
             & .true., sample_ok, sample_q)
          if (.not. sample_ok) return
          if (ipt == 1) then
             call COCAL_ID_bht_split_curvature(sample_q, avec, trk_first, split_ok)
          else
             call COCAL_ID_bht_split_curvature(sample_q, avec, trk, split_ok)
          end if
          if (.not. split_ok) return
          call COCAL_ID_bht_store_tracefree(sample_q, avec)
          qs(:,ipt) = sample_q
       end do

       ss = max(0.0d0, min(1.0d0, rr_coc/r_ah))
       h00 = 2.0d0*ss**3 - 3.0d0*ss**2 + 1.0d0
       h01 = -2.0d0*ss**3 + 3.0d0*ss**2
       h11 = ss**3 - ss**2
       do iv = 1, bht_nvars
          dqdr = (qs(iv,2) - qs(iv,1))/extra_points_dr
          bht_point_vals(iv) = h00*smooth_center(iv) + h01*qs(iv,1) + h11*r_ah*dqdr
       end do
    end if

    call COCAL_ID_bht_project_tracefree(bht_point_vals, project_ok)
    if (.not. project_ok) return
    if (trk_xi > 0.0d0 .and. r_ah > 0.0d0) then
       trk_inner = trk_first*dexp(-((max(0.0d0, r_ah - rr_coc))/(trk_xi*r_ah))**trk_npow)
    else
       trk_inner = 0.0d0
    end if
    call COCAL_ID_bht_add_trace(bht_point_vals, trk_inner)

    call COCAL_ID_bht_zero_hydro(bht_point_vals)
    ok = .true.
    call COCAL_ID_bht_validate_state(ok, bht_point_vals)
end subroutine COCAL_ID_bht_fill_point
