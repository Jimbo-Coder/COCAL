! The Carpet grid: the COCAL values at every point of the Cactus grid.
! Included by COCAL_ID_master.F90, which is the same in COCALX (CarpetX);
! only this file, the driver's, differs.

subroutine COCAL_ID_set_data(CCTK_ARGUMENTS)
  use COCAL_ID_state
  use COCAL_ID_eval
  use COCAL_ID_rns, only : rns_em, RNS_EM_POTENTIAL, RNS_EM_FARADAY
  implicit none
  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_PARAMETERS
  DECLARE_CCTK_FUNCTIONS
  ! IllinoisGRMHD stores Ax(i,j+1/2,k+1/2), Ay(i+1/2,j,k+1/2),
  ! Az(i+1/2,j+1/2,k) and Aphi(i+1/2,j+1/2,k+1/2): half-step offsets.
  real(8), parameter :: stagger(3,4) = reshape([0d0,1d0,1d0, 1d0,0d0,1d0, 1d0,1d0,0d0, 1d0,1d0,1d0], [3,4])
  real(8) :: s(NSTATE), so(NSTATE), offset(3)
  logical :: set_lapse, set_shift, set_hydro, set_bvec, set_avec, set_aphi, staggered, failed, ok
  integer :: i, j, k, n

  if (cocal_kind == 0) call CCTK_WARN(CCTK_WARN_ABORT, "COCAL_ID_set_data called before COCAL_ID_read_data.")
  set_lapse = CCTK_EQUALS(initial_lapse, "Cocal")
  set_shift = CCTK_EQUALS(initial_shift, "Cocal")
  set_hydro = CCTK_EQUALS(initial_hydro, "Cocal")
  set_bvec = CCTK_EQUALS(initial_Bvec, "Cocal")
  set_avec = CCTK_EQUALS(initial_Avec, "Cocal")
  set_aphi = CCTK_EQUALS(initial_Aphi, "Cocal")
  if (set_avec .neqv. set_aphi) &
    call CCTK_WARN(CCTK_WARN_ABORT, "COCAL_ID: initial_Avec and initial_Aphi must both be Cocal or neither.")
  if (set_avec .and. rns_em /= RNS_EM_POTENTIAL) &
    call CCTK_WARN(CCTK_WARN_ABORT, "COCAL_ID: Cocal Avec/Aphi need MRNS_WL data with COCAL_ID_read_magnetic_potential = yes.")
  if (set_bvec .and. rns_em /= RNS_EM_FARADAY) &
    call CCTK_WARN(CCTK_WARN_ABORT, "COCAL_ID: Cocal Bvec needs MRNS_WL data with COCAL_ID_read_magnetic_potential = no.")
  staggered = set_avec .and. COCAL_ID_offset_AvecAphi == 1

  failed = .false.
  !$OMP parallel do collapse(3) default(shared) schedule(static) &
  !$OMP& private(i,j,k,n,s,so,offset,ok) reduction(.or.:failed)
  do k = 1, cctk_lsh(3)
    do j = 1, cctk_lsh(2)
      do i = 1, cctk_lsh(1)
        call COCAL_ID_eval_point(x(i,j,k), y(i,j,k), z(i,j,k), s, ok)
        if (.not. ok) then
          failed = .true.
          cycle
        end if
        gxx(i,j,k) = s(S_G);     gxy(i,j,k) = s(S_G+1);   gxz(i,j,k) = s(S_G+2)
        gyy(i,j,k) = s(S_G+3);   gyz(i,j,k) = s(S_G+4);   gzz(i,j,k) = s(S_G+5)
        kxx(i,j,k) = s(S_K);     kxy(i,j,k) = s(S_K+1);   kxz(i,j,k) = s(S_K+2)
        kyy(i,j,k) = s(S_K+3);   kyz(i,j,k) = s(S_K+4);   kzz(i,j,k) = s(S_K+5)
        if (set_lapse) alp(i,j,k) = s(S_ALP)
        if (set_shift) then
          betax(i,j,k) = s(S_BETA); betay(i,j,k) = s(S_BETA+1); betaz(i,j,k) = s(S_BETA+2)
        end if
        if (set_hydro) then
          rho(i,j,k) = s(S_RHO)
          press(i,j,k) = s(S_PRESS)
          eps(i,j,k) = s(S_EPS)
          vel(i,j,k,:) = s(S_VEL:S_VEL+2)
          w_lorentz(i,j,k) = 1.0d0/sqrt(1.0d0 - COCAL_ID_speed2(s))
        end if
        if (set_bvec) Bvec(i,j,k,:) = s(S_BVEC:S_BVEC+2)
        if (set_avec) then
          Avec(i,j,k,:) = s(S_AVEC:S_AVEC+2)
          Aphi(i,j,k) = s(S_APHI)
          if (staggered) then
            do n = 1, 4
              offset = 0.5d0*stagger(:,n)*[CCTK_DELTA_SPACE(1), CCTK_DELTA_SPACE(2), CCTK_DELTA_SPACE(3)]
              call COCAL_ID_eval_point(x(i,j,k) + offset(1), y(i,j,k) + offset(2), z(i,j,k) + offset(3), so, ok)
              failed = failed .or. .not. ok
              if (n <= 3) then
                Avec(i,j,k,n) = so(S_AVEC + n - 1)
              else
                Aphi(i,j,k) = so(S_APHI)
              end if
            end do
          end if
        end if
      end do
    end do
  end do
  !$OMP end parallel do
  ! Invalid: non-finite, nonpositive psi or alpha, or v >= 1.
  if (failed) call CCTK_WARN(CCTK_WARN_ABORT, "COCAL_ID: a Cactus point is outside the COCAL data or invalid there.")
end subroutine COCAL_ID_set_data
