subroutine COCAL_ID_coordinate_patch_kit_grav_grid_coc2cac_mpt(igrid,impt,dir_path,readformat,compact)

use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
use COCAL_ID_phys_constant, only : pi
use COCAL_ID_grid_parameter
use COCAL_ID_coordinate_grav_r
use COCAL_ID_coordinate_grav_phi
use COCAL_ID_coordinate_grav_theta
use COCAL_ID_coordinate_grav_extended
use COCAL_ID_trigonometry_grav_theta
use COCAL_ID_trigonometry_grav_phi
implicit none
  integer, intent(in) :: igrid, impt
  character(len=*), intent(in) :: dir_path, readformat
  logical, intent(in) :: compact
  character(len=1) :: patch
  character(len=:), allocatable :: filename
  character(len=20) :: fmt
  character(len=256) :: message
  integer :: unit, ios, nr, nt, np, ir, it, ip, nr_finite
  logical :: exported
  real(long), parameter :: grid_tol = 1.0d-10

  write(patch,'(i1)') impt
  filename = trim(dir_path)//'/bnsgrids_3D_mpt'//patch//'.las'
  inquire(file=filename,exist=exported)
  nr_finite = nrg
  if (compact .and. (.not. exported .or. impt /= 3)) &
    call invalid_grid('Compact mode requires an exported outer (mpt3) grid.')

  if (exported) then
    if (readformat /= '20.12' .and. readformat /= '23.15') &
      call invalid_grid('Unknown real-number format.')
    fmt = '(1p,e'//trim(readformat)//')'
    open(newunit=unit,file=filename,status='old',action='read',iostat=ios)
    if (ios /= 0) call invalid_grid('Cannot open exported grid.')
    read(unit,'(5i5)',iostat=ios) nr, nt, np
    if (ios /= 0) call invalid_grid('Cannot read grid dimensions.')
    if (nr /= nrg .or. nt /= ntg .or. np /= npg) then
      write(message,'(a,3(1x,i0),a,3(1x,i0))') &
        'Grid dimensions differ from rnspar: exported nr/nt/np =', nr, nt, np, &
        '; expected =', nrg, ntg, npg
      call invalid_grid(trim(message))
    end if
    if (nr < 3 .or. nr > ubound(rg,1) .or. nt < 4 .or. nt > ubound(thg,1) .or. &
        np < 4 .or. np > ubound(phig,1)) call invalid_grid('Unsupported grid dimensions.')
    if (mod(nt,2) /= 0 .or. mod(np,4) /= 0) &
      call invalid_grid('Angular grids require even ntg and npg divisible by four.')
    do ir = 0, nr
      read(unit,fmt,iostat=ios) rg(ir)
      if (ios /= 0) call invalid_grid('Cannot read radial coordinates.')
    end do
    do it = 0, nt
      read(unit,fmt,iostat=ios) thg(it)
      if (ios /= 0) call invalid_grid('Cannot read theta coordinates.')
    end do
    do ip = 0, np
      read(unit,fmt,iostat=ios) phig(ip)
      if (ios /= 0) call invalid_grid('Cannot read phi coordinates.')
    end do
    close(unit)
    if (.not. all(ieee_is_finite(rg(0:nr))) .or. &
        .not. all(ieee_is_finite(thg(0:nt))) .or. &
        .not. all(ieee_is_finite(phig(0:np)))) call invalid_grid('Non-finite coordinates.')
    if (rg(0) < 0.0d0 .or. any(rg(1:nr) <= rg(0:nr-1))) &
      call invalid_grid('Radial coordinates must be nonnegative and increasing.')

    if (compact) then
      ! The last exported coordinate denotes infinity, not a usable data sample.
      ! Keep nrg for reading the file's full field arrays; limit all stencils below.
      nr_finite = nrg - 1
      if (nr_finite < 3 .or. rgin <= 0.0d0 .or. rgout <= rgin) &
        call invalid_grid('Invalid compact radial domain.')
      if (abs(rg(nrg)-rgout) > grid_tol*abs(rgout) .or. &
          abs(rg(0)-rgin) > grid_tol*abs(rgin) .or. any(rg(0:nr_finite) >= rgout)) &
        call invalid_grid('Compact endpoints do not match rgin/rgout.')
      rg(0:nr_finite) = rgin*(rgout-rgin)/(rgout-rg(0:nr_finite))
      if (.not. all(ieee_is_finite(rg(0:nr_finite)))) &
        call invalid_grid('Non-finite decompactified radii.')
    end if

    ! Rebuild auxiliary coordinates from the exported points, not a grid generator.
    rginv(0) = 0.0d0
    if (rg(0) > 0.0d0) rginv(0) = 1.0d0/rg(0)
    do ir = 1, nr_finite
      drg(ir) = rg(ir)-rg(ir-1)
      drginv(ir) = 1.0d0/drg(ir)
      hrg(ir) = 0.5d0*(rg(ir)+rg(ir-1))
      rginv(ir) = 1.0d0/rg(ir)
      hrginv(ir) = 1.0d0/hrg(ir)
    end do
    ! Reflection/quarter-turn index maps require the usual uniform angular grids.
    dthg = pi/dble(ntg)
    dthginv = 1.0d0/dthg
    do it = 0, ntg
      if (abs(thg(it)-dble(it)*dthg) > grid_tol) &
        call invalid_grid('Unsupported nonuniform theta grid.')
    end do
    hthg(1:ntg) = 0.5d0*(thg(1:ntg)+thg(0:ntg-1))
    dphig = 2.0d0*pi/dble(npg)
    dphiginv = 1.0d0/dphig
    do ip = 0, npg
      if (abs(phig(ip)-dble(ip)*dphig) > grid_tol) &
        call invalid_grid('Unsupported nonuniform phi grid.')
    end do
    hphig(1:npg) = 0.5d0*(phig(1:npg)+phig(0:npg-1))
  else
    ! Legacy data without exported coordinates retain their original generator.
    select case (igrid)
    case (1)
      call grid_r
    case (2)
      call grid_r_bhex('eBH')
    case (3)
      call grid_r_bns
    case (4)
      call grid_r_bqs
    case default
      call invalid_grid('Choose legacy grid 1, 2, 3 or 4.')
    end select
    call grid_theta
    call grid_phi
  end if

  if (impt <= 2) then
    if (nrf < 3 .or. nrf > nrg .or. ntf /= ntg .or. npf /= npg) &
      call invalid_grid('Fluid grid must use the gravitational angular grid and at least four radial points.')
    if (.not. ieee_is_finite(r_surf) .or. r_surf <= 0.0d0) &
      call invalid_grid('Invalid r_surf.')
    if (abs(rg(nrf)-r_surf) > grid_tol*max(1.0d0,abs(r_surf))) then
      write(message,'(a,i0,a,es20.12,a,es20.12)') &
        'rg(nrf) differs from r_surf: nrf=', nrf, '; rg(nrf)=', rg(nrf), '; r_surf=', r_surf
      call invalid_grid(trim(message))
    end if
  end if
  call trig_grav_theta
  call trig_grav_phi
  call grid_extended(nr_finite)

contains
  subroutine invalid_grid(message)
    character(len=*), intent(in) :: message
    write(*,'(a)') 'COCAL_ID: '//filename//': '//message
    error stop 'Invalid BNS grid'
  end subroutine invalid_grid
end subroutine COCAL_ID_coordinate_patch_kit_grav_grid_coc2cac_mpt
