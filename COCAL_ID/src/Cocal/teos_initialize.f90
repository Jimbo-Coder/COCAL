subroutine COCAL_ID_teos_initialize
!

use COCAL_ID_phys_constant        !g,c,solmas,nnteos
use COCAL_ID_def_teos_parameter   !prei,enei,kappa_crust,gamma_crust
implicit none
!
  real(8) :: facrho
  real(8) :: ene, pre, hh, rho
  integer :: ii, iphase, idum, i0
  real(8), external :: COCAL_ID_lagint_4th
  real(8)           :: x4(4), f4(4)
!
  facrho = (g/c**2)**3*solmas**2

  open(850,file='teos_parameter.dat',status='old')
  read(850,'(i5)') idum
  nphase=idum-1
  do ii=0,nphase
    read(850,'(1p,4e23.15)') ene, pre, hh, rho
    if (ii==0 .or. ii==nphase) write(6,'(i5,1p,4e23.15)') ii, ene, pre, hh, rho
    enei(ii) = ene
    rhoi(ii) = rho
    prei(ii) = pre
    hi(ii)   = hh
    qi(ii)   = prei(ii)/rhoi(ii)
  end do
  close(850)
  write(6,*) "i, rhoi, prei, hi, qi"
  write(6,'(i5,1p,4e23.15)') 0, rhoi(0), prei(0), hi(0), qi(0)
  write(6,'(i5,1p,4e23.15)') nphase, rhoi(nphase), prei(nphase), hi(nphase), qi(nphase)
!
  open(850,file='peos_parameter.dat',status='old')
  read(850,'(8x,1i5,es13.5)') idum, rhoini_cgs
  read(850,'(2es13.5)') kappa_crust, gamma_crust
  close(850)
  write(6,'(a21,1p,e23.15)') "******* rhoini_cgs = ", rhoini_cgs
  rhoini_gcm1 = facrho*rhoini_cgs

  call COCAL_ID_teos_lookup(rhoini_gcm1,rhoi,iphase)
  i0 = min0(max0(iphase-2,0),nphase-3)
  x4(1:4) = rhoi(i0:i0+3)
  f4(1:4) = qi(i0:i0+3)
  emdini_gcm1 = COCAL_ID_lagint_4th(x4,f4,rhoini_gcm1)
  write(6,'(a26,1p,2e23.15)') "(rhoini_gcm1,emdini_gcm1)=", rhoini_gcm1, emdini_gcm1
!
end subroutine COCAL_ID_teos_initialize

subroutine COCAL_ID_teos_initialize_cactus(dir_path)
!

use COCAL_ID_phys_constant        !g,c,solmas,nnteos
use COCAL_ID_def_teos_parameter   !prei,enei,kappa_crust,gamma_crust
implicit none
  character*400, intent(in) :: dir_path
  real(8) :: facrho
  real(8) :: ene, pre, hh, rho
  integer :: ii, iphase, idum, i0
  real(8), external :: COCAL_ID_lagint_4th
  real(8)           :: x4(4), f4(4)
!
  facrho = (g/c**2)**3*solmas**2

  open(850,file=trim(dir_path)//'/'//'teos_parameter.dat',status='old')
  read(850,'(i5)') idum
  nphase=idum-1
  do ii=0,nphase
    read(850,'(1p,4e23.15)') ene, pre, hh, rho
    if (ii==0 .or. ii==nphase) write(6,'(i5,1p,4e23.15)') ii, ene, pre, hh, rho
    enei(ii) = ene
    rhoi(ii) = rho
    prei(ii) = pre
    hi(ii)   = hh
    qi(ii)   = prei(ii)/rhoi(ii)
  end do
  close(850)
  write(6,*) "i, rhoi, prei, hi, qi"
  write(6,'(i5,1p,4e23.15)') 0, rhoi(0), prei(0), hi(0), qi(0)
  write(6,'(i5,1p,4e23.15)') nphase, rhoi(nphase), prei(nphase), hi(nphase), qi(nphase)
!
  open(850,file=trim(dir_path)//'/'//'peos_parameter.dat',status='old')
  read(850,'(8x,1i5,es13.5)') idum, rhoini_cgs
  read(850,'(2es13.5)') kappa_crust, gamma_crust
  close(850)
  write(6,'(a21,1p,e23.15)') "******* rhoini_cgs = ", rhoini_cgs
  rhoini_gcm1 = facrho*rhoini_cgs

  call COCAL_ID_teos_lookup(rhoini_gcm1,rhoi,iphase)
  i0 = min0(max0(iphase-2,0),nphase-3)
  x4(1:4) = rhoi(i0:i0+3)
  f4(1:4) = qi(i0:i0+3)
  emdini_gcm1 = COCAL_ID_lagint_4th(x4,f4,rhoini_gcm1)
  write(6,'(a26,1p,2e23.15)') "(rhoini_gcm1,emdini_gcm1)=", rhoini_gcm1, emdini_gcm1
!
end subroutine COCAL_ID_teos_initialize_cactus

subroutine COCAL_ID_teos_initialize_mpt_cactus(impt, dir_path)
!

use COCAL_ID_phys_constant        !g,c,solmas,nnteos
use COCAL_ID_def_teos_parameter   !prei,enei,kappa_crust,gamma_crust
implicit none
  integer, intent(in)       :: impt
  character*400, intent(in) :: dir_path
  character(len=1)          :: np(2) = (/'1', '2'/)
  real(8)                   :: facrho
  real(8)                   :: ene, pre, hh, rho
  integer                   :: ii, iphase, idum, i0
  real(8), external         :: COCAL_ID_lagint_4th
  real(8)                   :: x4(4), f4(4)
!
  facrho = (g/c**2)**3*solmas**2

  open(850,file=trim(dir_path)//'/'//'teos_parameter_mpt'//np(impt)//'.dat',status='old')
  read(850,'(i5)') idum
  nphase=idum-1
  do ii=0,nphase
    read(850,'(1p,4e23.15)') ene, pre, hh, rho
    if (ii==0 .or. ii==nphase) write(6,'(i5,1p,4e23.15)') ii, ene, pre, hh, rho
    enei(ii) = ene
    rhoi(ii) = rho
    prei(ii) = pre
    hi(ii)   = hh
    qi(ii)   = prei(ii)/rhoi(ii)
  end do
  close(850)
  write(6,*) "i, rhoi, prei, hi, qi"
  write(6,'(i5,1p,4e23.15)') 0, rhoi(0), prei(0), hi(0), qi(0)
  write(6,'(i5,1p,4e23.15)') nphase, rhoi(nphase), prei(nphase), hi(nphase), qi(nphase)
!
  open(850,file=trim(dir_path)//'/'//'peos_parameter_mpt'//np(impt)//'.dat',status='old')
  read(850,'(8x,1i5,es13.5)') idum, rhoini_cgs
  read(850,'(2es13.5)') kappa_crust, gamma_crust
  close(850)
  write(6,'(a21,1p,e23.15)') "******* rhoini_cgs = ", rhoini_cgs
  rhoini_gcm1 = facrho*rhoini_cgs

  call COCAL_ID_teos_lookup(rhoini_gcm1,rhoi,iphase)
  i0 = min0(max0(iphase-2,0),nphase-3)
  x4(1:4) = rhoi(i0:i0+3)
  f4(1:4) = qi(i0:i0+3)
  emdini_gcm1 = COCAL_ID_lagint_4th(x4,f4,rhoini_gcm1)
  write(6,'(a26,1p,2e23.15)') "(rhoini_gcm1,emdini_gcm1)=", rhoini_gcm1, emdini_gcm1
!
end subroutine COCAL_ID_teos_initialize_mpt_cactus
