subroutine IO_input_CF_surf_export(filenm,coc2cac_readformat,rs)
  use phys_constant, only : long, nnrg, nntg, nnpg
  implicit none
  integer :: ir, it, ip, nrtmp, nttmp, nptmp
  real(8), pointer ::  rs(:,:)
  character(len=*) :: filenm,coc2cac_readformat
!
! --- Star surface
  open(15,file=trim(filenm),status='old')
  if (coc2cac_readformat =="20.12") then
   read(15,'(5i5)')   nttmp,  nptmp
   do ip = 0, nptmp
      do it = 0, nttmp
         read(15,'(1p,6e20.12)')   rs(it,ip)
      end do
   end do
  else if (coc2cac_readformat == "23.15") then
      read(15,'(5i5)')   nttmp,  nptmp
      do ip = 0, nptmp
         do it = 0, nttmp
            read(15,'(1p,6e23.15)')   rs(it,ip)
         end do
      end do
  end if
  close(15)
!
end subroutine IO_input_CF_surf_export
