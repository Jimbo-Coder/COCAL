subroutine IO_input_gradvep_export(filenm,coc2cac_readformat, vepxf, vepyf, vepzf)                                                                                      
  use phys_constant, only : long
  implicit none
  integer :: ir, it, ip, nrtmp, nttmp, nptmp
  real(8), pointer :: vepxf(:,:,:), vepyf(:,:,:), vepzf(:,:,:)
  character(len=*) :: filenm,coc2cac_readformat
!
  write(6,*) "Reading vepxf, vepyf, vepzf..."
  open(13,file=trim(filenm),status='old')
  read(13,'(5i5)') nrtmp, nttmp, nptmp
  if (coc2cac_readformat == "20.12") then
   do ip = 0, nptmp
      do it = 0, nttmp
         do ir = 0, nrtmp
            read(13,'(1p,3e20.12)')  vepxf(ir,it,ip), vepyf(ir,it,ip), vepzf(ir,it,ip)
         end do
      end do
   end do
  else if (coc2cac_readformat == "23.15") then
   do ip = 0, nptmp
      do it = 0, nttmp
         do ir = 0, nrtmp
            read(13,'(1p,3e23.15)')  vepxf(ir,it,ip), vepyf(ir,it,ip), vepzf(ir,it,ip)
         end do
      end do
   end do
  end if
  close(13)
!
end subroutine IO_input_gradvep_export
