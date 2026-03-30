subroutine IO_input_matter_BHT_export(filenm,coc2cac_readformat, emdg, omeg, ome, ber, radi)
  use phys_constant, only : long
  implicit none
  integer :: ir, it, ip, nrtmp, nttmp, nptmp
  real(8), pointer :: emdg(:,:,:), omeg(:,:,:)
  real(8) :: ome, ber, radi, emdc
  character(len=*) :: filenm, coc2cac_readformat
!
  write(6,*) "Reading emdg, omeg..."   
! --- Matter
  open(12,file=trim(filenm),status='old')
  read(12,'(5i5)')  nrtmp, nttmp, nptmp
  if (coc2cac_readformat == "23.15") then
      do ip = 0, nptmp
         do it = 0, nttmp
            do ir = 0, nrtmp
            read(12,'(1p,6e23.15)') emdg(ir,it,ip),  omeg(ir,it,ip)
            end do
         end do
      end do
      read(12,'(1p,6e23.15)') ome, ber, radi
      !  read(12,'(1p,6e23.15)') emdc
   else if (coc2cac_readformat == "20.12") then
      do ip = 0, nptmp
         do it = 0, nttmp
            do ir = 0, nrtmp
            read(12,'(1p,6e20.12)') emdg(ir,it,ip),  omeg(ir,it,ip)
            end do
         end do
      end do
      read(12,'(1p,6e20.12)') ome, ber, radi
      !  read(12,'(1p,6e20.12)') emdc
  end if   
  close(12)
!
end subroutine IO_input_matter_BHT_export
