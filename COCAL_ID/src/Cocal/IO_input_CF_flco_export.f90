subroutine COCAL_ID_IO_input_CF_flco_export(filenm,coc2cac_readformat,emd,ome,ber,radi,xcm)

use COCAL_ID_phys_constant, only : long, nnrg, nntg, nnpg
implicit none
  integer :: ir, it, ip, nrtmp, nttmp, nptmp
  real(8), pointer :: emd(:,:,:)
  real(8) :: ome, ber, radi, xcm
  character(len=*) :: filenm,coc2cac_readformat
!
! --- Matter
    open(12,file=trim(filenm),status='old')
    read(12,'(5i5)')  nrtmp, nttmp, nptmp
    if (coc2cac_readformat == "20.12") then
      do ip = 0, nptmp
         do it = 0, nttmp
            do ir = 0, nrtmp
               read(12,'(1p,6e20.12)') emd(ir,it,ip)
            end do
         end do
      end do
      ! Legacy records have three fields; formatted input pads xcm with blanks (zero).
      read(12,'(1p,6e20.12)') ome, ber, radi, xcm
   else if (coc2cac_readformat == "23.15") then
      do ip = 0, nptmp
         do it = 0, nttmp
            do ir = 0, nrtmp
               read(12,'(1p,6e23.15)') emd(ir,it,ip)
            end do
         end do
      end do
      read(12,'(1p,6e23.15)') ome, ber, radi, xcm
   end if
    close(12)
!
end subroutine COCAL_ID_IO_input_CF_flco_export
