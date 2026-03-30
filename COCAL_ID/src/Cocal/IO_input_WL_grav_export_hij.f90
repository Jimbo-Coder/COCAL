subroutine IO_input_WL_grav_export_hij(filenm,coc2cac_readformat,hxxd,hxyd,hxzd,hyyd,hyzd,hzzd)
  use phys_constant, only : long
  implicit none
  integer :: irg, itg, ipg, nrtmp, nttmp, nptmp
  real(8), pointer :: hxxd(:,:,:), hxyd(:,:,:), hxzd(:,:,:), &
      &               hyyd(:,:,:), hyzd(:,:,:), hzzd(:,:,:)
  character(len=*) :: filenm, coc2cac_readformat
!
  write(6,*) "Reading hij..."
! --- Metric potentials.
  open(13,file=trim(filenm),status='old')
  read(13,'(5i5)')  nrtmp, nttmp, nptmp
  if (coc2cac_readformat == "23.15") then
    do ipg = 0, nptmp
      do itg = 0, nttmp
         do irg = 0, nrtmp
            read(13,'(1p,6e23.15)')  hxxd(irg,itg,ipg), &
            &                        hxyd(irg,itg,ipg), &
            &                        hxzd(irg,itg,ipg), &
            &                        hyyd(irg,itg,ipg), &
            &                        hyzd(irg,itg,ipg), &
            &                        hzzd(irg,itg,ipg)
         end do
      end do
    end do
   else if (coc2cac_readformat == "20.12") then
    do ipg = 0, nptmp
      do itg = 0, nttmp
         do irg = 0, nrtmp
            read(13,'(1p,6e20.12)')  hxxd(irg,itg,ipg), &
            &                        hxyd(irg,itg,ipg), &
            &                        hxzd(irg,itg,ipg), &
            &                        hyyd(irg,itg,ipg), &
            &                        hyzd(irg,itg,ipg), &
            &                        hzzd(irg,itg,ipg)
         end do
      end do
     end do
   end if
  close(13)
!
end subroutine IO_input_WL_grav_export_hij
