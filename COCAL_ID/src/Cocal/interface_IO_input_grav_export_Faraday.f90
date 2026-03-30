module interface_IO_input_grav_export_Faraday
  implicit none
  interface 
    subroutine IO_input_grav_export_Faraday(filenm,coc2cac_readformat,fxd,fyd,fzd,fxyd,fxzd,fyzd)
      real(8), pointer ::  fxd(:,:,:),  fyd(:,:,:),  fzd(:,:,:), &
          &               fxyd(:,:,:), fxzd(:,:,:), fyzd(:,:,:)
      character(len=*) :: filenm,coc2cac_readformat
    end subroutine IO_input_grav_export_Faraday
  end interface
end module interface_IO_input_grav_export_Faraday
