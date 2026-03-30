module interface_IO_input_gradvep_export 
  implicit none
  interface
    subroutine IO_input_gradvep_export(filenm,coc2cac_readformat, vepxf, vepyf, vepzf)
      real(8), pointer :: vepxf(:,:,:), vepyf(:,:,:), vepzf(:,:,:)
      character(len=*) :: filenm,coc2cac_readformat
    end subroutine IO_input_gradvep_export
  end interface
end module interface_IO_input_gradvep_export
