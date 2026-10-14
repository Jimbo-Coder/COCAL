module COCAL_ID_interface_IO_input_CF_flir_export
implicit none
  interface 
  
subroutine COCAL_ID_IO_input_CF_flir_export(filenm,coc2cac_readformat, emd,vep,ome,ber,radi,xcm,zcm)
      real(8), pointer :: emd(:,:,:), vep(:,:,:)
      real(8) ::  ome,ber,radi,xcm,zcm
      character(len=*) :: filenm,coc2cac_readformat
  
end subroutine COCAL_ID_IO_input_CF_flir_export
  end interface
end module COCAL_ID_interface_IO_input_CF_flir_export
