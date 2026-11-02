module COCAL_ID_interface_IO_input_matter_BHT_export
implicit none
  interface 
  
subroutine COCAL_ID_IO_input_matter_BHT_export(filenm,coc2cac_readformat,emdg,omeg,ome,ber,radi)
      real(8), pointer :: emdg(:,:,:), omeg(:,:,:)
      real(8) ::  ome,ber,radi
      character(len=*) :: filenm, coc2cac_readformat
  
end subroutine COCAL_ID_IO_input_matter_BHT_export
  end interface
end module COCAL_ID_interface_IO_input_matter_BHT_export