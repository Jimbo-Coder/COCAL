module COCAL_ID_interface_IO_input_CF_flsp_export
implicit none
  interface 
  
subroutine COCAL_ID_IO_input_CF_flsp_export(filenm,coc2cac_readformat,emd,vep,wxspf,wyspf,wzspf,ome,ber,radi,confpow,xcm)
      real(8), pointer :: emd(:,:,:), vep(:,:,:), wxspf(:,:,:), wyspf(:,:,:), wzspf(:,:,:)
      real(8) ::  ome,ber,radi, confpow, xcm
      character(len=*) :: filenm,coc2cac_readformat
  
end subroutine COCAL_ID_IO_input_CF_flsp_export
  end interface
end module COCAL_ID_interface_IO_input_CF_flsp_export
