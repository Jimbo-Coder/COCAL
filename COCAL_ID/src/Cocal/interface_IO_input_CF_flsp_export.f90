module interface_IO_input_CF_flsp_export
  implicit none
  interface 
    subroutine IO_input_CF_flsp_export(filenm,coc2cac_readformat,emd,vep,wxspf,wyspf,wzspf,ome,ber,radi,confpow)
      real(8), pointer :: emd(:,:,:), vep(:,:,:), wxspf(:,:,:), wyspf(:,:,:), wzspf(:,:,:)
      real(8) ::  ome,ber,radi, confpow, omespx, omespy, omespz
      character(len=*) :: filenm,coc2cac_readformat
    end subroutine IO_input_CF_flsp_export
  end interface
end module interface_IO_input_CF_flsp_export
