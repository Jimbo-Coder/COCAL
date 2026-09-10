module COCAL_ID_interface_interpo_lag4th_2Dsurf
implicit none
  interface 
  
subroutine COCAL_ID_interpo_lag4th_2Dsurf(val,fnc,tv,pv,ntg,npg, &
                                      thgex,phigex,itgex_th,ipgex_phi,ipgex_th)
      integer, intent(in) :: ntg, npg
      real(8), intent(out) :: val
      real(8), intent(in)  :: tv, pv
      real(8), pointer :: fnc(:,:)
      real(8), intent(in) :: thgex(-2:), phigex(-2:)
      integer, intent(in) :: itgex_th(-2:), ipgex_phi(-2:), ipgex_th(0:,-2:)
  
end subroutine COCAL_ID_interpo_lag4th_2Dsurf
  end interface
end module COCAL_ID_interface_interpo_lag4th_2Dsurf
