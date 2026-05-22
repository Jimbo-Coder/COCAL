subroutine COCAL_ID_lagint_4th_weights(x,v,w)

use COCAL_ID_phys_constant, only : long
implicit none
  real(long), intent(in)  :: x(4), v
  real(long), intent(out) :: w(4)
  real(long) :: dx12, dx13, dx14, dx23, dx24, dx34
  real(long) :: dx21, dx31, dx32, dx41, dx42, dx43
  real(long) :: xv1, xv2, xv3, xv4
!
      dx12 = x(1) - x(2)
      dx13 = x(1) - x(3)
      dx14 = x(1) - x(4)
      dx23 = x(2) - x(3)
      dx24 = x(2) - x(4)
      dx34 = x(3) - x(4)
      dx21 = - dx12
      dx31 = - dx13
      dx32 = - dx23
      dx41 = - dx14
      dx42 = - dx24
      dx43 = - dx34
      xv1 = v - x(1)
      xv2 = v - x(2)
      xv3 = v - x(3)
      xv4 = v - x(4)
      w(1) = xv2*xv3*xv4/(dx12*dx13*dx14)
      w(2) = xv1*xv3*xv4/(dx21*dx23*dx24)
      w(3) = xv1*xv2*xv4/(dx31*dx32*dx34)
      w(4) = xv1*xv2*xv3/(dx41*dx42*dx43)
!
end subroutine COCAL_ID_lagint_4th_weights

function COCAL_ID_lagint_4th_apply(w,y)

use COCAL_ID_phys_constant, only : long
implicit none
  real(long) :: COCAL_ID_lagint_4th_apply
  real(long), intent(in) :: w(4), y(4)
!
      COCAL_ID_lagint_4th_apply = w(1)*y(1) + w(2)*y(2) + w(3)*y(3) + w(4)*y(4)
!
end function COCAL_ID_lagint_4th_apply
