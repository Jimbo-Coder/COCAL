! One copy per node of the COCAL data (cocal_id_node_shared.c). The node
! leader reads the COCAL arrays and fits the metric spline into memory that
! every rank of the node maps; the other ranks only take the scalar constants
! from the leader. Without MPI, with COCAL_ID_share_node_memory = no, or when
! the MPI cannot create shared windows, every rank is its own leader.
module COCAL_ID_node_shared
  use iso_c_binding
  implicit none
  private
  public :: COCAL_ID_node_leader, COCAL_ID_node_barrier, COCAL_ID_node_bcast
  public :: COCAL_ID_node_alloc_2d, COCAL_ID_node_alloc_3d, COCAL_ID_node_alloc_4d
  public :: COCAL_ID_node_free_all

  interface
    integer(c_int) function node_leader_c() bind(c, name="COCAL_ID_node_leader")
      import :: c_int
    end function node_leader_c
    type(c_ptr) function node_alloc_c(n) bind(c, name="COCAL_ID_node_alloc")
      import :: c_ptr, c_long
      integer(c_long), value :: n
    end function node_alloc_c
    subroutine COCAL_ID_node_barrier() bind(c, name="COCAL_ID_node_barrier")
    end subroutine COCAL_ID_node_barrier
    subroutine node_bcast_c(buf, n) bind(c, name="COCAL_ID_node_bcast")
      import :: c_double, c_int
      real(c_double) :: buf(*)
      integer(c_int), value :: n
    end subroutine node_bcast_c
    subroutine COCAL_ID_node_free_all() bind(c, name="COCAL_ID_node_free_all")
    end subroutine COCAL_ID_node_free_all
  end interface

contains

  logical function COCAL_ID_node_leader()
    COCAL_ID_node_leader = node_leader_c() /= 0
  end function COCAL_ID_node_leader

  ! Collective over the node: a(0:n1,0:n2), the COCAL surface layout, zeroed.
  subroutine COCAL_ID_node_alloc_2d(a, n1, n2)
    real(8), pointer, intent(out) :: a(:,:)
    integer, intent(in) :: n1, n2
    real(8), pointer :: flat(:)
    integer(c_long) :: n
    n = int(n1+1,c_long)*(n2+1)
    call c_f_pointer(node_alloc_c(n), flat, [n])
    a(0:n1,0:n2) => flat
    if (COCAL_ID_node_leader()) a = 0
  end subroutine COCAL_ID_node_alloc_2d

  ! Collective over the node: a(0:n1,0:n2,0:n3), the COCAL grid layout, zeroed.
  subroutine COCAL_ID_node_alloc_3d(a, n1, n2, n3)
    real(8), pointer, intent(out) :: a(:,:,:)
    integer, intent(in) :: n1, n2, n3
    real(8), pointer :: flat(:)
    integer(c_long) :: n
    n = int(n1+1,c_long)*(n2+1)*(n3+1)
    call c_f_pointer(node_alloc_c(n), flat, [n])
    a(0:n1,0:n2,0:n3) => flat
    if (COCAL_ID_node_leader()) a = 0
  end subroutine COCAL_ID_node_alloc_3d

  ! Collective over the node: a(1:nq,0:nr,0:nt,0:np-1), the spline layout.
  ! Not zeroed: the caller fills every element.
  subroutine COCAL_ID_node_alloc_4d(a, nq, nr, nt, np)
    real(8), pointer, intent(out) :: a(:,:,:,:)
    integer, intent(in) :: nq, nr, nt, np
    real(8), pointer :: flat(:)
    integer(c_long) :: n
    n = int(nq,c_long)*(nr+1)*(nt+1)*np
    call c_f_pointer(node_alloc_c(n), flat, [n])
    a(1:nq,0:nr,0:nt,0:np-1) => flat
  end subroutine COCAL_ID_node_alloc_4d

  subroutine COCAL_ID_node_bcast(x)
    real(8), intent(inout) :: x
    real(8) :: buf(1)
    buf(1) = x
    call node_bcast_c(buf, 1_c_int)
    x = buf(1)
  end subroutine COCAL_ID_node_bcast

end module COCAL_ID_node_shared
