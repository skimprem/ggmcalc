! subroutine square_array(n, arr, a)
!     implicit none
!     integer, intent(in) :: n
!     real(8), intent(in) :: a
!     real(8), intent(inout) :: arr(n)
!     integer :: i
!     do i = 1, n
!         arr(i) = arr(i)**2
!     end do
!     print *, arr
!     print *, n
!     print *, a
! end subroutine square_array

subroutine square_array(arr, a, s)

    implicit none

    type mytype
        integer(4) :: ii
        real(8) :: rr
        character(10) :: ss
    end type

    real(8), intent(in) :: a
    real(8), intent(inout) :: arr(:)
    character(len=10), intent(in) :: s
    integer :: i
    type(mytype) :: conf

    do i = 1, size(arr)
        arr(i) = arr(i)**2
    end do

    print *, arr
    print *, a
    print *, s

    print *, conf

end subroutine square_array
