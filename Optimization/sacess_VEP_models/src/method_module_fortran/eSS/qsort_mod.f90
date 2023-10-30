! Recursive Fortran 95 quicksort routine
! sorts real numbers into ascending numerical order
! Author: Juli Rew, SCD Consulting (juliana@ucar.edu), 9/03
! Based on algorithm from Cormen et al., Introduction to Algorithms,
! 1997 printing

! Made F conformant by Walt Brainerd

module qsort_module
    USE scattersearchtypes
    IMPLICIT NONE
public :: QsortC
private :: Partition

type :: valuepos
    REAL(KIND=SELECTED_REAL_KIND(P=PRECISION_D,R=RANGE_D)) :: value
    integer :: position
end type valuepos

contains

recursive subroutine QsortC_ind(A,ind)
  REAL(KIND=SELECTED_REAL_KIND(P=PRECISION_D,R=RANGE_D)), intent(in out), dimension(:) :: A
  integer, intent(in out), dimension(:) :: ind
  REAL(KIND=SELECTED_REAL_KIND(P=PRECISION_D,R=RANGE_D)) :: temp
  type(valuepos), dimension(:), ALLOCATABLE :: data1
  integer :: i

  if (size(A) .GT. 2 ) then
        ALLOCATE(data1(size(A)))
  
        do i=1,size(A)
                data1(i) % value = A(i)
                data1(i) % position = i
        end do

        call QsortC(data1)

        do i=1,size(A)
                A(i) = data1(i) % value
                ind(i) = data1(i) % position
        end do

        DEALLOCATE(data1)
  else
      
        if (A(1) .LT. A(2) ) then
               ind(1) = 1
               ind(2) = 2  
        else
               temp = A(1)
               A(1) = A(2)
               A(2) = temp
               ind(1) = 2
               ind(2) = 1
        end if
  end if

end subroutine QsortC_ind

recursive subroutine QsortC(data1)
  type(valuepos), intent(in out), dimension(:) :: data1
  integer :: iq

  if(size(data1) > 1) then
     call Partition(data1, iq)
     call QsortC(data1(:iq-1))
     call QsortC(data1(iq:))
  endif
end subroutine QsortC


subroutine Partition(data1, marker)
  type(valuepos), intent(in out), dimension(:) :: data1
  integer, intent(out) :: marker
  integer :: i, j
  type(valuepos) :: temp
  type(valuepos) :: x      ! pivot point
  x = data1(1)
  i= 0
  j= size(data1) + 1

  do
     j = j-1
     do
        if (data1(j)%value <= x%value) exit
        j = j-1
     end do
     i = i+1
     do
        if (data1(i)%value >= x%value) exit
        i = i+1
     end do
     if (i < j) then
        ! exchange A(i) and A(j)
        temp = data1(i)
        data1(i) = data1(j)
        data1(j) = temp
     elseif (i == j) then
        marker = i+1
        return
     else
        marker = i
        return
     endif
  end do

  end subroutine Partition

end module qsort_module
