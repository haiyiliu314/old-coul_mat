program arraypass
implicit none

	integer, dimension (5) :: myarray
	integer :: i 

	call fillArray (myarray)
	call printArray (myarray)
        write(*,*) myarray  
end program arraypass

subroutine fillArray (a)
implicit none
	
	integer, dimension (5), intent (out) :: a

	! local variables
	integer :: i
	do i = 1, 5
		a(i) = i
	end do

end subroutine fillArray

subroutine printArray (a)
	
	integer, dimension (5) :: a
	integer :: i

	do i = 1, 5
		Print *, a(i)
	end do

end subroutine printArray
