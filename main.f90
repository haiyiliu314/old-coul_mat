! Insert Sort algorithm
!


program main
  implicit none

  integer, parameter            :: N = 10
  integer, dimension(N)         :: list
  integer, dimension(N)         :: list_copy

  integer                       :: i, j
  integer                       :: temp_storage

  ! Output file variables
 character(80) :: list_file

! --------------------------------------------------------------------

  list = (/ 4, 5, 2, 8, 9, 7, 6, 1, 3, 0 /)
  list_copy = list

  ! --- Unsorted list output --- !
  write(*,*) 'write file Unsorted.dat'
  write(list_file, '(A)') 'Unsorted.dat'
  open(unit=700,file=list_file)

  DO i = 1, N
    write(700,*)   list(i)
  END DO
  close(700)


  DO i = 2, N
    j = i
    DO WHILE ((j > 1) .AND. (list(j-1) < list(j)))
      temp_storage = list(j-1)
      list(j-1) = list(j)
      list(j) = temp_storage
      j = j-1
    END DO
    write(*,*) list
  END DO



  ! --- Sorted list output --- !
  write(*,*) 'write file Sorted.dat'
  write(list_file, '(A)') 'Sorted.dat'
  open(unit=700,file=list_file)

  DO i = 1, N
    write(700,*)   list(i)
  END DO
  close(700)

  

end program



