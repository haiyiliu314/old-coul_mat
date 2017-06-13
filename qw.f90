program main
  implicit none
  integer, parameter                        ::N = 1000, N1 = 100, N2 = 50
  integer, parameter                        ::LDA = N, LDVL = N, LDVR = N
  integer, parameter                        ::LWMAX = 10000
  double precision                          ::dN= 0, ymax = 100d0
  integer                                   ::INFO, LWORK,i,j=0, i1
  double precision                          ::RWORK(2*N), y( N )=0, y1( N ),s2( N )
  complex*16                                ::VL( LDVL, N ), VR( LDVR, N ),&
                                              W( N ), WORK( LWMAX ), A( LDA, N )=0 
  real*16, dimension(:), allocatable     ::eigen
  character                                 ::YES = 'V',NO = 'N'
  real*16,parameter                         ::m = 9.10938356d-31, epsilon = 8.854187817d-12,&
                                              hbar = 1.054571800d-34, e = 1.60217662d-19, &
                                              pi = 4*atan(1.)
  real*16, parameter                        ::aB = 4 * pi * epsilon * hbar ** 2/(m * e ** 2)
  real*16, parameter                        ::eB = e ** 2 / ( 8 * pi * epsilon * aB )
  character(80)                                  ::list_file
  character(len=100)                             ::format_V
  complex*16                                     ::test_arr(2,2)
  double precision                               :: test_arr1(2, 1)  
  dN=ymax/N
  do i=1,N
    y(i)=dN*i;
  end do
!  y=y
  y1=y
  LWORK = -1
  call potential(N,y,y1,A, 100)
  A=A*dN
  do i=1,N
      A(i,i)=A(i,i)+y(i)**2
  end do
  test_arr = 1;
  test_arr1 = 2;

  write(*,*) matmul(test_arr, test_arr1)
 ! write(*,*)  s2
  CALL ZGEEV(NO, YES, N, A, LDA, W, VL, LDVL,&
             VR, LDVR, WORK, LWORK, RWORK, INFO )
  LWORK = min( LWMAX, int( WORK( 1 ) ) )
  CALL ZGEEV(NO, YES, N, A, LDA, W, VL, LDVL,&
             VR, LDVR, WORK, LWORK, RWORK, INFO )
!  write(*,*)  aB,eB
  write(*,*)  "eigenvalue="
!  write(*,*)  W
  do i=1,N
    if(real(W(i))<0) then
    j=j+1
    end if 
  end do
  allocate(eigen(j))
  j=1
  do i=1,N
    if(real(W(i))<0) then
    eigen(j)=real(W(i))
    j=j+1
    end if 
  end do
  write(*,*) eigen
end program
subroutine potential(N1, y, y1, s1, N)
  implicit none
  integer, intent(in)                               ::N
  double precision, dimension(N)                    ::a,b,c
  double precision, dimension(N1)                   ::d
  integer,intent(in)                                ::N1
  double precision,dimension(N1), intent(in)        ::y, y1
  complex*16,dimension(N1,N1), intent(out)          ::s1
  integer                                           ::i,j
  real, parameter                                   ::pi=4*atan(1.)
  do j=1,N1
    d=0
    do i=1,N
      a(i)=i
      b(i)=cos((2*a(i)-1)*pi/(2*N))
      d=d+(1/(sqrt(y(j)**2+y1**2-2*y1*y(j)*b(i))))*pi/N*2
    end do
!    write(*,*), d*y
    s1(j,1:N1)=-cmplx(d*y/pi,0)
  end do
!write(*,*), N1, s1
end subroutine potential
