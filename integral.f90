program integral
  implicit none
  integer, parameter                                ::N=100
  double precision, dimension(N)                    ::a,b,c
  double precision                                  ::y=2,y1=1,s
  integer                                           ::i
  real, parameter                                   ::pi=4*atan(1.)
  do i=1,N
    a(i)=i
  end do
  b=cos((2*a-1)*pi/(2*N))
  s=sum(1/(sqrt(y**2+y1**2-2*y1*y*b)))*pi/N*2
  write(*,*) s
end program integral
