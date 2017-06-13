program main
  implicit none
  integer, parameter                        ::Na = 50, Nb = 50, Nc = 20, N = Na + Nb + Nc, N1 = 100, N2 = 100
  integer, parameter                        ::LDA = N, LDVL = N, LDVR = N
  integer, parameter                        ::LWMAX = 10000
  double precision, parameter               ::R = 100, Ra = 2, Rb = 4, Rc = R - Ra - Rb
  double precision, parameter               ::dNa = Ra / Na, dNb = Rb / Nb, dNc = Rc / Nc
  integer                                   ::INFO, LWORK,i,j=0, i1=0
  double precision                          ::RWORK(2*N), y( N )=0, y1( N ),s2( N ), y2(N1)
  complex*16                                ::A( LDA, N ), VL( LDVL, N ), VR( LDVR, N ),&
                                              W( N ), WORK( LWMAX )
  real*16, dimension(:), allocatable     ::eigen
  character                                 ::YES = 'V',NO = 'N'
  real*16,parameter                         ::m = 9.10938356d-31, epsilon = 8.854187817d-12,&
                                              hbar = 1.054571800d-34, e = 1.60217662d-19, &
                                              pi = 4*atan(1.)
  real*16, parameter                        ::aB = 4 * pi * epsilon * hbar ** 2/(m * e ** 2)
  real*16, parameter                        ::eB = e ** 2 / ( 8 * pi * epsilon * aB )
  double precision, dimension(N1)           ::a1,b,c, w1
  double precision, dimension(N,N)          ::d
  character(80)                             ::list_file

  do i= 1, Na
    y(i) = dNa * ( i - 0.5 );
  end do
  do i= ( Na + 1 ), Na + Nb
    y(i) = dNa * Na + dNb * ( i - Na - 0.5);
  end do
  do i = ( Na + Nb + 1 ), N
    y(i) = dNa * Na + dNb * Nb + dNc * (i-Na-Nb-0.5);
  end do
!  y=y
  y1=y
  LWORK = -1

!  write(list_file, '(A)') 'y.dat'
!  open(unit=701, file=list_file)
!  write(701,*) y
!  close(701)
!  call potential(N,y,y1,A)

  do i=1,N1
    a1(i)=i
  end do
b = -1 + 2. / (N1 + 1. ) * ( a1 - ( N1 + 1. ) / ( 2. * pi )*sin(2.*a1*pi/(N1+1.)))
w1 = 2. / ( N1 + 1. ) * ( 1. - cos( 2. * pi * a1 / ( N1 + 1.))) 
!write(*,*) b
  do j=1,Na
    y2=y1(j)-dNa/2+a1*dNa/N1
    do i=1,N1
      d(j,:)=d(j,:)+(1/(sqrt(y(j)**2+y1**2-2*y1*y(j)*b(i))))*2*y1/(sqrt( 1 - b(i) ** 2))*w1(i)
    end do
    d(j,j)=0
    do i1=1,N2
      d(j,j)=d(j,j) + sum( 1 / sqrt( y(j) ** 2 + y2(i1) ** 2 - 2 * y2(i1) * y(j) * b ) / N2 * 2 * y2(i1) / ( sqrt(1 - b ** 2 )) * w1)
    end do
!    write(*,*), d*y
    A(j,1:Na)=-cmplx(d(j,1:Na)/pi*dNa,0)
    A(j,(Na+1):(Na+Nb))=-cmplx(d(j,(Na+1):(Na+Nb))/pi*dNb,0)
    A(j,(Na+Nb+1):(N))=-cmplx(d(j,(Na+Nb+1):(N))/pi*dNc,0)
  end do

  do j=(Na+1), (Na+Nb)
      y2=y1(j)-dNb/2+a1*dNb/N1
    do i=1,N1
      d(j,:)=d(j,:)+(1/(sqrt(y(j)**2+y1**2-2*y1*y(j)*b(i))))*2*y1/(sqrt( 1 - b(i) ** 2))*w1(i)
    end do
    d(j,j)=0
    do i1=1,N2
      d(j,j)=d(j,j) + sum( 1 / sqrt( y(j) ** 2 + y2(i1) ** 2 - 2 * y2(i1) * y(j) * b ) / N2 * 2 * y2(i1) / ( sqrt(1 - b ** 2 )) * w1)
    end do
!    write(*,*), d*y
    A(j,1:Na)=-cmplx(d(j,1:Na)/pi*dNa,0)
    A(j,(Na+1):(Na+Nb))=-cmplx(d(j,(Na+1):(Na+Nb))/pi*dNb,0)
    A(j,(Na+Nb+1):(N))=-cmplx(d(j,(Na+Nb+1):(N))/pi*dNc,0)
  end do

  do j=(Na+Nb+1),N
      y2=y1(j)-dNc/2+a1*dNc/N1
    do i=1,N1
      d(j,:)=d(j,:)+(1/(sqrt(y(j)**2+y1**2-2*y1*y(j)*b(i))))*2*y1/(sqrt( 1 - b(i) ** 2))*w1(i)
    end do
    d(j,j)=0
    do i1=1,N2
      d(j,j)=d(j,j) + sum( 1 / sqrt( y(j) ** 2 + y2(i1) ** 2 - 2 * y2(i1) * y(j) * b ) / N2 * 2 * y2(i1) / ( sqrt(1 - b ** 2 )) * w1)
    end do
!    write(*,*), d*y
    A(j,1:Na)=-cmplx(d(j,1:Na)/pi*dNa,0)
    A(j,(Na+1):(Na+Nb))=-cmplx(d(j,(Na+1):(Na+Nb))/pi*dNb,0)
    A(j,(Na+Nb+1):(N))=-cmplx(d(j,(Na+Nb+1):(N))/pi*dNc,0)
  end do
!  A=A*dN
  do i=1,N
      A(i,i)=A(i,i)+y(i)**2
  end do
  write(*,*) b
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

