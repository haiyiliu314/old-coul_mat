program potential
  implicit none
  integer, parameter                        ::Na = 1, Nb = 1, Nc = 2, N = Na + Nb + Nc, N1 = 100
  integer, parameter                        ::LDA = N, LDVL = N, LDVR = N
  integer, parameter                        ::LWMAX = 10000
  double precision, parameter               ::R = 100, Ra = 2, Rb = 12.5, Rc = R - Ra - Rb
  double precision, parameter               ::dNa = Ra / Na, dNb = Rb / Nb, dNc = Rc / Nc
  integer                                   ::INFO, LWORK,i,j=0
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
  double precision, dimension(N1)           ::a1,b,c
  double precision, dimension(N)            ::d
  character(80)                             ::list_file

A=1
VL(1,:)=A(1,:)
write(*,*) VL
end program

