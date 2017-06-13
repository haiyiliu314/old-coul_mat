program eigen
  implicit none
      INTEGER,parameter          ::NIN=5, NOUT=6
!      PARAMETER        (NIN=5,NOUT=6)
      INTEGER,parameter          ::NB=64, NMAX=10
!      PARAMETER        (NB=64,NMAX=10)
      INTEGER,parameter          ::LDA=NMAX, LDVR=NMAX, LWORK=(1+NB)*NMAX
!      PARAMETER        (LDA=NMAX,LDVR=NMAX,LWORK=(1+NB)*NMAX)

      INTEGER                    ::I, IFAIL, INFO, J, LWKOPT, N
      
      COMPLEX *16                ::A(LDA,NMAX), DUMMY(1,1), VR(LDVR,NMAX), W(NMAX),&
                                   WORK(LWORK)

      DOUBLE PRECISION           ::RWORK(2*NMAX)
      N=10
      EXTERNAL         ZGEEV
      CALL ZGEEV('No left vectors','Vectors (right)',N,A,LDA,W,DUMMY,  &                        

      1,VR,LDVR,WORK,LWORK,RWORK,INFO)
end program eigen
