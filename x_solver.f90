! X_SOLVER.F90=====================================================================================
!
!
! =================================================================================================
! exciton solver for homogeneous problem
! find eigenvalues and -states
! =================================================================================================

subroutine x_solver(N_k, pi, k, V_b, E_1s)
    implicit none

    ! input variables
    integer, intent(in) :: N_k
    double precision, intent(in) :: pi
    double precision, dimension(N_k), intent(in) :: k
    double precision, dimension(N_k,N_k), intent(in) :: V_b
    double precision, intent(out) :: E_1s


    ! matrix to invert
    double precision, dimension(N_k,N_k) :: mat_m

    ! help variables to solve eigenvalue problem
    character :: JOBVL
    character :: JOBVR
    integer :: N
    double precision, allocatable, dimension(:,:) :: A
    integer :: LDA
    double precision, allocatable, dimension(:) :: WR
    double precision, allocatable, dimension(:) :: WI
    double precision, allocatable, dimension(:,:) :: VL
    integer :: LDVL
    double precision, allocatable, dimension(:,:) :: VR
    integer :: LDVR
    double precision, allocatable, dimension(:) :: WORK
    integer :: LWORK
    integer :: INFO


    ! loop indices
    integer :: i1
    integer :: i2



    ! create matrix to invert
    do i1 = 1, N_k
        do i2 = 1, N_k
            if (i1 == i2) then
                mat_m(i1,i2) = 2.0d0*pi*k(i1)**2 - V_b(i2,i1)
            else
                mat_m(i1,i2) = - V_b(i2,i1)
            end if
        end do
    end do


    ! set up parameters for matrix inversion
    ! compute left eigenvectors
    JOBVL = 'V'
    ! compute rigth eigenvectors
    JOBVR = 'V'
    ! order of matrizes
    N = N_k
    LDA = N
    ! matrix to invert
    allocate(A(LDA,N))
    A = mat_m/(2.0d0*pi)
    ! eigenvalues
    allocate(WR(N), WI(N))
    ! eigenvectors
    LDVL = N
    LDVR = N
    allocate(VL(LDVL,N), VR(LDVR,N))
    ! WORK
    LWORK = 4*N
    allocate(WORK(LWORK))


    ! compute eigenvalue problem
    call DGEEV(JOBVL, JOBVR, N, A, LDA, WR, WI, VL, LDVL, VR, LDVR, WORK, LWORK, INFO)

    ! info of DGEEV
    write(*,'("X solver (DGEEV) INFO=", I5)') INFO

    if (INFO == 0) then
        ! output eigenvalues (real part, imaginary part is zero for Hermitian matrices)
        open(unit=42, file='ev_WR.dat')
        ! open(unit=43, file='ev_WI.dat')
        do i1 = 1, N
            write(42, *) WR(i1)
        !     write(43, *) WI(i1)
        end do
        close(42)
        ! close(43)

        ! output eigenvectors
        open(unit=42, file='ev_VR.dat', access='stream')
        write(42) VR
        close(unit=42)
        open(unit=42, file='ev_VL.dat', access='stream')
        write(42) VL
        close(unit=42)


        ! find 1s binding energy as minimum of eigenvalues
        E_1s = minval(WR)
        write(*,'("E_1s / E_b = ", SE13.6E3)') E_1s
    end if


end subroutine x_solver
