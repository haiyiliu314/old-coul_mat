! INIT_COULOMB_MATRIX.F90==========================================================================
!
!
! =================================================================================================
! subroutine which computes the Coulomb matrix elements for the later Coulomb sums
!
! =================================================================================================
! INPUT:
! none
!
! =================================================================================================
! OUTPUT
! none
!
! =================================================================================================
! CHANGES

subroutine init_coulomb_matrix
    use params
    implicit none

    interface
        subroutine init_phi_grid(gauss_tscheby, N_phi, phi_min, phi_max, phi, dphi)
            use params, only: pi
            implicit none
            integer, intent(in) :: gauss_tscheby
            integer, intent(in) :: N_phi
            double precision, intent(in) :: phi_min
            double precision, intent(in) :: phi_max
            double precision, dimension(N_phi), intent(out) :: phi
            double precision, dimension(N_phi), intent(out) :: dphi
        end subroutine init_phi_grid
    end interface

    ! variables for phi grid
    double precision, dimension(N_phi) :: phi
    double precision, dimension(N_phi) :: dphi

    ! variables for fine grid
    double precision :: k_fine_min
    double precision :: k_fine_max
    double precision :: dk_fine
    double precision :: dk_fine2
    double precision, dimension(N_k_fine) :: k_fine
    double precision, dimension(N_k_fine2) :: k_fine2

    ! help variable to compute Coulomb matric elements
    double precision :: V_tmp

    ! loop indices
    integer :: i_1
    integer :: i_2
    integer :: j_1
    integer :: j_2



    !==================================================================================================
    ! initialize phi grid for phi integration
    call init_phi_grid(sw_phi_grid_property, N_phi, phi_min, phi_max, phi, dphi)


    !==================================================================================================
    ! compute phi integrated Coulomb matrix elements
    V_mat = 0.0d0


    if (sw_coulomb /= 0) then

        !============================================================================================
        ! compute Coulomb matrix without confinement (i_1=Q', i_2=Q)
        !$OMP PARALLEL DO PRIVATE(dk_fine, k_fine, dk_fine2, k_fine2, V_tmp, j_1, j_2) COLLAPSE(2)
        do i_2 = 1, N_k
            do i_1 = 1, N_k

                if (sw_coulomb_diagonal == 0) then
                    ! use two averaging fine sums (no averaging for N_k_fine=1)
                    dk_fine = dk(i_1) / dble(N_k_fine)
                    k_fine = (/ (k(i_1)-0.5d0*dk(i_1) + dk_fine*(dble(j_1)-0.5d0), j_1=1, N_k_fine) /)
                    dk_fine2 = dk(i_2) / dble(N_k_fine2)
                    k_fine2 = (/ (k(i_2)-0.5d0*dk(i_2) + dk_fine2*(dble(j_2)-0.5d0), j_2=1, N_k_fine2) /)

                    ! compute Coulomb matrix element
                    V_tmp = 0.0d0
                    do j_1 = 1, N_k_fine
                        do j_2 = 1, N_k_fine2
                            V_tmp = V_tmp + 4.0d0 * dk_fine2*dk_fine/dk(i_1)/dk(i_2) * k_fine(j_1) * sum(dphi/dsqrt(k_fine(j_1)**2 + k_fine2(j_2)**2 - 2.0d0*k_fine(j_1)*k_fine2(j_2)*dcos(phi) + eps_screen**2))
                        end do
                    end do

                    ! include quarature factor dk into Coulomb elements because they always appear together with the Coulomb elements
                    V_mat(i_1,i_2) = V_tmp * dk(i_1)

                else if (sw_coulomb_diagonal == 1) then
                    if (i_1 == i_2) then
                        ! use two averaging fine sums (no averaging for N_k_fine=1)
                        dk_fine = dk(i_1) / dble(N_k_fine)
                        k_fine = (/ (k(i_1)-0.5d0*dk(i_1) + dk_fine*(dble(j_1)-0.5d0), j_1=1, N_k_fine) /)
                        dk_fine2 = dk(i_2) / dble(N_k_fine2)
                        k_fine2 = (/ (k(i_2)-0.5d0*dk(i_2) + dk_fine2*(dble(j_2)-0.5d0), j_2=1, N_k_fine2) /)

                        ! compute Coulomb matrix element
                        V_tmp = 0.0d0
                        do j_1 = 1, N_k_fine
                            do j_2 = 1, N_k_fine2
                                V_tmp = V_tmp + 4.0d0 * dk_fine2*dk_fine/dk(i_1)/dk(i_2) * k_fine(j_1) * sum(dphi/dsqrt(k_fine(j_1)**2 + k_fine2(j_2)**2 - 2.0d0*k_fine(j_1)*k_fine2(j_2)*dcos(phi) + eps_screen**2))
                            end do
                        end do
                    else
                        V_tmp = 4.0d0 * k(i_1) * sum(dphi/dsqrt(k(i_1)**2 + k(i_2)**2 - 2.0d0*k(i_1)*k(i_2)*dcos(phi) + eps_screen**2))
                    end if

                    ! include quarature factor dk into Coulomb elements because they always appear together with the Coulomb elements
                    V_mat(i_1,i_2) = V_tmp * dk(i_1)


                else if (sw_coulomb_diagonal == 2) then
                    V_tmp = 4.0d0 * k(i_1) * sum(dphi/dsqrt(k(i_1)**2 + k(i_2)**2 - 2.0d0*k(i_1)*k(i_2)*dcos(phi) + eps_screen**2))

                    ! include quarature factor dk into Coulomb elements because they always appear together with the Coulomb elements
                    V_mat(i_1,i_2) = V_tmp * dk(i_1)
                end if

            end do
        end do
        !$OMP END PARALLEL DO

    end if


end subroutine init_coulomb_matrix
