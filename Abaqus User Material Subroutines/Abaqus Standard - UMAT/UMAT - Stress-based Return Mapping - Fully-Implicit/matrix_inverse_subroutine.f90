!=======================================================================
! User Material Subroutine for Isotropic Elasto-Visco-Plasticity
! CDD-based Material Model
! Fully-Implicit Constitutive Integration
!
! Author: Seyed Amir Hossein Motaman
! Steel Institute (IEHK), RWTH Aachen University
!
! References:
!    Motaman, S.A.H.; Prahl, U.; 2019.
!    Microstructural constitutive model for polycrystal viscoplasticity
!    in cold and warm regimes based on continuum dislocation dynamics.
!    Journal of the Mechanics and Physics of Solids 122, 205–243.
!    doi: https://doi.org/10.1016/j.jmps.2018.09.002
!
!    Motaman, S.A.H.; Schacht K.; Haase, C.; Prahl, U.; 2019.
!    Thermo-micro-mechanical simulation of metal forming processes.
!    International Journal of Solids and Structures.
!    doi: https://doi.org/10.1016/j.ijsolstr.2019.05.028
!=======================================================================

!***********************************************************************
! subroutine header
  subroutine matrix_inverse(A,C,n)

!-----------------------------------------------------------------------
!    use of global variables
     use controls, only: &
         pInt, &
         pReal

!-----------------------------------------------------------------------
!    declaration of subroutine's variables
     implicit none

     integer(pInt),               intent(in) :: &
         n                                                               ! dimension of input matrix to be inversed

     real(pReal), dimension(n,n), intent(in) :: &
         A                                                               ! input matrix

     real(pReal), dimension(n,n), intent(out) :: &
         C                                                               ! inverse of input matrix

!-----------------------------------------------------------------------
!    declaration of subroutines' local variables
     integer(pInt) :: &
         i, &
         j, &
         k

     real(pReal) :: &
         e

     real(pReal), dimension(n) :: &
         b, &
         d, &
         x

     real(pReal), dimension(n,n) :: &
         A_prime, &
         L, &
         U

!-----------------------------------------------------------------------
!    step 0: initializing for matrices L and U and b
     A_prime = A
     L       = 0.0e0
     U       = 0.0e0
     b       = 0.0e0

!    step 1: forward elimination
     do k = 1, n - 1
         do i = k + 1, n
             e = A_prime(i,k) / A_prime(k,k)
             L(i,k) = e
             do j = k + 1, n
                 A_prime(i,j) = A_prime(i,j) - e * A_prime(k,j)
             end do
         end do
     end do

!-----------------------------------------------------------------------
!    step 2: preparing L and U matrices 
!    L matrix is a matrix of the elimination coefficient with the diagonal elements equal to 1.0
     do i = 1, n
         L(i,i) = 1.0e0
     end do
!    calculating U matrix as the upper triangular part of A_prime
     do j = 1, n
         do i = 1, j
             U(i,j) = A_prime(i,j)
         end do
     end do

!-----------------------------------------------------------------------
!    step 3: computing columns of the inverse matrix C
     do k = 1, n
         b(k) = 1.0e0
         d(1) = b(1)
!        step 3a: solving Ld=b using the forward substitution
         do i = 2, n
             d(i) = b(i)
             do j = 1, i - 1
                 d(i) = d(i) - L(i,j) * d(j)
             end do
         end do
!        step 3b: solving Ux=d using the back substitution
         x(n) = d(n) / U(n,n)
         do i = n - 1, 1, -1
             x(i) = d(i)
             do j = n, i + 1, -1
                 x(i) = x(i) - U(i,j) * x(j)
             end do
             x(i) = x(i) / U(i,i)
         end do
!        step 3c: filling the solutions x(n) into column k of C
         do i = 1, n
             C(i,k) = x(i)
         end do
         b(k) = 0.0e0
     end do

  end subroutine matrix_inverse