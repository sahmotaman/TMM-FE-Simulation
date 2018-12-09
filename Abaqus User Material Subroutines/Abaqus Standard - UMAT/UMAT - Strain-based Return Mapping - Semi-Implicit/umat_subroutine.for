!=======================================================================
! User Material Subroutine for Isotropic Elasto-Visco-Plasticity
! CDD-based Material Model
! Semi-Implicit Constitutive Integration with Strain-based Return Mapping
!
! Author: Seyed Amir Hossein Motaman
! Steel Institute (IEHK), RWTH Aachen University
!
! References:
!    Motaman, S.A.H.; Prahl, U.; 2019.
!    Microstructural constitutive model for polycrystal viscoplasticity
!    in cold and warm regimes based on continuum dislocation dynamics.
!    Journal of the Mechanics and Physics of Solids 122, 205–243.
!    doi: 10.1016/j.jmps.2018.09.002
!
!    Motaman, S.A.H.; Schacht K.; Haase, C.; Prahl, U.; 2019.
!    Thermo-micro-mechanical simulation of bulk metal forming processes.
!=======================================================================

!***********************************************************************
! including other modules and subroutines
#include 'controls_module.f90'
#include 'strings_module.f90'
#include 'material_properties_module.f90'
#include 'plasticity_subroutine.f90'


!***********************************************************************
! subroutine header
  subroutine umat(sigma,      sv,          C,         sse,       delta_w_p, &
                  scd,        q_dot_p,     dsigma_dT, dq_deps,   dq_dT, &
                  eps,        delta_eps,   t,         delta_t,   T_n, &
                  delta_temp, pfv,         delta_pfv, cmname,    n_dir, &
                  n_shr,      n_total,     n_sv,      props,     n_props, &
                  coords,     delta_R,     pnewdt,    l_element, F_n, &
                  F,          i_el,        i_ip,      i_layer,   i_sp, &
                  i_step,     i_increment)

!-----------------------------------------------------------------------
! pre-processing
!-----------------------------------------------------------------------
!    use of global variables
     use controls
     use material_properties

!-----------------------------------------------------------------------
!    declaration of subroutine's variables
     implicit none

     character(len = 64),                     intent(in) :: &
         cmname                                                          ! uses-specified material name, left justified

     integer(pInt),                           intent(in) :: &
         n_dir, &                                                        ! number of independent direct components of strain and stress tensors
         n_shr, &                                                        ! number of independent shear  components of strain and stress tensors
         n_total, &                                                      ! number of independent total  components of strain and stress tensors
         n_sv, &                                                         ! number of user-defined state variables
         n_props, &                                                      ! number of user-defined material parameters
         i_el, &                                                         ! element number
         i_ip, &                                                         ! integration point number
         i_layer, &                                                      ! layer number
         i_sp, &                                                         ! section point number within the current layer
         i_step, &                                                       ! step number
         i_increment                                                     ! increment number

     real(pReal),                             intent(in) :: &
         delta_t, &                                                      ! time increment [s]
         T_n, &                                                          ! temperature at the beginning of the time increment [C]
         delta_temp, &                                                   ! increment of temperature [C]
         l_element                                                       ! characteristic element length [mm]

     real(pReal), dimension(1),               intent(in) :: &
         pfv, &                                                          ! array of interpolated values of predefined field variables at the beginning of the time increment  
         delta_pfv                                                       ! array of increments of predefined field variables

     real(pReal), dimension(2),               intent(in) :: &
         t                                                               ! t(1): step time and t(2): total time, at the beginning of the time increment [s]

     real(pReal), dimension(3),               intent(in) :: &
         coords                                                          ! coordinates of the integration point [mm]

     real(pReal), dimension(n_total),         intent(in) :: &
         eps, &                                                          ! (total) strain tensor (voigt notation, vector form) at beginning of the time increment [-]
         delta_eps                                                       ! (total) strain increment tensor (voigt notation, vector form) [-]

     real(pReal), dimension(n_props),         intent(in) :: &
         props                                                           ! material constants array

     real(pReal), dimension(3,3),             intent(in) :: &
         delta_R, &                                                      ! Incremental rotation matrix [-]
         F_n, &                                                          ! deformation gradient tensor at the beginning of the time increment [-]
         F                                                               ! deformation gradient tensor at the end       of the time increment [-]

     real(pReal),                             intent(inout) :: &
         sse, &                                                          ! specific elastic strain energy
         delta_w_p, &                                                    ! specific plastic dissipation
         scd, &                                                          ! specific creep dissipation
         q_dot_p, &                                                      ! volumetric heat generation per unit time at the end of the time increment caused by mechanical working of the material
         dq_dT, &                                                        ! variation of volumetric heat generation per unit time w.r.t. temperature
         pnewdt                                                          ! ratio of suggested new time increment to the time increment being used (delta_t_next/delta_t_current)

     real(pReal), dimension(n_total),         intent(inout) :: &
         sigma                                                           ! stress tensor (voigt notation, vector form) at the beginning of the time increment, needs to be updated at the end of the time increment [MPa]

     real(pReal), dimension(n_sv),            intent(inout) :: &
         sv                                                              ! (solution-dependent) state variables array at the beginning of the time increment, needs to be updated at the end of the time increment

     real(pReal), dimension(n_total),         intent(out) :: &
         dsigma_dT, &                                                    ! derivative of stress tensor w.r.t. temperature (voigt notation, vector form) [MPa/T]
         dq_deps                                                         ! derivative of heat generation w.r.t. strain tensor (voigt notation, vector form)

     real(pReal), dimension(n_total,n_total), intent(out) :: &
         C                                                               ! consistent/algorithmic tangent stiffness tensor or dsigma_deps (constitutive/material Jacobian) (voigt notation, vector form) [MPa]

!-----------------------------------------------------------------------
!    declaration of subroutines' local variables
     integer(pInt) :: &
         i_material, &                                                   ! index number associated with each material
         i, &                                                            ! loop counter variable i
         j, &                                                            ! loop counter variable j
         k_NR                                                            ! newton-raphson (NR) iterative loop index

     real(pReal) :: &
         rho_hat_cm_n, &                                                 ! normalized cell mobile   dislocation density at the beginning of the time increment [-]
         rho_hat_ci_n, &                                                 ! normalized cell immobile dislocation density at the beginning of the time increment [-]
         rho_hat_wi_n, &                                                 ! normalized wall immobile dislocation density at the beginning of the time increment [-]
         rho_hat_cm, &                                                   ! normalized cell mobile   dislocation density at the end       of the time increment [-]
         rho_hat_ci, &                                                   ! normalized cell immobile dislocation density at the end       of the time increment [-]
         rho_hat_wi, &                                                   ! normalized wall immobile dislocation density at the end       of the time increment [-]
         eps_bar_p, &                                                    ! equivalent plastic strain [-]
         eps_bar_dot_p, &                                                ! equivalent plastic strain rate [s^-1]
         T_hat_n, &                                                      ! normalized absolute temperature at the beginning of the time increment [-]
         nu, &                                                           ! poisson's ratio [-]
         G, &                                                            ! shear modulus [MPa]
         K, &                                                            ! bulk modulus [MPa]
         sigma_bar_trial, &                                              ! equivalent trial stress [MPa]
         delta_eps_bar_p, &                                              ! equivalent plastic strain increment [-]
         sigma_p, &                                                      ! plastic stress [MPa]
         H_p, &                                                          ! plastic tangent modulus (derivative of plastic stress w.r.t. equivalent plastic strain increment) [MPa]
         beta, &                                                         ! dissipation factor [-]
         sigma_y, &                                                      ! yield stress [MPa]
         sigma_y_trial, &                                                ! trial yield stress [MPa]
         sigma_v0, &                                                     ! Reference viscous stress at current temperature [MPa]
         m_v, &                                                          ! strain rate sensitivity associated with viscous stress at current temperature [-]
         sigma_v, &                                                      ! viscous stress [MPa]
         R, &                                                            ! residual function to be solved using the newton-raphson method [-]
         G_eff, &                                                        ! effective modulus mu [MPa]
         H_vp, &                                                         ! viscoplastic tangent modulus (derivative of yield stress w.r.t. equivalent plastic strain increment) [MPa]
         H_eff                                                           ! effective viscoplastic tangent modulus [MPa]

     real(pReal), dimension(n_total) :: &
         eps_e, &                                                        ! elastic strains tensor (vector form) [-]
         eps_p, &                                                        ! plastic strains tensor (vector form) [-]
         sigma_h, &                                                      ! hydrostatic stress tensor (vector form) [MPa]
         sigma_d_trial, &                                                ! trial deviatoric stress tensor (vector form) [MPa]
         delta_eps_p, &                                                  ! plastic strain increment tensor [-]
         N, &                                                            ! flow direction tensor (vector form) [-]
         sigma_d                                                         ! deviatoric stress tensor (vector form) [MPa]

!-----------------------------------------------------------------------
!    assigning values to (initializing) global variables associated with material parameters
     if (.not. material_parameters_assigned) then
         call MutexInit(1)
         call MutexLock(1)
             if (.not. material_parameters_assigned) then
                 call assign_material_parameters()
             end if
         call MutexUnlock(1)
     end if

!-----------------------------------------------------------------------
!    reading the material index number associated with the current material point
     i_material = int(props(1))

!-----------------------------------------------------------------------
!    initializing/retrieving state variables
!
!    state variables array:
!    2D case (*DEPVAR = n_sv = 13)
!    sv(1)     : cell mobile   dislocation density [m^-2]
!    sv(2)     : cell immobile dislocation density [m^-2]
!    sv(3)     : wall immobile dislocation density [m^-2]
!    sv(4:7)   : elastic strain tensor (vector form) [-]
!    sv(8:11)  : plastic strain tensor (vector form) [-]
!    sv(12)    : equivalent plastic strain [-]
!    sv(13)    : equivalent plastic strain rate [s^-1]
!
!    3D case (*DEPVAR = n_sv = 17)
!    sv(1)     : cell mobile   dislocation density [m^-2]
!    sv(2)     : cell immobile dislocation density [m^-2]
!    sv(3)     : wall immobile dislocation density [m^-2]
!    sv(4:9)   : elastic strain tensor (vector form) [-]
!    sv(10:15) : plastic strain tensor (vector form) [-]
!    sv(16)    : equivalent plastic strain [-]
!    sv(17)    : equivalent plastic strain rate [s^-1]

     if (t(2) <= delta_t) then
!        initializing state variables for undeformed state
         rho_hat_cm_n  = rho_hat_cm0(i_material)
         rho_hat_ci_n  = rho_hat_ci0(i_material)
         rho_hat_wi_n  = rho_hat_wi0(i_material)
         eps_e         = 0.0e0
         eps_p         = 0.0e0
         eps_bar_p     = 0.0e0
         eps_bar_dot_p = 0.0e0
     else
!        retrieving state variables at the beginning of the time increment (end of the last time step) [-]
         rho_hat_cm_n  = sv(1) / rho_0(i_material)   ! {Eq. 84(b,c)}
         rho_hat_ci_n  = sv(2) / rho_0(i_material)   ! {Eq. 84(b,c)}
         rho_hat_wi_n  = sv(3) / rho_0(i_material)   ! {Eq. 84(b,c)}
         call rotsig    (sv(4          ), delta_R, eps_e, 2, n_dir, n_shr)
         call rotsig    (sv(n_total + 4), delta_R, eps_p, 2, n_dir, n_shr)
         eps_bar_p     = sv(2 * n_total + 4)
         eps_bar_dot_p = sv(2 * n_total + 5)
     end if

!-----------------------------------------------------------------------
! processing
!-----------------------------------------------------------------------
!    calculation of normalized relative temperatures   {Eq. 133(b)}
     T_hat_n = (T_n - T_abs0) / (T_0(i_material) - T_abs0)

!    calculation of elastic stiffness moduli
     nu = nu_0(i_material) * (1.0e0 + r_nu(i_material) * (T_hat_n - 1.0e0)**s_nu(i_material))   ! {Eq. 161}
     G  = G_0 (i_material) * (1.0e0 + r_G (i_material) * (T_hat_n - 1.0e0)**s_G (i_material))   ! {Eq. 160}
     K  = 2.0e0 / 3.0e0 * (1.0e0 + nu) / (1.0e0 - 2.0e0 * nu) * G                               ! {Eq. Box 1.1(b)}

!    constructing elastic stiffness tensor   {Eq. Box 1.1}
     C = 0.0e0
     do i = 1, n_dir
         do j = 1, n_dir
             C(i,j) = K - 2.0e0 / 3.0e0 * G
         end do
     end do
     do i = 1, n_total
         C(i,i) = C(i,i) + 2.0e0 * G
     end do

!    calculation of trial stress tensor   {Eq. Box 1.2, Eq. Box 1.5}
     do i = 1, n_total
         do j = 1, n_total
             sigma(i) = sigma(i) + C(i,j) * delta_eps(j)
         end do
     end do

!    calculation of hydrostatic stress tensor   {Eq. Box 1.3(b,c)}
     sigma_h = 0.0e0
     sigma_h(1:n_dir) = sum(sigma(1:n_dir)) / 3.0e0

!    calculation of trial deviatoric stress tensor   {Eq. Box 1.3}
     sigma_d_trial = sigma - sigma_h

!    calculation of equivalent trial stress   {Eq. Box 1.4(b)}
     sigma_bar_trial = 0.0e0
     do i = 1, n_dir
         sigma_bar_trial = sigma_bar_trial + sigma_d_trial(i)**2
     end do
     do i = n_dir + 1, n_total
         sigma_bar_trial = sigma_bar_trial + 2.0e0 * sigma_d_trial(i)**2
     end do
     sigma_bar_trial = sqrt(3.0e0 / 2.0e0 * sigma_bar_trial)

!    initializing (trial) plastic strain increments   {Eq. Box 1.5(d)}
     delta_eps_bar_p = 0.0e0
     delta_eps_p     = 0.0e0

!    calculation of trial plastic stress   {Eq. 154}
     call plasticity(i_material, delta_t, T_hat_n, delta_eps_bar_p, rho_hat_cm_n, rho_hat_ci_n, rho_hat_wi_n, rho_hat_cm, rho_hat_ci, rho_hat_wi, sigma_p, H_p, beta)
     sigma_y       = sigma_p
     sigma_y_trial = sigma_y

!    calculating temperature-dependent reference viscous stress and strain rate sensitivity parameter
     sigma_v0 = sigma_v00(i_material) * (1.0e0 + r_v (i_material) * (T_hat_n - 1.0e0)**s_v (i_material))   ! {Eq. 134(b)}
     m_v      = m_v0     (i_material) * (1.0e0 + r_mv(i_material) * (T_hat_n - 1.0e0)**s_mv(i_material))   ! {Eq. 135}

!    checking for plastic yielding: plastic yielding occurs if equivalent trial stress is greater than trial yield stress
     if(sigma_bar_trial > sigma_y_trial) then   ! {Eq. 128}
!        initializing variables for multivariate return mapping
         delta_eps_bar_p = delta_t * eps_bar_dot_p   ! {Eq. 155}
         call plasticity(i_material, delta_t, T_hat_n, delta_eps_bar_p, rho_hat_cm_n, rho_hat_ci_n, rho_hat_wi_n, rho_hat_cm, rho_hat_ci, rho_hat_wi, sigma_p, H_p, beta)

!        return mapping loop: solving for equivalent plastic strain increment using newton-raphson iterative method
         do k_NR = 1, k_NR_max
             sigma_y = sigma_bar_trial - 3.0e0 * G * delta_eps_bar_p                           ! {Eq. 152(c)}
             sigma_v = sigma_y - sigma_p                                                       ! {Eq. 152(b)}
             eps_bar_dot_p = eps_bar_dot_0(i_material) * (sigma_v / sigma_v0)**(1.0e0 / m_v)   ! {Eq. 152}
             R = delta_eps_bar_p - delta_t * eps_bar_dot_p                                     ! {Eq. 153}
             if (abs(R) < chi * delta_t * xi_mean * eps_bar_dot_0(i_material)) exit            ! {Eq. 156}
             delta_eps_bar_p = delta_eps_bar_p - R / (1.0e0 + (3.0e0 * G + H_p) / (m_v * sigma_v) * delta_t * eps_bar_dot_p)   ! {Eqs. 147, 157}
             call plasticity(i_material, delta_t, T_hat_n, delta_eps_bar_p, rho_hat_cm_n, rho_hat_ci_n, rho_hat_wi_n, rho_hat_cm, rho_hat_ci, rho_hat_wi, sigma_p, H_p, beta)
         end do

!        checking for divergence: reducing time step size and writing warning message
         if (abs(R) > chi * delta_t * xi_mean * eps_bar_dot_0(i_material)) then
             pnewdt = 0.5e0
             write(6,"('warning: return mapping algorithm for element ', i5, ', integration point ', i2, ' did not converge after ', i2,' iterations')") i_el, i_ip, k_NR
             return
         end if

!        calculation of flow direction tensor   {Eq. Box 1.6(c)}
         N = 3.0e0 / 2.0e0 * sigma_d_trial / sigma_bar_trial

!        updating stress tensor
         sigma_d = 2.0e0 / 3.0e0 * sigma_y * N   ! {Eq. Box 1.6(b)}
         sigma   = sigma_d + sigma_h             ! {Eq. Box 1.6}

!        calculation of effective moduli
         G_eff = sigma_y / sigma_bar_trial * G                              ! {Eq. Box 1.7(b)}
         H_vp  = m_v / delta_eps_bar_p * sigma_v + H_p                      ! {Eqs. 149, 150}
         H_eff = 4.0e0 / 3.0e0 * (G / (1.0e0 + 3.0e0 * G / H_vp) - G_eff)   ! {Eq. Box 1.7}

!        updating consistent tangent stiffness tensor   {Eq. Box 1.8}
         C = 0.0e0
         do i = 1, n_dir
             do j = 1, n_dir
                 C(i,j) = K - 2.0e0 / 3.0e0 * G_eff
             end do
         end do
         do i = 1, n_total
             C(i,i) = C(i,i) + 2.0e0 * G_eff
         end do
         do i = 1, n_total
             do j = 1, n_total
                 C(i,j) = C(i,j) + H_eff * N(i) * N(j)
             end do
         end do

!        updating plastic strain increment tensor   {Eq. Box 1.6(c)}
         delta_eps_p = delta_eps_bar_p * N
     end if

!    calculation of plastic dissipation work and plastic heat generation   {Eqs. Box 1.9, 105}
     delta_w_p = delta_eps_bar_p * sigma_y
     q_dot_p   = beta * delta_w_p / delta_t

!    calculation of strains   {Eq. 105}
     eps_e = eps_e + delta_eps - delta_eps_p
     eps_p = eps_p + delta_eps_p
     eps_bar_p     = eps_bar_p + delta_eps_bar_p

!-----------------------------------------------------------------------
! post-processing
!-----------------------------------------------------------------------
!    updating state variables
     sv(1) = rho_hat_cm * rho_0(i_material)   ! {Eq. 84(b,c)}
     sv(2) = rho_hat_ci * rho_0(i_material)   ! {Eq. 84(b,c)}
     sv(3) = rho_hat_wi * rho_0(i_material)   ! {Eq. 84(b,c)}

     sv(4 : n_total + 3)               = eps_e
     sv(n_total + 4 : 2 * n_total + 3) = eps_p
     sv(2 * n_total + 4)               = eps_bar_p
     sv(2 * n_total + 5)               = eps_bar_dot_p

     return
  end subroutine umat