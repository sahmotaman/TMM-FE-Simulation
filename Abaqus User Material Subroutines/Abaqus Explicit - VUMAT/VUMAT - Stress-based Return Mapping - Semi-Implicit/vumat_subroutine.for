!=======================================================================
! User Material Subroutine for Isotropic Elasto-Visco-Plasticity
! CDD-based Material Model
! Semi-Implicit Constitutive Integration with Stress-based Return Mapping
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
! including other modules and subroutines
#include 'controls_module.f90'
#include 'strings_module.f90'
#include 'material_properties_module.f90'
#include 'viscoplasticity_subroutine.f90'


!***********************************************************************
! subroutine header
  subroutine vumat(n_block,     n_dir,   n_shr,         n_sv, &
                   n_fv,        n_props, anneal_status, t_step, &
                   t_total,     delta_t, cmname,        coords, &
                   l_element,   props,   rho_mass,      delta_eps, &
                   delta_W_rel, T_n,     U_n,           F_n, &
                   fv_n,        sigma_n, sv_n,          w_n, &
                   q_p_n,       T,       U,             F, &
                   fv,          sigma,   sv,            w, &
                   q_p)

!-----------------------------------------------------------------------
! pre-processing
!-----------------------------------------------------------------------
!    use of global variables
     use controls
     use material_properties

!-----------------------------------------------------------------------
!    declaration of subroutine's variables
     implicit none

     character(len = 64),                                    intent(in) :: &
         cmname                                                          ! uses-specified material name, left justified

     integer(pInt),                                          intent(in) :: &
         n_dir, &                                                        ! number of independent direct components of strain and stress tensors
         n_shr, &                                                        ! number of independent shear  components of strain and stress tensors
         n_sv, &                                                         ! number of user-defined state variables
         n_fv, &                                                         ! number of user-defined field variables
         n_props, &                                                      ! number of user-defined material parameters
         anneal_status                                                   ! status flag indicating whether the routine is being called during an annealing process (0:normal mechanics increment, 1:annealing increment)

     integer(pInt), dimension(*),                            intent(in) :: &
         n_block                                                         ! number of 1: material points in current call, 2: (total) material points, 3: layers, 4: section points, 5: elements

     real(pReal),                                            intent(in) :: &
         t_step, &                                                       ! elapsed simulation step  time [s]
         t_total, &                                                      ! elapsed simulation total time [s]
         delta_t                                                         ! time increment [s]

     real(pReal),   dimension(n_props),                      intent(in) :: &
         props                                                           ! material constants array

     real(pReal),   dimension(n_block(1)),                   intent(in) :: &
         l_element, &                                                    ! characteristic element length [mm]
         rho_mass, &                                                     ! current density at the material points in the midstep configuration (not affected by mass scaling) [kg/mm^3]
         w_n, &                                                          ! specific internal energy (total elasto-plastic work) (per unit mass) at each material point at the beginning of the time increment [uJ.kg^-1]
         q_p_n, &                                                        ! specific dissipated inelastic energy (per unit mass) at each material point at the beginning of the time increment [uJ.kg^-1]
         T_n, &                                                          ! temperature at each material point at the beginning of the time increment [C]
         T                                                               ! temperature at each material point at the end       of the time increment [C]

     real(pReal),   dimension(n_block(1),*),                 intent(in) :: &
         coords                                                          ! material point coordinates [mm]

     real(pReal),   dimension(n_block(1),n_sv),              intent(in) :: &
         sv_n                                                            ! (solution-dependent) state variables array at each material point at the beginning of the increment

     real(pReal),   dimension(n_block(1),n_fv),              intent(in) :: &
         fv_n, &                                                         ! array of user-defined field variables at each material point at the beginning of the increment
         fv                                                              ! array of user-defined field variables at each material point at the end       of the increment

     real(pReal),   dimension(n_block(1),n_shr),             intent(in) :: &
         delta_W_rel                                                     ! incremental relative rotation/spin tensor (voigt notation, vector form) at each material point defined in the corotational system [-]

     real(pReal),   dimension(n_block(1),n_dir + n_shr),     intent(in) :: &
         delta_eps, &                                                    ! (total) strain increment tensor (voigt notation, vector form) at each material point [-]
         U_n, &                                                          ! stretch tensor (voigt notation, vector form) at each material point at the beginning of the time increment [-]
         U, &                                                            ! stretch tensor (voigt notation, vector form) at each material point at the end       of the time increment [-]
         sigma_n                                                         ! stress  tensor (voigt notation, vector form) at each material point at the beginning of the time increment [MPa]

     real(pReal),   dimension(n_block(1),n_dir + 2 * n_shr), intent(in) :: &
         F_n, &                                                          ! deformation gradient tensor (voigt notation, vector form) at each material point at the beginning of the time increment [-]
         F                                                               ! deformation gradient tensor (voigt notation, vector form) at each material point at the end       of the time increment [-]

     real(pReal),   dimension(n_block(1)),                   intent(out) :: &
         w, &                                                            ! specific internal energy (total elasto-plastic work) (per unit mass) at each material point at the end of the time increment [uJ.kg^-1]
         q_p                                                             ! specific dissipated inelastic energy (per unit mass) at each material point at the end of the time increment [uJ.kg^-1]

     real(pReal),   dimension(n_block(1),n_sv),              intent(out) :: &
         sv                                                              ! (solution-dependent) state variables array at each material point at the end of the increment

     real(pReal),   dimension(n_block(1),n_dir + n_shr),     intent(out) :: &
         sigma                                                           ! stress  tensor (voigt notation, vector form) at each material point at the end of the time increment [MPa]

!-----------------------------------------------------------------------
!    declaration of subroutines' local variables
     logical :: &
         trial_flag                                                      ! flag to determine whether the call to the yield function is for calculation of trial stress or not

     integer(pInt) :: &
         i_material, &                                                   ! index number associated with each material
         i_mp, &                                                         ! material point index number
         k_NR                                                            ! newton-raphson (NR) iterative loop index

     real(pReal) :: &
         rho_hat_cm_n, &                                                 ! normalized cell mobile   dislocation density at the beginning of the time increment [-]
         rho_hat_ci_n, &                                                 ! normalized cell immobile dislocation density at the beginning of the time increment [-]
         rho_hat_wi_n, &                                                 ! normalized wall immobile dislocation density at the beginning of the time increment [-]
         rho_hat_cm, &                                                   ! normalized cell mobile   dislocation density at the end       of the time increment [-]
         rho_hat_ci, &                                                   ! normalized cell immobile dislocation density at the end       of the time increment [-]
         rho_hat_wi, &                                                   ! normalized wall immobile dislocation density at the end       of the time increment [-]
         eps_bar_p_n, &                                                  ! accumulated equivalent plastic strain at the beginning of the time increment [-]
         eps_bar_p, &                                                    ! accumulated equivalent plastic strain at the end       of the time increment [-]
         eps_bar_dot_p_n, &                                              ! equivalent plastic strain rate at the beginning of the time increment [s^-1]
         eps_bar_dot_p, &                                                ! equivalent plastic strain rate at the end       of the time increment [s^-1]
         T_hat_n, &                                                      ! normalized absolute temperature at the beginning of the time increment [-] [-]
         nu, &                                                           ! poisson's ratio [-]
         G, &                                                            ! shear modulus [MPa]
         K, &                                                            ! bulk modulus [MPa]
         delta_eps_v, &                                                  ! volumetric strain increment [-]
         sigma_h, &                                                      ! hydrostatic stress [MPa]
         sigma_bar_trial, &                                              ! equivalent trial stress [MPa]
         eps_dot_p_min, &                                                ! minimum plastic strain rate [s^-1]
         eps_bar_dot_p_corr, &                                           ! corrected equivalent plastic strain rate [s^-1]
         delta_eps_bar_p, &                                              ! equivalent plastic strain increment [-]
         sigma_y, &                                                      ! yield stress [MPa]
         sigma_y_trial, &                                                ! trial yield stress [MPa]
         H_vp, &                                                         ! viscoplastic tangent modulus (derivative of yield stress w.r.t. equivalent plastic strain increment) [MPa]
         beta, &                                                         ! dissipation factor [-]
         R, &                                                            ! residual/yield function to be solved using the newton-raphson method [MPa]
         delta_w                                                         ! incremental specific internal energy (incremental total elasto-plastic work) (per unit mass) at each material point [uJ.kg^-1]

     real(pReal), dimension(n_dir + n_shr) :: &
         sigma_d_trial                                                   ! trial deviatoric stress tensor (vector form) [MPa]

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

!    loop for material points
     do i_mp = 1, n_block(1)
!        state variables array:
!        *DEPVAR = n_sv = 5
!        sv(1) : cell mobile   dislocation density [m^-2]
!        sv(2) : cell immobile dislocation density [m^-2]
!        sv(3) : wall immobile dislocation density [m^-2]
!        sv(4) : accumulated equivalent plastic strain [-]
!        sv(5) : equivalent plastic strain rate [-]

         if (t_total <= delta_t) then
!            initializing state variables for undeformed state
             rho_hat_cm_n    = rho_hat_cm0(i_material)
             rho_hat_ci_n    = rho_hat_ci0(i_material)
             rho_hat_wi_n    = rho_hat_wi0(i_material)
             eps_bar_p_n     = 0.0e0
             eps_bar_dot_p_n = 0.0e0
         else
!            retrieving state variables at the beginning of the time increment (end of the last time step) [-]
             rho_hat_cm_n    = sv_n(i_mp,1) / rho_0(i_material)   ! {Eq. 59(b,c)}
             rho_hat_ci_n    = sv_n(i_mp,2) / rho_0(i_material)   ! {Eq. 59(b,c)}
             rho_hat_wi_n    = sv_n(i_mp,3) / rho_0(i_material)   ! {Eq. 59(b,c)}
             eps_bar_p_n     = sv_n(i_mp,4)
             eps_bar_dot_p_n = sv_n(i_mp,5)
         end if

!-----------------------------------------------------------------------
! processing
!-----------------------------------------------------------------------
!        calculation of normalized relative temperatures   {Eq. 93(d)}
         T_hat_n = (T_n(i_mp) - T_abs0) / (T_0(i_material) - T_abs0)

!        calculation of elastic stiffness moduli
         nu = nu_0(i_material) * (1.0e0 + r_nu(i_material) * (T_hat_n - 1.0e0)**s_nu(i_material))   ! {Eq. 118}
         G  = G_0 (i_material) * (1.0e0 + r_G (i_material) * (T_hat_n - 1.0e0)**s_G (i_material))   ! {Eq. 117}
         K  = 2.0e0 / 3.0e0 * (1.0e0 + nu) / (1.0e0 - 2.0e0 * nu) * G                               ! {Eq. Box 1.1(b)}

!        calculation of volumetric strain increment   {Eq. Box 1.2(b)}
         delta_eps_v = delta_eps(i_mp,1) + delta_eps(i_mp,2) + delta_eps(i_mp,3)

!        calculation of trial stress tensor   {Eq. Box 1.2, Eq. Box 1.5}
         sigma(i_mp,1) = sigma_n(i_mp,1) + 2.0e0 * G * delta_eps(i_mp,1) + (K - 2.0e0 / 3.0e0 * G) * delta_eps_v
         sigma(i_mp,2) = sigma_n(i_mp,2) + 2.0e0 * G * delta_eps(i_mp,2) + (K - 2.0e0 / 3.0e0 * G) * delta_eps_v
         sigma(i_mp,3) = sigma_n(i_mp,3) + 2.0e0 * G * delta_eps(i_mp,3) + (K - 2.0e0 / 3.0e0 * G) * delta_eps_v
         sigma(i_mp,4) = sigma_n(i_mp,4) + 2.0e0 * G * delta_eps(i_mp,4)
         if (n_shr > 1) then
             sigma(i_mp,5) = sigma_n(i_mp,5) + 2.0e0 * G * delta_eps(i_mp,5)
             sigma(i_mp,6) = sigma_n(i_mp,6) + 2.0e0 * G * delta_eps(i_mp,6)
         end if

!        calculation of hydrostatic stress   {Eq. Box 1.3(c)}
         sigma_h = (sigma(i_mp,1) + sigma(i_mp,2) + sigma(i_mp,3)) / 3.0e0

!        calculation of trial deviatoric stress tensor   {Eq. Box 1.3}
         sigma_d_trial(1) = sigma(i_mp,1) - sigma_h
         sigma_d_trial(2) = sigma(i_mp,2) - sigma_h
         sigma_d_trial(3) = sigma(i_mp,3) - sigma_h
         sigma_d_trial(4) = sigma(i_mp,4)
         if (n_shr > 1) then
             sigma_d_trial(5) = sigma(i_mp,5)
             sigma_d_trial(6) = sigma(i_mp,6)
         end if

!        calculation of equivalent trial stress   {Eq. Box 1.4(b)}
         if (n_shr == 1) then
             sigma_bar_trial = sqrt(3.0e0 / 2.0e0 * (sigma_d_trial(1)**2 + sigma_d_trial(2)**2 + sigma_d_trial(3)**2 + 2.0e0 * sigma_d_trial(4)**2))
         else
             sigma_bar_trial = sqrt(3.0e0 / 2.0e0 * (sigma_d_trial(1)**2 + sigma_d_trial(2)**2 + sigma_d_trial(3)**2 + 2.0e0 * (sigma_d_trial(4)**2 + sigma_d_trial(5)**2 + sigma_d_trial(6)**2)))
         end if

!        calculation of corrected equivalent strain rate   {Eq. 99}
         eps_dot_p_min = xi_min * eps_bar_dot_0(i_material)
         if (eps_bar_dot_p_n > eps_dot_p_min) then
             eps_bar_dot_p_corr = eps_bar_dot_p_n
         else
             eps_bar_dot_p_corr = eps_dot_p_min
         end if

!        initializing (trial) equivalent plastic strain increment   {Eq. Box 1.5(d)}
         delta_eps_bar_p = 0.0e0

!        calculation of trial plastic stress   {Eq. 101}
         trial_flag = .true.
         call viscoplasticity(trial_flag, i_material, delta_t, T_hat_n, eps_bar_dot_p_corr, delta_eps_bar_p, rho_hat_cm_n, rho_hat_ci_n, rho_hat_wi_n, rho_hat_cm, rho_hat_ci, rho_hat_wi, sigma_y, H_vp, beta)
         sigma_y_trial = sigma_y

!        checking for plastic yielding: plastic yielding occurs if equivalent trial stress is greater than trial plastic stress
         if(sigma_bar_trial > sigma_y_trial) then   ! {Eq. 128}
!            initializing variables for return mapping
             trial_flag = .false.
             delta_eps_bar_p = eps_bar_dot_p_corr * delta_t   ! {Eq. 145}

!            return mapping loop: solving for equivalent plastic strain increment using newton-raphson iterative method
             do k_NR = 1, k_NR_max
                 call viscoplasticity(trial_flag, i_material, delta_t, T_hat_n, eps_bar_dot_p_corr, delta_eps_bar_p, rho_hat_cm_n, rho_hat_ci_n, rho_hat_wi_n, rho_hat_cm, rho_hat_ci, rho_hat_wi, sigma_y, H_vp, beta)   ! {Eq. Box 1.7(b)}
                 R = sigma_bar_trial - sigma_y - 3.0e0 * G * delta_eps_bar_p   ! {Eq. 96}
                 if (abs(R) < chi * sigma_y_trial) exit                        ! {Eq. 103}
                 delta_eps_bar_p = delta_eps_bar_p + R / (H_vp + 3.0e0 * G)    ! {Eq. 105}
             end do

!            checking for divergence: reducing time step size and writing warning message
             if (abs(R) > chi * sigma_y_trial) then
                 write(6,"('warning: return mapping algorithm did not converge after ', i2,' iterations')") k_NR
                 return
             end if

!            updating stress tensor   {Eq. Box 1.6}
             sigma(i_mp,1) = sigma_y / sigma_bar_trial * sigma_d_trial(1) + sigma_h
             sigma(i_mp,2) = sigma_y / sigma_bar_trial * sigma_d_trial(2) + sigma_h
             sigma(i_mp,3) = sigma_y / sigma_bar_trial * sigma_d_trial(3) + sigma_h
             sigma(i_mp,4) = sigma_y / sigma_bar_trial * sigma_d_trial(4)
             if (n_shr > 1) then
                 sigma(i_mp,5) = sigma_y / sigma_bar_trial * sigma_d_trial(5)
                 sigma(i_mp,6) = sigma_y / sigma_bar_trial * sigma_d_trial(6)
             end if
         end if

!        calculation of total specific elasto-plastic work
         if (n_shr == 1) then
             delta_w = ((sigma_n(i_mp,1) + sigma(i_mp,1)) * delta_eps(i_mp,1) &
                     +  (sigma_n(i_mp,2) + sigma(i_mp,2)) * delta_eps(i_mp,2) &
                     +  (sigma_n(i_mp,3) + sigma(i_mp,3)) * delta_eps(i_mp,3) &
                     +  (sigma_n(i_mp,4) + sigma(i_mp,4)) * delta_eps(i_mp,4)) / (2.0e0 * rho_mass(i_mp))
         else
             delta_w = ((sigma_n(i_mp,1) + sigma(i_mp,1)) * delta_eps(i_mp,1) &
                     +  (sigma_n(i_mp,2) + sigma(i_mp,2)) * delta_eps(i_mp,2) &
                     +  (sigma_n(i_mp,3) + sigma(i_mp,3)) * delta_eps(i_mp,3) &
                     +  (sigma_n(i_mp,4) + sigma(i_mp,4)) * delta_eps(i_mp,4) &
                     +  (sigma_n(i_mp,5) + sigma(i_mp,5)) * delta_eps(i_mp,5) &
                     +  (sigma_n(i_mp,6) + sigma(i_mp,6)) * delta_eps(i_mp,6)) / (2.0e0 * rho_mass(i_mp))
         end if
         w(i_mp) = w_n(i_mp) + delta_w   ! {Eq. 68}

!        calculation of plastic dissipation work and plastic heat generation   {Eqs. Box 1.9, 68}
         q_p(i_mp) = q_p_n(i_mp) + beta * sigma_y * delta_eps_bar_p / rho_mass(i_mp)

!        calculation of equivalent plastic strain rate and equivalent plastic strain
         eps_bar_dot_p = delta_eps_bar_p / delta_t       ! {Eq. Box 1.9}
         eps_bar_p     = eps_bar_p_n + delta_eps_bar_p   ! {Eq. 68}

!-----------------------------------------------------------------------
! post-processing
!-----------------------------------------------------------------------
!        updating state variables
         sv(i_mp,1) = rho_hat_cm * rho_0(i_material)   ! {Eq. 59(b,c)}
         sv(i_mp,2) = rho_hat_ci * rho_0(i_material)   ! {Eq. 59(b,c)}
         sv(i_mp,3) = rho_hat_wi * rho_0(i_material)   ! {Eq. 59(b,c)}
         sv(i_mp,4) = eps_bar_p
         sv(i_mp,5) = eps_bar_dot_p
     end do

     return
  end subroutine vumat