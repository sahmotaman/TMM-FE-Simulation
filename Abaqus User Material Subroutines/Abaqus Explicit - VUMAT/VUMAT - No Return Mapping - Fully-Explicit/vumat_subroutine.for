!=======================================================================
! User Material Subroutine for Isotropic Elasto-Visco-Plasticity
! CDD-based Material Model
! Fully-Explicit Constitutive Integration
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
         i_mp                                                            ! material point index number

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
         sigma_h_n, &                                                    ! hydrostatic stress at the beginning of the time increment [MPa]
         sigma_bar_n, &                                                  ! equivalent  stress at the beginning of the time increment [MPa]
         T_hat_n, &                                                      ! normalized absolute temperature at the beginning of the time increment [-]
         G_alpha_c, &                                                    ! combination of shear modulus and interaction coefficient associated with cell immobile dislocations at current temperature [-]
         G_alpha_w, &                                                    ! combination of shear modulus and interaction coefficient associated with wall immobile dislocations at current temperature [-]
         sigma_p_n, &                                                    ! plastic stress at the beginning of the time increment [MPa]
         sigma_v_n, &                                                    ! viscous stress at the beginning of the time increment [MPa]
         m_v, &                                                          ! strain rate sensitivity associated with viscous stress at current temperature [-]
         sigma_v0, &                                                     ! viscous stress at reference strain rate [MPa]
         delta_eps_bar_p, &                                              ! equivalent plastic strain increment [-]
         nu, &                                                           ! poisson's ratio [-]
         G, &                                                            ! shear modulus [MPa]
         K, &                                                            ! bulk modulus [MPa]
         delta_eps_v, &                                                  ! volumetric strain increment [-]
         m_an_cm, &                                                      ! strain rate sensitivity parameter associated with dynamic annihilation   of cell mobile   dislocations [-]
         m_an_ci, &                                                      ! strain rate sensitivity parameter associated with dynamic annihilation   of cell immobile dislocations [-]
         m_an_wi, &                                                      ! strain rate sensitivity parameter associated with dynamic annihilation   of wall immobile dislocations [-]
         m_tr_cm, &                                                      ! strain rate sensitivity parameter associated with dynamic trapping       of cell mobile   dislocations [-]
         m_nc_wi, &                                                      ! strain rate sensitivity parameter associated with dynamic nucleation     of wall immobile dislocations [-]
         m_rm_ci, &                                                      ! strain rate sensitivity parameter associated with dynamic remobilization of cell immobile dislocations [-]
         m_rm_wi, &                                                      ! strain rate sensitivity parameter associated with dynamic remobilization of wall immobile dislocations [-]
         c_gn_cm, &                                                      ! material parameter associated with dynamic generation     of cell mobile   dislocations at current temperature [-]
         c_an_cm, &                                                      ! material parameter associated with dynamic annihilation   of cell mobile   dislocations at current temperature [-]
         c_an_ci, &                                                      ! material parameter associated with dynamic annihilation   of cell immobile dislocations at current temperature [-]
         c_an_wi, &                                                      ! material parameter associated with dynamic annihilation   of wall immobile dislocations at current temperature [-]
         c_ac_ci, &                                                      ! material parameter associated with dynamic accumulation   of cell immobile dislocations at current temperature [-]
         c_ac_wi, &                                                      ! material parameter associated with dynamic accumulation   of wall immobile dislocations at current temperature [-]
         c_tr_cm, &                                                      ! material parameter associated with dynamic trapping       of cell mobile   dislocations at current temperature [-]
         c_nc_wi, &                                                      ! material parameter associated with dynamic nucleation     of cell immobile dislocations at current temperature [-]
         c_rm_ci, &                                                      ! material parameter associated with dynamic remobilization of cell immobile dislocations at current temperature [-]
         c_rm_wi, &                                                      ! material parameter associated with dynamic remobilization of wall immobile dislocations at current temperature [-]
         del_rho_hat_gn_cm_n, &                                          ! normalized dynamic generation     rate of cell mobile    dislocations w.r.t. equivalent plastic strain at the beginning of the time increment [-]
         del_rho_hat_an_cm_n, &                                          ! normalized dynamic annihilation   rate of cell mobile    dislocations w.r.t. equivalent plastic strain at the beginning of the time increment [-]
         del_rho_hat_an_ci_n, &                                          ! normalized dynamic annihilation   rate of cell immmobile dislocations w.r.t. equivalent plastic strain at the beginning of the time increment [-]
         del_rho_hat_an_wi_n, &                                          ! normalized dynamic annihilation   rate of wall immmobile dislocations w.r.t. equivalent plastic strain at the beginning of the time increment [-]
         del_rho_hat_ac_ci_n, &                                          ! normalized dynamic accumulation   rate of cell immmobile dislocations w.r.t. equivalent plastic strain at the beginning of the time increment [-]
         del_rho_hat_ac_wi_n, &                                          ! normalized dynamic accumulation   rate of wall immmobile dislocations w.r.t. equivalent plastic strain at the beginning of the time increment [-]
         del_rho_hat_tr_cm_n, &                                          ! normalized dynamic trapping       rate of cell mobile    dislocations w.r.t. equivalent plastic strain at the beginning of the time increment [-]
         del_rho_hat_nc_wi_n, &                                          ! normalized dynamic nucleation     rate of wall immmobile dislocations w.r.t. equivalent plastic strain at the beginning of the time increment [-]
         del_rho_hat_rm_ci_n, &                                          ! normalized dynamic remobilization rate of cell immmobile dislocations w.r.t. equivalent plastic strain at the beginning of the time increment [-]
         del_rho_hat_rm_wi_n, &                                          ! normalized dynamic remobilization rate of wall immmobile dislocations w.r.t. equivalent plastic strain at the beginning of the time increment [-]
         del_rho_hat_cm_n, &                                             ! derivative of cell mobile   dislocation density                       w.r.t. equivalent plastic strain at the beginning of the time increment [-]
         del_rho_hat_ci_n, &                                             ! derivative of cell immobile dislocation density                       w.r.t. equivalent plastic strain at the beginning of the time increment [-]
         del_rho_hat_wi_n, &                                             ! derivative of wall immobile dislocation density                       w.r.t. equivalent plastic strain at the beginning of the time increment [-]
         delta_w, &                                                      ! incremental specific internal energy (incremental total elasto-plastic work) (per unit mass) at each material point [uJ.kg^-1]
         beta                                                            ! dissipation factor [-]

     real(pReal), dimension(n_dir + n_shr) :: &
         sigma_d_n, &                                                    ! trial deviatoric stress tensor (voigt notation, vector form) at the beginning of the time increment [MPa]
         delta_eps_p, &                                                  ! plastic strain increment tensor (voigt notation, vector form) [-]
         delta_eps_e                                                     ! elastic strain increment tensor (voigt notation, vector form) [-]

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
!        calculation of hydrostatic stress   {Eq. Box 1.3(c)}
         sigma_h_n = (sigma_n(i_mp,1) + sigma_n(i_mp,2) + sigma_n(i_mp,3)) / 3.0e0

!        calculation of trial deviatoric stress tensor   {Eq. Box 1.3}
         sigma_d_n(1) = sigma_n(i_mp,1) - sigma_h_n
         sigma_d_n(2) = sigma_n(i_mp,2) - sigma_h_n
         sigma_d_n(3) = sigma_n(i_mp,3) - sigma_h_n
         sigma_d_n(4) = sigma_n(i_mp,4)
         if (n_shr > 1) then
             sigma_d_n(5) = sigma_n(i_mp,5)
             sigma_d_n(6) = sigma_n(i_mp,6)
         end if

!        calculation of equivalent trial stress   {Eq. Box 1.4(b)}
         if (n_shr == 1) then
             sigma_bar_n = sqrt(3.0e0 / 2.0e0 * (sigma_d_n(1)**2 + sigma_d_n(2)**2 + sigma_d_n(3)**2 + 2.0e0 * sigma_d_n(4)**2))
         else
             sigma_bar_n = sqrt(3.0e0 / 2.0e0 * (sigma_d_n(1)**2 + sigma_d_n(2)**2 + sigma_d_n(3)**2 + 2.0e0 * (sigma_d_n(4)**2 + sigma_d_n(5)**2 + sigma_d_n(6)**2)))
         end if

!        calculation of normalized relative temperatures   {Eq. 93(d)}
         T_hat_n = (T_n(i_mp) - T_abs0) / (T_0(i_material) - T_abs0)

!        calculation of plastic/critical stress
         G_alpha_c = G_0(i_material) * alpha_c0(i_material) * (1.0e0 + r_G_alpha_c(i_material) * (T_hat_n - 1.0e0)**s_G_alpha_c(i_material))   ! {Eq. 93(c)}
         G_alpha_w = G_0(i_material) * alpha_w0(i_material) * (1.0e0 + r_G_alpha_w(i_material) * (T_hat_n - 1.0e0)**s_G_alpha_w(i_material))   ! {Eq. 93(c)}
         sigma_p_n = M(i_material) * b(i_material) * G_alpha_c * sqrt(rho_0(i_material) * rho_hat_ci_n) + M(i_material) * b(i_material) * G_alpha_w * sqrt(rho_0(i_material) * rho_hat_wi_n)   ! {Eqs. 92(b) and 93}

!        calculation of viscous stress
         if (sigma_bar_n > sigma_p_n) then
             sigma_v_n = sigma_bar_n - sigma_p_n
         else
             sigma_v_n = 0.0e0
         end if

!        calculation of equivalent plastic strain increment
         m_v      = m_v0     (i_material) * (1.0e0 + r_mv(i_material) * (T_hat_n - 1.0e0)**s_mv(i_material))   ! {Eq. 94(d)}
         sigma_v0 = sigma_v00(i_material) * (1.0e0 + r_v (i_material) * (T_hat_n - 1.0e0)**s_v (i_material))   ! {Eq. 94(b)}
         eps_bar_dot_p   = eps_bar_dot_0(i_material) * (sigma_v_n / sigma_v0)**(1.0e0 / m_v)                   ! {Eq. 109}
         delta_eps_bar_p = delta_t * eps_bar_dot_p                                                             ! {Eq. Box 1.10}

!        calculation of plastic and elastic strain increment tensors
         if (sigma_bar_n == 0.0e0) then
             delta_eps_p = 0.0e0
         else
             delta_eps_p = 3.0e0 / 2.0e0 * sigma_d_n / sigma_bar_n * delta_eps_bar_p   ! {Eq. 83(b)}
         end if
         delta_eps_e = delta_eps(i_mp,:) - delta_eps_p   ! {Eq. 13}

!        calculation of elastic stiffness moduli
         nu = nu_0(i_material) * (1.0e0 + r_nu(i_material) * (T_hat_n - 1.0e0)**s_nu(i_material))   ! {Eq. 118}
         G  = G_0 (i_material) * (1.0e0 + r_G (i_material) * (T_hat_n - 1.0e0)**s_G (i_material))   ! {Eq. 117}
         K  = 2.0e0 / 3.0e0 * (1.0e0 + nu) / (1.0e0 - 2.0e0 * nu) * G                               ! {Eq. Box 1.1(b)}

!        calculation of volumetric strain increment   {Eq. Box 1.2(b)}
         delta_eps_v = delta_eps(i_mp,1) + delta_eps(i_mp,2) + delta_eps(i_mp,3)

!        calculation of stress tensor   {Eq. Box 1.2}
         sigma(i_mp,1) = sigma_n(i_mp,1) + 2.0e0 * G * delta_eps_e(1) + (K - 2.0e0 / 3.0e0 * G) * delta_eps_v
         sigma(i_mp,2) = sigma_n(i_mp,2) + 2.0e0 * G * delta_eps_e(2) + (K - 2.0e0 / 3.0e0 * G) * delta_eps_v
         sigma(i_mp,3) = sigma_n(i_mp,3) + 2.0e0 * G * delta_eps_e(3) + (K - 2.0e0 / 3.0e0 * G) * delta_eps_v
         sigma(i_mp,4) = sigma_n(i_mp,4) + 2.0e0 * G * delta_eps_e(4)
         if (n_shr > 1) then
             sigma(i_mp,5) = sigma_n(i_mp,5) + 2.0e0 * G * delta_eps_e(5)
             sigma(i_mp,6) = sigma_n(i_mp,6) + 2.0e0 * G * delta_eps_e(6)
         end if

!        calculation of strain rate sensitivity parameters associated with different dynamic dislocation processes   {Eq. 95(b)}
         m_an_cm = m_an_cm0(i_material) * (1.0e0 + r_m_an_cm(i_material) * (T_hat_n - 1.0e0)**s_m_an_cm(i_material))
         m_an_ci = m_an_ci0(i_material) * (1.0e0 + r_m_an_ci(i_material) * (T_hat_n - 1.0e0)**s_m_an_ci(i_material))
         m_an_wi = m_an_wi0(i_material) * (1.0e0 + r_m_an_wi(i_material) * (T_hat_n - 1.0e0)**s_m_an_wi(i_material))
         m_tr_cm = m_tr_cm0(i_material) * (1.0e0 + r_m_tr_cm(i_material) * (T_hat_n - 1.0e0)**s_m_tr_cm(i_material))
         m_nc_wi = m_nc_wi0(i_material) * (1.0e0 + r_m_nc_wi(i_material) * (T_hat_n - 1.0e0)**s_m_nc_wi(i_material))
         m_rm_ci = m_rm_ci0(i_material) * (1.0e0 + r_m_rm_ci(i_material) * (T_hat_n - 1.0e0)**s_m_rm_ci(i_material))
         m_rm_wi = m_rm_wi0(i_material) * (1.0e0 + r_m_rm_wi(i_material) * (T_hat_n - 1.0e0)**s_m_rm_wi(i_material))

!        calculation of probability amplitudes associated with different dynamic dislocation processes   {Eqs. 63, 95}
         c_gn_cm = c_gn_cm0(i_material)
         c_an_cm = c_an_cm0(i_material) * (1.0e0 + r_an_cm(i_material)   * (T_hat_n - 1.0e0)**s_an_cm(i_material)) * (delta_eps_bar_p / delta_t / eps_bar_dot_0(i_material))**m_an_cm
         c_an_ci = c_an_ci0(i_material) * (1.0e0 + r_an_ci(i_material)   * (T_hat_n - 1.0e0)**s_an_ci(i_material)) * (delta_eps_bar_p / delta_t / eps_bar_dot_0(i_material))**m_an_ci
         c_an_wi = c_an_wi0(i_material) * (1.0e0 + r_an_wi(i_material)   * (T_hat_n - 1.0e0)**s_an_wi(i_material)) * (delta_eps_bar_p / delta_t / eps_bar_dot_0(i_material))**m_an_wi
         c_ac_ci = c_ac_ci0(i_material)
         c_ac_wi = c_ac_wi0(i_material)
         c_tr_cm = c_tr_cm0(i_material) * (1.0e0 + r_tr_cm(i_material)   * (T_hat_n - 1.0e0)**s_tr_cm(i_material)) * (delta_eps_bar_p / delta_t / eps_bar_dot_0(i_material))**m_tr_cm
         c_nc_wi = c_nc_wi0(i_material) * (1.0e0 + r_nc_wi(i_material)   * (T_hat_n - 1.0e0)**s_nc_wi(i_material)) * (delta_eps_bar_p / delta_t / eps_bar_dot_0(i_material))**m_nc_wi
         c_rm_ci = c_rm_ci0(i_material) * (1.0e0 + r_rm_ci(i_material)   * (T_hat_n - 1.0e0)**s_rm_ci(i_material)) * (delta_eps_bar_p / delta_t / eps_bar_dot_0(i_material))**m_rm_ci
         c_rm_wi = c_rm_wi0(i_material) * (1.0e0 + r_rm_wi(i_material)   * (T_hat_n - 1.0e0)**s_rm_wi(i_material)) * (delta_eps_bar_p / delta_t / eps_bar_dot_0(i_material))**m_rm_wi

!        calculation of dimensionless rates (w.r.t. equivalent plastic strain) or kinetics of different dynamic dislocation processes at the beginning of the time increment   {Eq. 60}
         del_rho_hat_gn_cm_n = M(i_material) * c_gn_cm * rho_hat_cm_n / sqrt(rho_hat_ci_n + rho_hat_wi_n)
         del_rho_hat_an_cm_n = M(i_material) * c_an_cm * rho_hat_cm_n * rho_hat_cm_n
         del_rho_hat_an_ci_n = M(i_material) * c_an_ci * rho_hat_ci_n * rho_hat_cm_n
         del_rho_hat_an_wi_n = M(i_material) * c_an_wi * rho_hat_wi_n * rho_hat_cm_n
         del_rho_hat_ac_ci_n = M(i_material) * c_ac_ci * sqrt(rho_hat_ci_n) * rho_hat_cm_n
         del_rho_hat_ac_wi_n = M(i_material) * c_ac_wi * sqrt(rho_hat_wi_n) * rho_hat_cm_n
         del_rho_hat_tr_cm_n = M(i_material) * c_tr_cm * sqrt(rho_hat_cm_n) * rho_hat_cm_n
         del_rho_hat_nc_wi_n = M(i_material) * c_nc_wi * sqrt(rho_hat_ci_n) * rho_hat_ci_n * rho_hat_cm_n
         del_rho_hat_rm_ci_n = M(i_material) * c_rm_ci * rho_hat_ci_n
         del_rho_hat_rm_wi_n = M(i_material) * c_rm_wi * rho_hat_wi_n

!        calculation of dimensionless rates (w.r.t. equivalent plastic strain) or kinetics of different types of dislocation densities at the beginning of the time increment   {Eq. 59}
         del_rho_hat_cm_n = del_rho_hat_gn_cm_n + del_rho_hat_rm_ci_n + del_rho_hat_rm_wi_n &
                          - (2.0e0 * del_rho_hat_an_cm_n + del_rho_hat_an_ci_n + del_rho_hat_an_wi_n + del_rho_hat_ac_ci_n + del_rho_hat_ac_wi_n + del_rho_hat_tr_cm_n)
         del_rho_hat_ci_n = del_rho_hat_ac_ci_n + del_rho_hat_tr_cm_n - (del_rho_hat_an_ci_n + del_rho_hat_nc_wi_n + del_rho_hat_rm_ci_n)
         del_rho_hat_wi_n = del_rho_hat_ac_wi_n + del_rho_hat_nc_wi_n - (del_rho_hat_an_wi_n + del_rho_hat_rm_wi_n)

!        calculation of microstructural state variables (different types of dislocation density) at the end of the time increment   {Eq. 90}
         rho_hat_cm = rho_hat_cm_n + delta_eps_bar_p * del_rho_hat_cm_n
         rho_hat_ci = rho_hat_ci_n + delta_eps_bar_p * del_rho_hat_ci_n
         rho_hat_wi = rho_hat_wi_n + delta_eps_bar_p * del_rho_hat_wi_n

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

!        calculation of plastic dissipation work and plastic heat generation
         beta      = (2.0e0 * (del_rho_hat_an_cm_n + del_rho_hat_an_ci_n + del_rho_hat_an_wi_n) / del_rho_hat_gn_cm_n)**kappa(i_material)   ! {Eq. 65}
         q_p(i_mp) = q_p_n(i_mp) + beta * sigma_bar_n * delta_eps_bar_p / rho_mass(i_mp)   ! {Eqs. Box 1.9, 68}

!        calculation of equivalent plastic strain rate and equivalent plastic strain
         eps_bar_p = eps_bar_p_n + delta_eps_bar_p   ! {Eq. 68}

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