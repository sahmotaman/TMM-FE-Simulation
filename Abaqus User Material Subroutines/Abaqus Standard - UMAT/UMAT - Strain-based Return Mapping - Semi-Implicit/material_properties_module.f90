!=======================================================================
! User Material Subroutine for Isotropic Elasto-Visco-Plasticity
! CDD-based Material Model
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
! module header
  module material_properties

!-----------------------------------------------------------------------
!    use of global variables
     use controls
     use strings

!-----------------------------------------------------------------------
!    declaration of modules' private variables
     implicit none

     character(len = 64), dimension(n_materials, n_parameters), private :: &
         parameter_name                                                  ! array containing the names of material parameters/properties

     real(pReal),         dimension(n_materials,n_parameters),  private :: &
         material_parameters                                             ! array containing constitutive parameters/properties associated with each material

!-----------------------------------------------------------------------
!    declaration of variables for material parameters (modules' public variables)
     logical,                                     public, protected :: &
         material_parameters_assigned = .false.                          ! flag determining whether the material parameters for a specific material type is assigned/initialized or not

     real(pReal),         dimension(n_materials), public, protected :: &
         T_0, &                                                          ! reference temperature [C]
         nu_0, &                                                         ! reference Poisson's ratio [-]
         r_nu, &                                                         ! temperature sensitivity coefficient associated to Poisson's ratio [-]
         s_nu, &                                                         ! temperature sensitivity exponent    associated to Poisson's ratio [-]
         G_0, &                                                          ! reference Shear modulus [MPa]
         r_G, &                                                          ! temperature sensitivity coefficient associated to shear modulus [-]
         s_G, &                                                          ! temperature sensitivity exponent    associated to shear modulus [-]
         rho_0, &                                                        ! reference dislocation density [m^-2]
         rho_hat_cm0, &                                                  ! normalized initial cell mobile   dislocation density [-]
         rho_hat_ci0, &                                                  ! normalized initial cell immobile dislocation density [-]
         rho_hat_wi0, &                                                  ! normalized initial wall immobile dislocation density [-]
         alpha_c0, &                                                     ! reference interaction strength associated with cell immobile dislocations [-]
         r_G_alpha_c, &                                                  ! temperature sensitivity coefficient associated with combination of shear modulus and interaction strength of cell immobile dislocations [-]
         s_G_alpha_c, &                                                  ! temperature sensitivity exponent    associated with combination of shear modulus and interaction strength of cell immobile dislocations [-]
         alpha_w0, &                                                     ! reference interaction strength associated with wall immobile dislocations [-]
         r_G_alpha_w, &                                                  ! temperature sensitivity coefficient associated with combination of shear modulus and interaction strength of wall immobile dislocations [-]
         s_G_alpha_w, &                                                  ! temperature sensitivity exponent    associated with combination of shear modulus and interaction strength of wall immobile dislocations [-]
         sigma_v00, &                                                    ! reference viscous stress at reference temperature and strain rate [MPa]
         r_v, &                                                          ! temperature sensitivity coefficient associated with viscous stress [-]
         s_v, &                                                          ! temperature sensitivity exponent    associated with viscous stress [-]
         m_v0, &                                                         ! reference strain rate sensitivity parameter associated with viscous stress [-]
         r_mv, &                                                         ! temperature sensitivity coefficient associated with strain rate sensitivity of viscous stress [-]
         s_mv, &                                                         ! temperature sensitivity exponent    associated with strain rate sensitivity of viscous stress [-]
         eps_bar_dot_0, &                                                ! reference equivalent strain rate [s^-1]
         M, &                                                            ! taylor factor [-]
         b, &                                                            ! burgers length (magnitude of Burgers vector) [m]
         kappa, &                                                        ! material constant associated with dissipation factor [-]
         c_gn_cm0, &                                                     ! reference material parameter associated with dynamic generation   of cell mobile dislocations [-]
         c_an_cm0, &                                                     ! reference material parameter associated with dynamic annihilation of cell mobile dislocations [-]
         r_an_cm, &                                                      ! temperature sensitivity coefficient associated with dynamic annihilation of cell mobile dislocations [-]
         s_an_cm, &                                                      ! temperature sensitivity exponent    associated with dynamic annihilation of cell mobile dislocations [-]
         m_an_cm0, &                                                     ! reference strain rate sensitivity parameter associated with dynamic annihilation of cell mobile dislocations [-]
         r_m_an_cm, &                                                    ! temperature sensitivity coefficient associated with strain rate sensitivity of dynamic annihilation of cell mobile dislocations [-]
         s_m_an_cm, &                                                    ! temperature sensitivity exponent    associated with strain rate sensitivity of dynamic annihilation of cell mobile dislocations [-]
         c_an_ci0, &                                                     ! reference material parameter associated with dynamic annihilation of cell immobile dislocations [-]
         r_an_ci, &                                                      ! temperature sensitivity coefficient associated with dynamic annihilation of cell immobile dislocations [-]
         s_an_ci, &                                                      ! temperature sensitivity exponent    associated with dynamic annihilation of cell immobile dislocations [-]
         m_an_ci0, &                                                     ! reference strain rate sensitivity parameter associated with dynamic annihilation of cell immobile dislocations [-]
         r_m_an_ci, &                                                    ! temperature sensitivity coefficient associated with strain rate sensitivity of dynamic annihilation of cell immobile dislocations [-]
         s_m_an_ci, &                                                    ! temperature sensitivity exponent    associated with strain rate sensitivity of dynamic annihilation of cell immobile dislocations [-]
         c_an_wi0, &                                                     ! reference material parameter associated with dynamic annihilation of wall immobile dislocations [-]
         r_an_wi, &                                                      ! temperature sensitivity coefficient associated with dynamic annihilation of wall immobile dislocations [-]
         s_an_wi, &                                                      ! temperature sensitivity exponent    associated with dynamic annihilation of wall immobile dislocations [-]
         m_an_wi0, &                                                     ! reference strain rate sensitivity parameter associated with dynamic annihilation of wall immobile dislocations [-]
         r_m_an_wi, &                                                    ! temperature sensitivity coefficient associated with strain rate sensitivity of dynamic annihilation of wall immobile dislocations [-]
         s_m_an_wi, &                                                    ! temperature sensitivity exponent    associated with strain rate sensitivity of dynamic annihilation of wall immobile dislocations [-]
         c_ac_ci0, &                                                     ! reference material parameter associated with dynamic accumulation of cell immobile dislocations [-]
         c_ac_wi0, &                                                     ! reference material parameter associated with dynamic accumulation of wall immobile dislocations [-]
         c_tr_cm0, &                                                     ! reference material parameter associated with dynamic trapping of cell mobile dislocations [-]
         r_tr_cm, &                                                      ! temperature sensitivity coefficient associated with dynamic trapping of cell mobile dislocations [-]
         s_tr_cm, &                                                      ! temperature sensitivity exponent    associated with dynamic trapping of cell mobile dislocations [-]
         m_tr_cm0, &                                                     ! reference strain rate sensitivity parameter associated with dynamic trapping of cell mobile dislocations [-]
         r_m_tr_cm, &                                                    ! temperature sensitivity coefficient associated with strain rate sensitivity of dynamic trapping of cell mobile dislocations [-]
         s_m_tr_cm, &                                                    ! temperature sensitivity exponent    associated with strain rate sensitivity of dynamic trapping of cell mobile dislocations [-]
         c_nc_wi0, &                                                     ! reference material parameter associated with dynamic trapping of cell immobile dislocations [-]
         r_nc_wi, &                                                      ! temperature sensitivity coefficient associated with dynamic nucleation of wall immobile dislocations [-]
         s_nc_wi, &                                                      ! temperature sensitivity exponent    associated with dynamic nucleation of wall immobile dislocations [-]
         m_nc_wi0, &                                                     ! reference strain rate sensitivity parameter associated with dynamic nucleation of wall immobile dislocations [-]
         r_m_nc_wi, &                                                    ! temperature sensitivity coefficient associated with strain rate sensitivity of dynamic nucleation of wall immobile dislocations [-]
         s_m_nc_wi, &                                                    ! temperature sensitivity exponent    associated with strain rate sensitivity of dynamic nucleation of wall immobile dislocations [-]
         c_rm_ci0, &                                                     ! reference material parameter associated with dynamic remobilization of cell immobile dislocations [-]
         r_rm_ci, &                                                      ! temperature sensitivity coefficient associated with dynamic remobilization of cell immobile dislocations [-]
         s_rm_ci, &                                                      ! temperature sensitivity exponent    associated with dynamic remobilization of cell immobile dislocations [-]
         m_rm_ci0, &                                                     ! reference strain rate sensitivity parameter associated with dynamic remobilization of cell immobile dislocations [-]
         r_m_rm_ci, &                                                    ! temperature sensitivity coefficient associated with strain rate sensitivity of dynamic remobilization of cell immobile dislocations [-]
         s_m_rm_ci, &                                                    ! temperature sensitivity exponent    associated with strain rate sensitivity of dynamic remobilization of cell immobile dislocations [-]
         c_rm_wi0, &                                                     ! reference material parameter associated with dynamic remobilization of wall immobile dislocations [-]
         r_rm_wi, &                                                      ! temperature sensitivity coefficient associated with dynamic remobilization of wall immobile dislocations [-]
         s_rm_wi, &                                                      ! temperature sensitivity exponent    associated with dynamic remobilization of wall immobile dislocations [-]
         m_rm_wi0, &                                                     ! reference strain rate sensitivity parameter associated with dynamic remobilization of wall immobile dislocations [-]
         r_m_rm_wi, &                                                    ! temperature sensitivity coefficient associated with strain rate sensitivity of dynamic remobilization of wall immobile dislocations [-]
         s_m_rm_wi                                                       ! temperature sensitivity exponent    associated with strain rate sensitivity of dynamic remobilization of wall immobile dislocations [-]

     contains

!-----------------------------------------------------------------------
!        subroutine for reading constitutive parameters/properties of each material
         subroutine read_material_parameters()

!-----------------------------------------------------------------------
!            declaration of subroutines' local variables
             implicit none

             character(len = 512) :: &
                 output_dir, &                                           ! path of output/working/job directory
                 path, &                                                 ! directory path of the input file
                 line                                                    ! string containing the read line from the corresponding input file

             character(len = 64), dimension(5) :: &
                 tokens                                                  ! array of strings containing the individual strings on the line separated by delimiter characters

             integer(pInt) :: &
                 length_ouput_dir, &                                     ! length of the character string output_dir
                 i_material, &                                           ! index number associated with each material
                 io_status, &                                            ! error flag
                 i, &                                                    ! loop counter variable i
                 i_parameter, &                                          ! index number associated with each constitutive parameter
                 brackets_position, &                                    ! brackets (), [], {}, <> position
                 n_tokens                                                ! number of found strings on the line separated by delimiter characters

!-----------------------------------------------------------------------
!            opening the input file
             call getoutdir(output_dir, length_ouput_dir)
             path = trim(output_dir) // trim(directory_separator) // trim(material_data_file_name)
             open(unit = 20, file = path, status = 'old', action = 'read', iostat = io_status)

!-----------------------------------------------------------------------
!            reading the input file and storing the parameters
             material_parameters = 0.0
             i_material = 0
             if (io_status == 0) then
                 i_parameter = 0
                 do i = 1, 500
                     brackets_position = 0

!                    reading line from the file
                     call readline(20, line, io_status)

!                    check for (), [], {}, <>
                     call match(line, 1, brackets_position)
                     if (brackets_position > 1) then

!                        parsing the line using delimiter characters
                         call parse(line, delimiters, tokens, n_tokens)

!                        converting the corresponding string number to real number
                         call value(tokens(3), i_material, io_status)

                         i_parameter = 1
                         n_tokens    = 0
                         tokens      = ''
                         cycle
                     end if

!                    check for end of file
                     if (io_status < 0) then
                         close(unit = 20)
                         exit

!                    error: read line error
                     else if (io_status > 0) then
                         write(6,*) 'error: could not read line from the material data file.'
                         close(unit = 20)

!                    no end of file and no read line error
                     else if (i_parameter /= 0) then

!                        parsing the lines under the brackets
                         call parse(line, delimiters, tokens, n_tokens)

!                        storing parsed values
                         parameter_name(i_material,i_parameter) = trim(tokens(1))
                         call value(tokens(2), material_parameters(i_material,i_parameter), io_status)

                         i_parameter = i_parameter + 1

!                    error: no brackets found
                     else
                         write(6,*) 'error: no material definition was found.'
                         close(unit = 20)
                     end if
                 end do

!            error: file could not be opened
             else
                 write(6,*) 'error: file could not be opened.'
             end if

         end subroutine read_material_parameters

!-----------------------------------------------------------------------
!        subroutine for assigning values to each constitutive parameter for each material
         subroutine assign_material_parameters()

!-----------------------------------------------------------------------
!            declaration of subroutines' local variables
             implicit none

             integer(pInt) :: &
                 i_material, &                                           ! index number associated with each material
                 i_parameter                                             ! index number associated with each constitutive parameter

!-----------------------------------------------------------------------
!            reading material parameters of each material type from the material data input file
             call read_material_parameters()

!-----------------------------------------------------------------------
!            storing constitutive parameters of each material
             do i_material  = 1, n_materials
                 do i_parameter = 1, n_parameters
                     select case (parameter_name(i_material,i_parameter))
                         case ('T_0')
                             T_0(i_material)           = material_parameters(i_material,i_parameter)
                         case ('nu_0')
                             nu_0(i_material)          = material_parameters(i_material,i_parameter)
                         case ('r_nu')
                             r_nu(i_material)          = material_parameters(i_material,i_parameter)
                         case ('s_nu')
                             s_nu(i_material)          = material_parameters(i_material,i_parameter)
                         case ('G_0')
                             G_0(i_material)           = material_parameters(i_material,i_parameter)
                         case ('r_G')
                             r_G(i_material)           = material_parameters(i_material,i_parameter)
                         case ('s_G')
                             s_G(i_material)           = material_parameters(i_material,i_parameter)
                         case ('rho_0')
                             rho_0(i_material)         = material_parameters(i_material,i_parameter)
                         case ('rho_hat_cm0')
                             rho_hat_cm0(i_material)   = material_parameters(i_material,i_parameter)
                         case ('rho_hat_ci0')
                             rho_hat_ci0(i_material)   = material_parameters(i_material,i_parameter)
                         case ('rho_hat_wi0')
                             rho_hat_wi0(i_material)   = material_parameters(i_material,i_parameter)
                         case ('alpha_c0')
                             alpha_c0(i_material)      = material_parameters(i_material,i_parameter)
                         case ('r_G_alpha_c')
                             r_G_alpha_c(i_material)   = material_parameters(i_material,i_parameter)
                         case ('s_G_alpha_c')
                             s_G_alpha_c(i_material)   = material_parameters(i_material,i_parameter)
                         case ('alpha_w0')
                             alpha_w0(i_material)      = material_parameters(i_material,i_parameter)
                         case ('r_G_alpha_w')
                             r_G_alpha_w(i_material)   = material_parameters(i_material,i_parameter)
                         case ('s_G_alpha_w')
                             s_G_alpha_w(i_material)   = material_parameters(i_material,i_parameter)
                         case ('sigma_v00')
                             sigma_v00(i_material)     = material_parameters(i_material,i_parameter)
                         case ('r_v')
                             r_v(i_material)           = material_parameters(i_material,i_parameter)
                         case ('s_v')
                             s_v(i_material)           = material_parameters(i_material,i_parameter)
                         case ('m_v0')
                             m_v0(i_material)          = material_parameters(i_material,i_parameter)
                         case ('r_mv')
                             r_mv(i_material)          = material_parameters(i_material,i_parameter)
                         case ('s_mv')
                             s_mv(i_material)          = material_parameters(i_material,i_parameter)
                         case ('eps_bar_dot_0')
                             eps_bar_dot_0(i_material) = material_parameters(i_material,i_parameter)
                         case ('M')
                             M(i_material)             = material_parameters(i_material,i_parameter)
                         case ('b')
                             b(i_material)             = material_parameters(i_material,i_parameter)
                         case ('c_gn_cm0')
                             c_gn_cm0(i_material)      = material_parameters(i_material,i_parameter)
                         case ('c_an_cm0')
                             c_an_cm0(i_material)      = material_parameters(i_material,i_parameter)
                         case ('r_an_cm')
                             r_an_cm(i_material)       = material_parameters(i_material,i_parameter)
                         case ('s_an_cm')
                             s_an_cm(i_material)       = material_parameters(i_material,i_parameter)
                         case ('m_an_cm0')
                             m_an_cm0(i_material)      = material_parameters(i_material,i_parameter)
                         case ('r_m_an_cm')
                             r_m_an_cm(i_material)     = material_parameters(i_material,i_parameter)
                         case ('s_m_an_cm')
                             s_m_an_cm(i_material)     = material_parameters(i_material,i_parameter)
                         case ('c_an_ci0')
                             c_an_ci0(i_material)      = material_parameters(i_material,i_parameter)
                         case ('r_an_ci')
                             r_an_ci(i_material)       = material_parameters(i_material,i_parameter)
                         case ('s_an_ci')
                             s_an_ci(i_material)       = material_parameters(i_material,i_parameter)
                         case ('m_an_ci0')
                             m_an_ci0(i_material)      = material_parameters(i_material,i_parameter)
                         case ('r_m_an_ci')
                             r_m_an_ci(i_material)     = material_parameters(i_material,i_parameter)
                         case ('s_m_an_ci')
                             s_m_an_ci(i_material)     = material_parameters(i_material,i_parameter)
                         case ('c_an_wi0')
                             c_an_wi0(i_material)      = material_parameters(i_material,i_parameter)
                         case ('r_an_wi')
                             r_an_wi(i_material)       = material_parameters(i_material,i_parameter)
                         case ('s_an_wi')
                             s_an_wi(i_material)       = material_parameters(i_material,i_parameter)
                         case ('m_an_wi0')
                             m_an_wi0(i_material)      = material_parameters(i_material,i_parameter)
                         case ('r_m_an_wi')
                             r_m_an_wi(i_material)     = material_parameters(i_material,i_parameter)
                         case ('s_m_an_wi')
                             s_m_an_wi(i_material)     = material_parameters(i_material,i_parameter)
                         case ('c_ac_ci0')
                             c_ac_ci0(i_material)      = material_parameters(i_material,i_parameter)
                         case ('c_ac_wi0')
                             c_ac_wi0(i_material)      = material_parameters(i_material,i_parameter)
                         case ('c_tr_cm0')
                             c_tr_cm0(i_material)      = material_parameters(i_material,i_parameter)
                         case ('r_tr_cm')
                             r_tr_cm(i_material)       = material_parameters(i_material,i_parameter)
                         case ('s_tr_cm')
                             s_tr_cm(i_material)       = material_parameters(i_material,i_parameter)
                         case ('m_tr_cm0')
                             m_tr_cm0(i_material)      = material_parameters(i_material,i_parameter)
                         case ('r_m_tr_cm')
                             r_m_tr_cm(i_material)     = material_parameters(i_material,i_parameter)
                         case ('s_m_tr_cm')
                             s_m_tr_cm(i_material)     = material_parameters(i_material,i_parameter)
                         case ('c_nc_wi0')
                             c_nc_wi0(i_material)      = material_parameters(i_material,i_parameter)
                         case ('r_nc_wi')
                             r_nc_wi(i_material)       = material_parameters(i_material,i_parameter)
                         case ('s_nc_wi')
                             s_nc_wi(i_material)       = material_parameters(i_material,i_parameter)
                         case ('m_nc_wi0')
                             m_nc_wi0(i_material)      = material_parameters(i_material,i_parameter)
                         case ('r_m_nc_wi')
                             r_m_nc_wi(i_material)     = material_parameters(i_material,i_parameter)
                         case ('s_m_nc_wi')
                             s_m_nc_wi(i_material)     = material_parameters(i_material,i_parameter)
                         case ('c_rm_ci0')
                             c_rm_ci0(i_material)      = material_parameters(i_material,i_parameter)
                         case ('r_rm_ci')
                             r_rm_ci(i_material)       = material_parameters(i_material,i_parameter)
                         case ('s_rm_ci')
                             s_rm_ci(i_material)       = material_parameters(i_material,i_parameter)
                         case ('m_rm_ci0')
                             m_rm_ci0(i_material)      = material_parameters(i_material,i_parameter)
                         case ('r_m_rm_ci')
                             r_m_rm_ci(i_material)     = material_parameters(i_material,i_parameter)
                         case ('s_m_rm_ci')
                             s_m_rm_ci(i_material)     = material_parameters(i_material,i_parameter)
                         case ('c_rm_wi0')
                             c_rm_wi0(i_material)      = material_parameters(i_material,i_parameter)
                         case ('r_rm_wi')
                             r_rm_wi(i_material)       = material_parameters(i_material,i_parameter)
                         case ('s_rm_wi')
                             s_rm_wi(i_material)       = material_parameters(i_material,i_parameter)
                         case ('m_rm_wi0')
                             m_rm_wi0(i_material)      = material_parameters(i_material,i_parameter)
                         case ('r_m_rm_wi')
                             r_m_rm_wi(i_material)     = material_parameters(i_material,i_parameter)
                         case ('s_m_rm_wi')
                             s_m_rm_wi(i_material)     = material_parameters(i_material,i_parameter)
                         case ('kappa')
                             kappa(i_material)         = material_parameters(i_material,i_parameter)
                     end select
                 end do
             end do
             material_parameters_assigned = .true.

         end subroutine assign_material_parameters
 
  end module material_properties