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
!    doi: 10.1016/j.jmps.2018.09.002
!
!    Motaman, S.A.H.; Schacht K.; Haase, C.; Prahl, U.; 2019.
!    Thermo-micro-mechanical simulation of bulk metal forming processes.
!=======================================================================

!***********************************************************************
! module header
  module controls

!-----------------------------------------------------------------------
!    setting precision controlling parameters
     implicit none

!    setting precision for real and integer kinds
     integer,                     parameter, public :: & 
         pReal = 8, &
         pInt  = 4

!-----------------------------------------------------------------------
!    assignment of numerics (return mapping) controlling parameters
     real(pReal),                 parameter, public :: &
         chi      = 1.0e-6, &                                            ! normalized tolerance for return mapping newton-raphson loop
         xi_min   = 1.0e-2, &                                            ! the parameter controlling the minimum equivalent plastic strain rate [s^-1]
         xi_mean  = 1.0e2                                                ! the parameter controlling the mean    equivalent plastic strain rate [s^-1]

     integer(pInt),               parameter, public :: &
         k_NR_max = 20                                                   ! maximum loop iterations in newton-raphson (NR) scheme

!-----------------------------------------------------------------------
!    assignment of parameters for extraction of material constants
     character(len = 64),         parameter, public :: &
         directory_separator     = '/', &                                ! directory separator character
         material_data_file_name = 'material_data.inp', &                ! the name of phase data input file
         delimiters              = '= :	,'                               ! delimiter characters for parsing input file lines

     integer(pInt),               parameter, public :: &
         n_materials  = 2, &                                             ! maximum number of different materials
         n_parameters = 80                                               ! maximum number of properties of different materials

!-----------------------------------------------------------------------
!    assignment of other parameters
     real(pReal),                 parameter, public :: &
         T_abs0 = -273.15e0                                              ! absolute zero temperature [C]

  end module controls