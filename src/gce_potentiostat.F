!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2026 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Integrates the fractional electron number used by dynamic GCE potentiostats.
! **************************************************************************************************
MODULE gce_potentiostat

   USE input_constants,                 ONLY: gce_method_direct,&
                                              gce_method_langevin,&
                                              gce_method_nhc
   USE kinds,                           ONLY: dp
   USE qs_environment_types,            ONLY: qs_environment_type
   USE qs_mo_types,                     ONLY: get_mo_set
   USE scf_control_types,               ONLY: gce_type
#include "./base/base_uses.f90"

   IMPLICIT NONE

   PRIVATE

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'gce_potentiostat'

   PUBLIC :: gce_potentiostat_finish, &
             gce_potentiostat_prepare

CONTAINS

! **************************************************************************************************
!> \brief Advances the first part of a dynamic potentiostat step before the fixed-N SCF calculation.
!> \param qs_env Quickstep environment
! **************************************************************************************************
   SUBROUTINE gce_potentiostat_prepare(qs_env)

      TYPE(qs_environment_type), POINTER                 :: qs_env

      REAL(KIND=dp)                                      :: c, force_n, random_normal, sigma
      TYPE(gce_type), POINTER                            :: gce

      gce => qs_env%scf_control%gce
      IF (.NOT. gce%do_gce) RETURN
      IF (gce%method == gce_method_direct) RETURN

      IF (gce%last_step < 0 .OR. gce%last_time < 0.0_dp .OR. gce%nelec < 0.0_dp) THEN
         IF (gce%nelec < 0.0_dp) gce%nelec = reference_electron_number(qs_env)
         CALL check_electron_number(qs_env, gce%nelec)
         gce%last_step = qs_env%sim_step
         gce%last_time = qs_env%sim_time
         RETURN
      END IF

      IF (qs_env%sim_step == gce%last_step) RETURN
      IF (qs_env%sim_step < gce%last_step) THEN
         ! MD restart evaluates the initial force before restoring its step and time in qs_env.
         IF (qs_env%sim_step == 0 .AND. qs_env%sim_time == 0.0_dp) RETURN
         CPABORT("The MD step moved backwards relative to the GCE potentiostat state.")
      END IF

      gce%integration_dt = qs_env%sim_time - gce%last_time
      IF (gce%integration_dt <= 0.0_dp) THEN
         CPABORT("The GCE potentiostat requires a positive elapsed MD time.")
      END IF

      force_n = gce%prev_workfunction - gce%target_workfunction

      SELECT CASE (gce%method)
      CASE (gce_method_nhc)
         CALL nose_hoover_chain_step(gce, gce%integration_dt)
         gce%electron_velocity = gce%electron_velocity + &
                                 0.5_dp*gce%integration_dt*force_n/gce%electron_mass
         gce%nelec = gce%nelec + gce%integration_dt*gce%electron_velocity

      CASE (gce_method_langevin)
         ! BAOAB: B-A-O-A before the SCF, followed by the final B in the finish routine.
         gce%electron_velocity = gce%electron_velocity + &
                                 0.5_dp*gce%integration_dt*force_n/gce%electron_mass
         gce%nelec = gce%nelec + 0.5_dp*gce%integration_dt*gce%electron_velocity

         c = EXP(-gce%langevin_gamma*gce%integration_dt)
         sigma = SQRT((1.0_dp - c*c)*gce%temperature/gce%electron_mass)
         random_normal = gce%rng%next()
         gce%random_counter = gce%random_counter + 1
         gce%electron_velocity = c*gce%electron_velocity + sigma*random_normal

         gce%nelec = gce%nelec + 0.5_dp*gce%integration_dt*gce%electron_velocity

      CASE DEFAULT
         CPABORT("Invalid dynamic GCE potentiostat method.")
      END SELECT

      CALL check_electron_number(qs_env, gce%nelec)

   END SUBROUTINE gce_potentiostat_prepare

! **************************************************************************************************
!> \brief Completes the force half-step after a converged fixed-N SCF calculation.
!> \param qs_env Quickstep environment
! **************************************************************************************************
   SUBROUTINE gce_potentiostat_finish(qs_env)

      TYPE(qs_environment_type), POINTER                 :: qs_env

      REAL(KIND=dp)                                      :: force_n
      TYPE(gce_type), POINTER                            :: gce

      gce => qs_env%scf_control%gce
      IF (gce%integration_dt <= 0.0_dp) RETURN

      force_n = gce%prev_workfunction - gce%target_workfunction
      gce%electron_velocity = gce%electron_velocity + &
                              0.5_dp*gce%integration_dt*force_n/gce%electron_mass
      IF (gce%method == gce_method_nhc) THEN
         CALL nose_hoover_chain_step(gce, gce%integration_dt)
      END IF

      gce%last_step = qs_env%sim_step
      gce%last_time = qs_env%sim_time
      gce%integration_dt = 0.0_dp

   END SUBROUTINE gce_potentiostat_finish

! **************************************************************************************************
!> \brief Symmetric third-order Yoshida propagation of one Nose-Hoover chain step.
!> \param gce GCE control and dynamic state
!> \param dt elapsed MD time
! **************************************************************************************************
   SUBROUTINE nose_hoover_chain_step(gce, dt)

      TYPE(gce_type), POINTER                            :: gce
      REAL(KIND=dp), INTENT(IN)                          :: dt

      INTEGER                                            :: i, j, m
      REAL(KIND=dp)                                      :: aa, chain_mass, kinetic_twice, rescale, &
                                                            tmp, wdt
      REAL(KIND=dp), ALLOCATABLE                         :: nhc_forces(:)
      REAL(KIND=dp), DIMENSION(3)                        :: yoshida_dt

      m = gce%chain_length
      ALLOCATE (nhc_forces(m))
      chain_mass = gce%temperature*gce%timecon*gce%timecon
      kinetic_twice = gce%electron_mass*gce%electron_velocity*gce%electron_velocity
      rescale = 1.0_dp

      tmp = 1.0_dp/(2.0_dp - 2.0_dp**(1.0_dp/3.0_dp))
      yoshida_dt = [tmp*dt, (1.0_dp - 2.0_dp*tmp)*dt, tmp*dt]
      nhc_forces(1) = (kinetic_twice - gce%temperature)/chain_mass
      DO j = 2, m
         nhc_forces(j) = (chain_mass*gce%nhc_velocities(j - 1)**2 - gce%temperature)/chain_mass
      END DO

      DO i = 1, 3
         wdt = yoshida_dt(i)
         gce%nhc_velocities(m) = gce%nhc_velocities(m) + 0.25_dp*nhc_forces(m)*wdt

         DO j = 1, m - 1
            aa = EXP(-0.125_dp*wdt*gce%nhc_velocities(m - j + 1))
            gce%nhc_velocities(m - j) = gce%nhc_velocities(m - j)*aa*aa + &
                                        0.25_dp*wdt*nhc_forces(m - j)*aa
         END DO

         aa = EXP(-0.5_dp*wdt*gce%nhc_velocities(1))
         rescale = rescale*aa
         nhc_forces(1) = (rescale*rescale*kinetic_twice - gce%temperature)/chain_mass

         DO j = 1, m - 1
            aa = EXP(-0.125_dp*wdt*gce%nhc_velocities(j + 1))
            gce%nhc_velocities(j) = gce%nhc_velocities(j)*aa*aa + &
                                    0.25_dp*wdt*nhc_forces(j)*aa
            nhc_forces(j + 1) = &
               (chain_mass*gce%nhc_velocities(j)*gce%nhc_velocities(j) - gce%temperature)/chain_mass
         END DO

         gce%nhc_velocities(m) = gce%nhc_velocities(m) + 0.25_dp*nhc_forces(m)*wdt
      END DO

      gce%electron_velocity = gce%electron_velocity*rescale

   END SUBROUTINE nose_hoover_chain_step

! **************************************************************************************************
!> \brief Returns the reference electron count defined by CHARGE and MULTIPLICITY.
!> \param qs_env Quickstep environment
!> \return total reference electron number
! **************************************************************************************************
   FUNCTION reference_electron_number(qs_env) RESULT(nelec)

      TYPE(qs_environment_type), POINTER                 :: qs_env
      REAL(KIND=dp)                                      :: nelec

      INTEGER                                            :: ispin, nelectron

      nelec = 0.0_dp
      DO ispin = 1, SIZE(qs_env%mos)
         CALL get_mo_set(qs_env%mos(ispin), nelectron=nelectron)
         nelec = nelec + REAL(nelectron, KIND=dp)
      END DO

   END FUNCTION reference_electron_number

! **************************************************************************************************
!> \brief Checks whether the requested electron number fits in the allocated MO space.
!> \param qs_env Quickstep environment
!> \param nelec requested electron number
! **************************************************************************************************
   SUBROUTINE check_electron_number(qs_env, nelec)

      TYPE(qs_environment_type), POINTER                 :: qs_env
      REAL(KIND=dp), INTENT(IN)                          :: nelec

      INTEGER                                            :: ispin, nmo
      REAL(KIND=dp)                                      :: capacity, maxocc

      capacity = 0.0_dp
      DO ispin = 1, SIZE(qs_env%mos)
         CALL get_mo_set(qs_env%mos(ispin), nmo=nmo, maxocc=maxocc)
         capacity = capacity + REAL(nmo, KIND=dp)*maxocc
      END DO

      IF (nelec <= 0.0_dp .OR. nelec >= capacity) THEN
         CPABORT("GCE electron number outside MO capacity; increase ADDED_MOS.")
      END IF

   END SUBROUTINE check_electron_number

END MODULE gce_potentiostat
