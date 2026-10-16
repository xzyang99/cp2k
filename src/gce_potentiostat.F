!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2026 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Dynamic grand canonical potentiostats: the electron number N is a fictitious particle
!>        with mass m driven by the force WF - TWF (work function minus its target) and coupled
!>        to a Nose-Hoover chain thermostat. N is propagated once per MD step with
!>        velocity Verlet, split around the SCF that is run at fixed N.
! **************************************************************************************************
MODULE gce_potentiostat

   USE input_constants,                 ONLY: gce_method_direct,&
                                              gce_method_nhc
   USE kinds,                           ONLY: dp
   USE qs_environment_types,            ONLY: get_qs_env,&
                                              qs_environment_type
   USE qs_mo_types,                     ONLY: get_mo_set,&
                                              mo_set_type
   USE scf_control_types,               ONLY: gce_type,&
                                              scf_control_type
#include "./base/base_uses.f90"

   IMPLICIT NONE

   PRIVATE

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'gce_potentiostat'

   PUBLIC :: gce_potentiostat_finish, gce_potentiostat_prepare

CONTAINS

! **************************************************************************************************
!> \brief Propagates N to the current MD time before the SCF, using the work function of the
!>        previous step: thermostat half step, half kick and drift.
!> \param qs_env ...
! **************************************************************************************************
   SUBROUTINE gce_potentiostat_prepare(qs_env)
      TYPE(qs_environment_type), POINTER                 :: qs_env

      INTEGER                                            :: ispin, nelectron
      REAL(KIND=dp)                                      :: dt
      TYPE(gce_type), POINTER                            :: gce
      TYPE(mo_set_type), DIMENSION(:), POINTER           :: mos
      TYPE(scf_control_type), POINTER                    :: scf_control

      CALL get_qs_env(qs_env, scf_control=scf_control, mos=mos)
      gce => scf_control%gce
      IF (.NOT. gce%do_gce .OR. gce%method == gce_method_direct) RETURN

      IF (gce%nelec < 0.0_dp) THEN
         gce%nelec = 0.0_dp
         DO ispin = 1, SIZE(mos)
            CALL get_mo_set(mos(ispin), nelectron=nelectron)
            gce%nelec = gce%nelec + REAL(nelectron, dp)
         END DO
      END IF

      dt = qs_env%sim_time - gce%last_time
      IF (gce%last_time < 0.0_dp) THEN
         ! First force evaluation of a new trajectory
         gce%last_time = qs_env%sim_time
      ELSE IF (dt > 0.0_dp) THEN
         gce%dt = dt
         SELECT CASE (gce%method)
         CASE (gce_method_nhc)
            CALL nhc_half_step(gce)
            CALL kick(gce)
            gce%nelec = gce%nelec + dt*gce%electron_velocity
         END SELECT
      ELSE IF (dt < 0.0_dp .AND. qs_env%sim_time /= 0.0_dp) THEN
         ! dt < 0 at time zero is the initial force evaluation of a restarted MD run
         CPABORT("The MD time moved backwards with respect to the GCE potentiostat.")
      END IF

      CALL check_electron_number(gce%nelec, mos)

   END SUBROUTINE gce_potentiostat_prepare

! **************************************************************************************************
!> \brief Completes the step after the SCF with the half kick from the new work function.
!> \param qs_env ...
! **************************************************************************************************
   SUBROUTINE gce_potentiostat_finish(qs_env)
      TYPE(qs_environment_type), POINTER                 :: qs_env

      TYPE(gce_type), POINTER                            :: gce
      TYPE(scf_control_type), POINTER                    :: scf_control

      CALL get_qs_env(qs_env, scf_control=scf_control)
      gce => scf_control%gce
      IF (gce%dt <= 0.0_dp) RETURN

      CALL kick(gce)
      IF (gce%method == gce_method_nhc) CALL nhc_half_step(gce)
      gce%last_time = qs_env%sim_time
      gce%dt = 0.0_dp

   END SUBROUTINE gce_potentiostat_finish

! **************************************************************************************************
!> \brief Half kick of the electron velocity by the force WF - TWF
!> \param gce ...
! **************************************************************************************************
   SUBROUTINE kick(gce)
      TYPE(gce_type), POINTER                            :: gce

      gce%electron_velocity = gce%electron_velocity + &
                              0.5_dp*gce%dt*(gce%prev_workfunction - gce%target_workfunction)/gce%electron_mass

   END SUBROUTINE kick

! **************************************************************************************************
!> \brief Propagates the Nose-Hoover chain and rescales the electron velocity over half a time
!>        step, using a third-order Suzuki-Yoshida decomposition.
!> \param gce ...
! **************************************************************************************************
   SUBROUTINE nhc_half_step(gce)
      TYPE(gce_type), POINTER                            :: gce

      INTEGER                                            :: i, j, m
      REAL(KIND=dp)                                      :: aa, kinetic_twice, q, scale, w, wdt
      REAL(KIND=dp), DIMENSION(gce%chain_length)         :: g
      REAL(KIND=dp), DIMENSION(3)                        :: weights

      m = gce%chain_length
      q = gce%temperature*gce%timecon**2
      kinetic_twice = gce%electron_mass*gce%electron_velocity**2
      scale = 1.0_dp

      w = 1.0_dp/(2.0_dp - 2.0_dp**(1.0_dp/3.0_dp))
      weights = [w, 1.0_dp - 2.0_dp*w, w]
      g(1) = (kinetic_twice - gce%temperature)/q
      DO j = 2, m
         g(j) = (q*gce%nhc_velocities(j - 1)**2 - gce%temperature)/q
      END DO

      DO i = 1, 3
         wdt = weights(i)*gce%dt
         gce%nhc_velocities(m) = gce%nhc_velocities(m) + 0.25_dp*wdt*g(m)
         DO j = m - 1, 1, -1
            aa = EXP(-0.125_dp*wdt*gce%nhc_velocities(j + 1))
            gce%nhc_velocities(j) = gce%nhc_velocities(j)*aa*aa + 0.25_dp*wdt*g(j)*aa
         END DO

         aa = EXP(-0.5_dp*wdt*gce%nhc_velocities(1))
         scale = scale*aa
         g(1) = (scale*scale*kinetic_twice - gce%temperature)/q

         DO j = 1, m - 1
            aa = EXP(-0.125_dp*wdt*gce%nhc_velocities(j + 1))
            gce%nhc_velocities(j) = gce%nhc_velocities(j)*aa*aa + 0.25_dp*wdt*g(j)*aa
            g(j + 1) = (q*gce%nhc_velocities(j)**2 - gce%temperature)/q
         END DO
         gce%nhc_velocities(m) = gce%nhc_velocities(m) + 0.25_dp*wdt*g(m)
      END DO

      gce%electron_velocity = gce%electron_velocity*scale

   END SUBROUTINE nhc_half_step

! **************************************************************************************************
!> \brief Aborts if the electron number does not fit into the MO space.
!> \param nelec ...
!> \param mos ...
! **************************************************************************************************
   SUBROUTINE check_electron_number(nelec, mos)
      REAL(KIND=dp), INTENT(IN)                          :: nelec
      TYPE(mo_set_type), DIMENSION(:), POINTER           :: mos

      INTEGER                                            :: ispin, nmo
      REAL(KIND=dp)                                      :: capacity, maxocc

      capacity = 0.0_dp
      DO ispin = 1, SIZE(mos)
         CALL get_mo_set(mos(ispin), nmo=nmo, maxocc=maxocc)
         capacity = capacity + REAL(nmo, dp)*maxocc
      END DO
      IF (nelec <= 0.0_dp .OR. nelec >= capacity) &
         CPABORT("The GCE electron number does not fit into the MO space; increase ADDED_MOS.")

   END SUBROUTINE check_electron_number

END MODULE gce_potentiostat
