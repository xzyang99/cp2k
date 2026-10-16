!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2026 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!
! **************************************************************************************************
!> \author  Ziwei Chai
!> \date    09.04.2026
!> \version 1.0
! **************************************************************************************************
MODULE gce_methods

   USE cp_control_types,                ONLY: dft_control_type,&
                                              pcc_control_type
   USE input_constants,                 ONLY: pcc_mode_double,&
                                              pcc_mode_single
   USE kinds,                           ONLY: dp
   USE message_passing,                 ONLY: mp_para_env_type
   USE pw_methods,                      ONLY: pw_axpy,&
                                              pw_integrate_function,&
                                              pw_scale,&
                                              pw_transfer,&
                                              pw_zero
   USE pw_grid_types,                   ONLY: pw_grid_type
   USE pw_pool_types,                   ONLY: pw_pool_type
   USE pw_types,                        ONLY: pw_c1d_gs_type,&
                                              pw_r3d_rs_type
   USE scf_control_types,               ONLY: gce_type
#include "./base/base_uses.f90"

   IMPLICIT NONE

   PRIVATE

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'gce_methods'

   PUBLIC :: gce_reference_layer, planar_averaged_v_hartree_3d, planar_counter_charge

CONTAINS

! **************************************************************************************************
!> \brief map the surface normal to the two in-plane axes and one normal axis
!> \param surf_normal ...
!> \param a ...
!> \param b ...
!> \param c ...
!> \author Ziwei Chai
! **************************************************************************************************
   SUBROUTINE get_planar_axes(surf_normal, a, b, c)

      INTEGER, INTENT(IN)                                :: surf_normal
      INTEGER, INTENT(OUT)                               :: a, b, c

      SELECT CASE (surf_normal)
      CASE (3)
         a = 1
         b = 2
         c = 3
      CASE (1)
         a = 2
         b = 3
         c = 1
      CASE (2)
         a = 3
         b = 1
         c = 2
      CASE DEFAULT
         CPABORT("Invalid planar surface normal.")
      END SELECT

   END SUBROUTINE get_planar_axes

! **************************************************************************************************
!> \brief grid layer along the surface normal at which the GCE reference potential is sampled
!> \param pw_grid real-space grid
!> \param c index of the surface normal
!> \param ref_position distance from the lower cell boundary; negative selects the upper cell boundary
!> \return index of the reference layer
! **************************************************************************************************
   FUNCTION gce_reference_layer(pw_grid, c, ref_position) RESULT(iref)

      TYPE(pw_grid_type), INTENT(IN)                     :: pw_grid
      INTEGER, INTENT(IN)                                :: c
      REAL(KIND=dp), INTENT(IN)                          :: ref_position
      INTEGER                                            :: iref

      IF (ref_position < 0.0_dp) THEN
         iref = pw_grid%bounds(2, c)
      ELSE
         IF (ref_position >= pw_grid%dr(c)*REAL(pw_grid%npts(c), dp)) &
            CPABORT("GCE%REFERENCE_POSITION lies outside the cell.")
         iref = pw_grid%bounds(1, c) + MODULO(NINT(ref_position/pw_grid%dr(c)), pw_grid%npts(c))
      END IF

   END FUNCTION gce_reference_layer

! **************************************************************************************************
!> \brief calculate the planar averaged real space potential (e.g. Hartree potential) along the
!>        surface normal of a slab model and the corresponding reference (electrostatic) energy
!>        at the GCE reference layer.
!> \param v_rspace ...
!> \param dft_control ...
!> \param gce ...
!> \param para_env ...
!> \author Ziwei Chai
! **************************************************************************************************
   SUBROUTINE planar_averaged_v_hartree_3d(v_rspace, dft_control, gce, para_env)

      TYPE(pw_r3d_rs_type), INTENT(IN)                   :: v_rspace
      TYPE(dft_control_type), POINTER                    :: dft_control
      TYPE(gce_type), POINTER                            :: gce
      TYPE(mp_para_env_type), POINTER                    :: para_env

      CHARACTER(LEN=*), PARAMETER :: routineN = 'planar_averaged_v_hartree_3d'

      INTEGER                                            :: a, b, c, handle, iref, x, y, z
      REAL(KIND=dp)                                      :: my_ref_esp, ngrids_in_plane
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:)           :: profile, profile_local

      CALL timeset(routineN, handle)

      IF (gce%do_gce) THEN

         ! a and b cell vectors are parallel to the surface, c vector is not.
         CALL get_planar_axes(dft_control%pcc_control%surf_normal, a, b, c)

         iref = gce_reference_layer(v_rspace%pw_grid, c, gce%ref_position)
         my_ref_esp = 0.0_dp
         z = v_rspace%pw_grid%npts(a)*v_rspace%pw_grid%npts(b)
         ngrids_in_plane = REAL(z, dp)

!$OMP    PARALLEL DO DEFAULT(NONE) &
!$OMP                PRIVATE(x,y,z) &
!$OMP                SHARED(v_rspace,A,b,c,iref)&
!$OMP                REDUCTION(+:my_ref_esp)
         DO y = v_rspace%pw_grid%bounds_local(1, b), v_rspace%pw_grid%bounds_local(2, b)
            DO x = v_rspace%pw_grid%bounds_local(1, a), v_rspace%pw_grid%bounds_local(2, a)
               DO z = v_rspace%pw_grid%bounds_local(1, c), v_rspace%pw_grid%bounds_local(2, c)
                  IF (z == iref) THEN
                     IF (c == 3) THEN
                        my_ref_esp = v_rspace%array(x, y, z) + my_ref_esp
                     ELSE IF (c == 2) THEN
                        my_ref_esp = v_rspace%array(y, z, x) + my_ref_esp
                     ELSE
                        my_ref_esp = v_rspace%array(z, x, y) + my_ref_esp
                     END IF
                  END IF
               END DO
            END DO
         END DO
!$OMP    END PARALLEL DO

         CALL para_env%sum(my_ref_esp)
         my_ref_esp = my_ref_esp/ngrids_in_plane
         gce%ref_esp = my_ref_esp

      END IF

      IF (dft_control%do_paep) THEN

         CALL get_planar_axes(dft_control%paep_control%surf_normal, a, b, c)

         z = v_rspace%pw_grid%npts(a)*v_rspace%pw_grid%npts(b)
         ngrids_in_plane = REAL(z, dp)

         ALLOCATE (profile(v_rspace%pw_grid%bounds(1, c):v_rspace%pw_grid%bounds(2, c)))
         profile(:) = 0.0_dp

!$OMP    PARALLEL    DEFAULT(NONE) &
!$OMP                PRIVATE(x,y,z,profile_local) &
!$OMP                SHARED(v_rspace,A,b,c) &
!$OMP                SHARED(profile)

         ALLOCATE (profile_local(v_rspace%pw_grid%bounds(1, c):v_rspace%pw_grid%bounds(2, c)))
         profile_local(:) = 0.0_dp

!$OMP DO
         DO z = v_rspace%pw_grid%bounds_local(1, c), v_rspace%pw_grid%bounds_local(2, c)
            DO y = v_rspace%pw_grid%bounds_local(1, b), v_rspace%pw_grid%bounds_local(2, b)
               DO x = v_rspace%pw_grid%bounds_local(1, a), v_rspace%pw_grid%bounds_local(2, a)
                  IF (c == 3) profile_local(z) = profile_local(z) + v_rspace%array(x, y, z)
                  IF (c == 2) profile_local(z) = profile_local(z) + v_rspace%array(y, z, x)
                  IF (c == 1) profile_local(z) = profile_local(z) + v_rspace%array(z, x, y)
               END DO
            END DO
         END DO
!$OMP END DO

!$OMP CRITICAL
         profile(:) = profile(:) + profile_local(:)
!$OMP END CRITICAL

         DEALLOCATE (profile_local)

!$OMP    END PARALLEL

         CALL para_env%sum(profile(:))
         profile(:) = profile(:)/ngrids_in_plane

         DEALLOCATE (profile)

      END IF

      CALL timestop(handle)

   END SUBROUTINE planar_averaged_v_hartree_3d

! **************************************************************************************************
!> \brief add the planar counter charge density to the total charge density
!> \param rho_tot_gspace ...
!> \param pcc_env ...
!> \param auxbas_pw_pool ...
!> \author Ziwei Chai
! **************************************************************************************************
   SUBROUTINE planar_counter_charge(rho_tot_gspace, pcc_env, auxbas_pw_pool)

      TYPE(pw_c1d_gs_type), INTENT(INOUT)                :: rho_tot_gspace
      TYPE(pcc_control_type), POINTER                    :: pcc_env
      TYPE(pw_pool_type), INTENT(IN), POINTER            :: auxbas_pw_pool

      CHARACTER(LEN=*), PARAMETER :: routineN = 'planar_counter_charge'

      INTEGER                                            :: a, b, c, handle, lcenter, ucenter, x, y, &
                                                            z
      REAL(KIND=dp)                                      :: dist, dz, gau_a, length, tot_charge, val
      TYPE(pw_c1d_gs_type)                               :: pcc_density_g
      TYPE(pw_r3d_rs_type)                               :: pcc_density_r

      CALL timeset(routineN, handle)

      ! a and b cell vectors are parallel to the surface, c vector is not.
      CALL get_planar_axes(pcc_env%surf_normal, a, b, c)

      CALL auxbas_pw_pool%create_pw(pcc_density_r)
      CALL auxbas_pw_pool%create_pw(pcc_density_g)
      CALL pw_zero(pcc_density_r)
      CALL pw_zero(pcc_density_g)

      tot_charge = pw_integrate_function(rho_tot_gspace, isign=-1)

      pcc_env%charge = tot_charge

      dz = pcc_density_r%pw_grid%dr(c)

      length = dz*REAL(pcc_density_r%pw_grid%npts(c), dp)

      SELECT CASE (pcc_env%mode)
      CASE (pcc_mode_double)
         ucenter = pcc_density_r%pw_grid%bounds(2, c) - INT(pcc_env%dist_edge/dz)
         lcenter = pcc_density_r%pw_grid%bounds(1, c) + INT(pcc_env%dist_edge/dz)
      CASE (pcc_mode_single)
         IF (pcc_env%plane_position < 0.0_dp .OR. pcc_env%plane_position >= length) &
            CPABORT("PLANAR_COUNTER_CHARGE%PLANE_POSITION lies outside the cell.")
         lcenter = pcc_density_r%pw_grid%bounds(1, c) + INT(pcc_env%plane_position/dz)
      END SELECT
      gau_a = 1.0_dp
!$OMP    PARALLEL DO DEFAULT(NONE) &
!$OMP                PRIVATE(dist,x,y,z,val) &
!$OMP                SHARED(pcc_density_r,gau_a,ucenter,lcenter,length,A,b,c,dz,pcc_env)
      DO z = pcc_density_r%pw_grid%bounds_local(1, c), pcc_density_r%pw_grid%bounds_local(2, c)
         IF (pcc_env%mode == pcc_mode_double) THEN
            val = gau_a*EXP(-(dz*REAL(z - ucenter, dp))**2.0_dp/(2.0_dp*pcc_env%gau_c**2.0_dp))
            val = val + gau_a*EXP(-(dz*REAL(z - lcenter, dp))**2.0_dp/(2.0_dp*pcc_env%gau_c**2.0_dp))
         ELSE
            ! minimum image distance, so that the single plane is periodic along the normal
            dist = ABS(dz*REAL(z - lcenter, dp))
            dist = MIN(dist, length - dist)
            val = gau_a*EXP(-dist**2.0_dp/(2.0_dp*pcc_env%gau_c**2.0_dp))
         END IF
         IF (val < 1.0E-15_dp) val = 0.0_dp
         DO y = pcc_density_r%pw_grid%bounds_local(1, b), pcc_density_r%pw_grid%bounds_local(2, b)
            DO x = pcc_density_r%pw_grid%bounds_local(1, a), pcc_density_r%pw_grid%bounds_local(2, a)
               IF (c == 3) pcc_density_r%array(x, y, z) = val
               IF (c == 2) pcc_density_r%array(y, z, x) = val
               IF (c == 1) pcc_density_r%array(z, x, y) = val
            END DO
         END DO
      END DO
!$OMP    END PARALLEL DO

      tot_charge = tot_charge/pw_integrate_function(pcc_density_r, isign=1)
      CALL pw_scale(pcc_density_r, tot_charge)
      CALL pw_transfer(pcc_density_r, pcc_density_g)
      CALL pw_axpy(pcc_density_g, rho_tot_gspace)
      CALL pw_transfer(rho_tot_gspace, pcc_density_r)

      tot_charge = pw_integrate_function(rho_tot_gspace, isign=-1)

      CALL auxbas_pw_pool%give_back_pw(pcc_density_r)
      CALL auxbas_pw_pool%give_back_pw(pcc_density_g)

      CALL timestop(handle)

   END SUBROUTINE planar_counter_charge

END MODULE gce_methods
