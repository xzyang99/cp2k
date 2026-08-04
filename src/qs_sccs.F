!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2026 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief   Self-consistent continuum solvation (SCCS) model implementation
!> \author  Matthias Krack (MK)
!> \version 1.0
!> \par Literature:
!>      - J.-L. Fattebert and F. Gygi,
!>        Density functional theory for efficient ab initio molecular dynamics
!>        simulations in solution, J. Comput. Chem. 23, 662-666 (2002)
!>      - O. Andreussi, I. Dabo, and N. Marzari,
!>        Revised self-consistent continuum solvation in electronic-structure
!>        calculations, J. Chem. Phys. 136, 064102-20 (2012)
!> \par History:
!>      - Creation (10.10.2013,MK)
!>      - Derivatives using finite differences (26.11.2013,MK)
!>      - Cube file dump of the dielectric function (19.12.2013,MK)
!>      - Cube file dump of the polarisation potential (20.12.2013,MK)
!>      - Calculation of volume and surface of the cavity (21.12.2013,MK)
!>      - Functional derivative of the cavitation energy (28.12.2013,MK)
!>      - Update printout (11.11.2022,MK)
! **************************************************************************************************

MODULE qs_sccs

   USE atomic_kind_types,               ONLY: atomic_kind_type,&
                                              get_atomic_kind
   USE cell_types,                      ONLY: cell_type,&
                                              get_cell
   USE cp_control_types,                ONLY: dft_control_type,&
                                              sccs_control_type
   USE cp_log_handling,                 ONLY: cp_get_default_logger,&
                                              cp_logger_type
   USE cp_output_handling,              ONLY: cp_p_file,&
                                              cp_print_key_finished_output,&
                                              cp_print_key_should_output,&
                                              cp_print_key_unit_nr,&
                                              low_print_level
   USE cp_realspace_grid_cube,          ONLY: cp_pw_to_cube
   USE cp_realspace_grid_init,          ONLY: init_input_type
   USE cp_subsys_types,                 ONLY: cp_subsys_get,&
                                              cp_subsys_type
   USE cp_units,                        ONLY: cp_unit_from_cp2k
   USE input_constants,                 ONLY: sccs_andreussi,&
                                              sccs_derivative_cd3,&
                                              sccs_derivative_cd5,&
                                              sccs_derivative_cd7,&
                                              sccs_derivative_fft,&
                                              sccs_fattebert_gygi,&
                                              sccs_nonlinear_electrolyte,&
                                              sccs_saa_andreussi
   USE input_section_types,             ONLY: section_get_ivals,&
                                              section_get_lval,&
                                              section_vals_get_subs_vals,&
                                              section_vals_type
   USE kinds,                           ONLY: default_path_length,&
                                              default_string_length,&
                                              dp,&
                                              int_8
   USE mathconstants,                   ONLY: fourpi,&
                                              pi,&
                                              twopi
   USE message_passing,                 ONLY: mp_para_env_type
   USE particle_list_types,             ONLY: particle_list_type
   USE physcon,                         ONLY: a_bohr,&
                                              n_avogadro
   USE pw_env_types,                    ONLY: pw_env_get,&
                                              pw_env_type
   USE pw_methods,                      ONLY: &
        pw_axpy, pw_copy, pw_derive, pw_integral_ab, pw_integrate_function, pw_laplace, pw_scale, &
        pw_transfer, pw_zero
   USE pw_poisson_methods,              ONLY: pw_func_u_convolution,&
                                              pw_poisson_solve
   USE pw_poisson_types,                ONLY: pw_poisson_analytic,&
                                              pw_poisson_mt,&
                                              pw_poisson_type
   USE pw_pool_types,                   ONLY: pw_pool_p_type,&
                                              pw_pool_type
   USE pw_types,                        ONLY: pw_c1d_gs_type,&
                                              pw_r3d_rs_type
   USE qs_energy_types,                 ONLY: qs_energy_type
   USE qs_environment_types,            ONLY: get_qs_env,&
                                              qs_environment_type
   USE qs_kind_types,                   ONLY: get_qs_kind,&
                                              qs_kind_type
   USE qs_rho_types,                    ONLY: qs_rho_get,&
                                              qs_rho_type
   USE qs_scf_types,                    ONLY: qs_scf_env_type
   USE realspace_grid_types,            ONLY: realspace_grid_desc_type,&
                                              realspace_grid_input_type,&
                                              realspace_grid_type,&
                                              rs_grid_create,&
                                              rs_grid_create_descriptor,&
                                              rs_grid_release,&
                                              rs_grid_release_descriptor
   USE rs_methods,                      ONLY: derive_fdm_cd3,&
                                              derive_fdm_cd5,&
                                              derive_fdm_cd7
   USE scf_control_types,               ONLY: scf_control_type
#include "./base/base_uses.f90"

   IMPLICIT NONE

   PRIVATE

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'qs_sccs'
   REAL(KIND=dp), PARAMETER, PRIVATE     :: nle_exp_cutoff = 100.0_dp, &
                                            nle_n_min = 1.0E-4_dp, &
                                            nle_strict_scf_factor = 5.0_dp

   PUBLIC :: print_sccs_results, sccs

CONTAINS

! **************************************************************************************************
!> \brief   Self-consistent continuum solvation (SCCS) model implementation
!> \param qs_env ...
!> \param rho_tot_gspace ...
!> \param v_hartree_gspace ...
!> \param v_sccs ...
!> \param h_stress ...
!> \par History:
!>      - Creation (10.10.2013,MK)
!> \author  Matthias Krack (MK)
!> \version 1.0
! **************************************************************************************************

   SUBROUTINE sccs(qs_env, rho_tot_gspace, v_hartree_gspace, v_sccs, h_stress)

      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(pw_c1d_gs_type), INTENT(INOUT)                :: rho_tot_gspace, v_hartree_gspace
      TYPE(pw_r3d_rs_type), INTENT(INOUT)                :: v_sccs
      REAL(KIND=dp), DIMENSION(3, 3), INTENT(OUT), &
         OPTIONAL                                        :: h_stress

      CHARACTER(LEN=*), PARAMETER                        :: routineN = 'sccs'
      REAL(KIND=dp), PARAMETER                           :: epstol = 1.0E-8_dp

      CHARACTER(LEN=4*default_string_length)             :: message, my_pos_cube
      CHARACTER(LEN=default_path_length)                 :: cube_path, filename, mpi_filename, &
                                                            print_path
      INTEGER :: cube_unit, handle, i, iat, iatom, ikind, ispin, iter, j, k, nat, natom, nkind, &
         nle_max_outer, nle_pcg_iterations, nspin, output_unit, print_level
      INTEGER(KIND=int_8)                                :: ngpts
      INTEGER, DIMENSION(3)                              :: cell_periodic, lb, ub
      LOGICAL :: append_cube, calculate_stress_tensor, do_kpoints, do_nle, mpi_io, nle_converged, &
         nle_strict_solve, should_output
      REAL(KIND=dp) :: alpha_zeta, beta_phi, cavity_surface, cavity_volume, cell_volume, &
         charge_correction, delta_eta, delta_zeta, dphi2, dvol, e_nle_diel, e_nle_field, &
         e_nle_ion, e_nle_shift, e_tot, electrolyte_charge, eps_sccs_current, epsilon_solvent, f, &
         f0, ion_access_volume, ion_bulk_density, ion_capacity, ion_charge_target, ion_denom, &
         ion_diameter, ion_free_density, ion_packing, ion_phi_shift, kT, nle_residual_rms, &
         polarisation_charge, r, R_solv, response_charge, rho_delta, rho_delta_avg, rho_delta_max, &
         rho_iter_new, solute_charge, tot_rho_elec, tot_rho_solute, zeff
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:)           :: zcharge
      REAL(KIND=dp), DIMENSION(3)                        :: abc, lxyz, uxyz
      TYPE(atomic_kind_type), DIMENSION(:), POINTER      :: atomic_kind_set
      TYPE(cell_type), POINTER                           :: cell
      TYPE(cp_logger_type), POINTER                      :: logger
      TYPE(cp_subsys_type), POINTER                      :: cp_subsys
      TYPE(dft_control_type), POINTER                    :: dft_control
      TYPE(mp_para_env_type), POINTER                    :: para_env
      TYPE(particle_list_type), POINTER                  :: particles
      TYPE(pw_c1d_gs_type)                               :: nle_kernel_b, nle_kernel_cav, &
                                                            nle_kernel_diel, nle_kernel_ion, &
                                                            nle_kernel_solv
      TYPE(pw_env_type), POINTER                         :: pw_env
      TYPE(pw_poisson_type), POINTER                     :: poisson_env
      TYPE(pw_pool_p_type), DIMENSION(:), POINTER        :: pw_pools
      TYPE(pw_pool_type), POINTER                        :: auxbas_pw_pool
      TYPE(pw_r3d_rs_type) :: d_s_rhoel, deps_elec, dtf, eps_elec, ff, nle_lambda_cav, &
         nle_lambda_diel, nle_lambda_ion, nle_rho_bound, nle_rho_ion, nle_scav, nle_sdiel, &
         nle_sdiel_eff, nle_sion, nle_ssolv, nle_svdw, nle_xcav, nle_xdiel, nle_xion, nle_xsolv, &
         nle_xvdw, p_e_interface, s, t, tem_convolution, tem_func, u
      TYPE(pw_r3d_rs_type), DIMENSION(3)                 :: dln_eps_elec, dphi_tot, nle_eloc
      TYPE(pw_r3d_rs_type), DIMENSION(:), POINTER        :: rho_pw_r
      TYPE(pw_r3d_rs_type), POINTER                      :: rho_pw_r_sccs
      TYPE(qs_energy_type), POINTER                      :: energy
      TYPE(qs_kind_type), DIMENSION(:), POINTER          :: qs_kind_set
      TYPE(qs_rho_type), POINTER                         :: rho
      TYPE(qs_scf_env_type), POINTER                     :: scf_env
      TYPE(sccs_control_type), POINTER                   :: sccs_control
      TYPE(scf_control_type), POINTER                    :: scf_control
      TYPE(section_vals_type), POINTER                   :: input

      CALL timeset(routineN, handle)

      NULLIFY (auxbas_pw_pool)
      NULLIFY (atomic_kind_set)
      NULLIFY (cell)
      NULLIFY (cp_subsys)
      NULLIFY (dft_control)
      NULLIFY (energy)
      NULLIFY (input)
      NULLIFY (logger)
      NULLIFY (para_env)
      NULLIFY (particles)
      NULLIFY (poisson_env)
      NULLIFY (pw_env)
      NULLIFY (pw_pools)
      NULLIFY (qs_kind_set)
      NULLIFY (rho)
      NULLIFY (sccs_control)
      NULLIFY (scf_env)
      NULLIFY (scf_control)

      ! Load data from Quickstep environment
      CALL get_qs_env(qs_env=qs_env, &
                      atomic_kind_set=atomic_kind_set, &
                      cp_subsys=cp_subsys, &
                      do_kpoints=do_kpoints, &
                      dft_control=dft_control, &
                      energy=energy, &
                      input=input, &
                      natom=natom, &
                      nkind=nkind, &
                      para_env=para_env, &
                      pw_env=pw_env, &
                      qs_kind_set=qs_kind_set, &
                      rho=rho, &
                      scf_env=scf_env, &
                      scf_control=scf_control)
      CALL cp_subsys_get(cp_subsys, particles=particles, cell=cell)

      sccs_control => dft_control%sccs_control
      do_nle = (sccs_control%method_id == sccs_nonlinear_electrolyte)
      IF (do_nle) CALL get_cell(cell=cell, periodic=cell_periodic)

      ALLOCATE (zcharge(natom))
      DO ikind = 1, nkind
         CALL get_qs_kind(qs_kind_set(ikind), zeff=zeff)
         CALL get_atomic_kind(atomic_kind_set(ikind), natom=nat)
         DO iatom = 1, nat
            iat = atomic_kind_set(ikind)%atom_list(iatom)
            zcharge(iat) = zeff
         END DO
      END DO

      CPASSERT(ASSOCIATED(qs_env))

      IF (do_kpoints) THEN
         CPWARN("SCCS with k-points has not yet been fully validated")
      END IF

      IF (PRESENT(h_stress)) THEN
         IF (do_nle) THEN
            CPABORT("The stress tensor is not implemented for SCCS METHOD NONLINEAR_ELECTROLYTE.")
         END IF
         calculate_stress_tensor = .TRUE.
         h_stress(:, :) = 0.0_dp
         CPWARN("The stress tensor for SCCS has not yet been fully validated")
      ELSE
         calculate_stress_tensor = .FALSE.
      END IF

      ! Get access to the PW grid pool
      CALL pw_env_get(pw_env, &
                      auxbas_pw_pool=auxbas_pw_pool, &
                      pw_pools=pw_pools, &
                      poisson_env=poisson_env)

      IF (do_nle) THEN
         IF (.NOT. ALL(cell_periodic == 1)) THEN
            CPABORT("SCCS METHOD NONLINEAR_ELECTROLYTE requires CELL PERIODIC XYZ.")
         END IF
         IF (.NOT. ALL(poisson_env%parameters%periodic == 1)) THEN
            CPABORT("SCCS METHOD NONLINEAR_ELECTROLYTE requires POISSON PERIODIC XYZ.")
         END IF
      END IF

      IF (do_nle .AND. (.NOT. sccs_control%nle_force_strict) .AND. (scf_env%iter_count == 1)) THEN
         sccs_control%nle_scf_strict = .FALSE.
         sccs_control%nle_scf_near_count = 0
         sccs_control%nle_scf_last_iter = scf_env%iter_count
      END IF

      CALL pw_zero(v_sccs)

      ! Calculate no SCCS contribution, if the requested SCF convergence threshold is not reached yet
      IF (.NOT. sccs_control%sccs_activated) THEN
         IF ((sccs_control%eps_scf > 0.0_dp) .AND. (.NOT. sccs_control%nle_force_strict)) THEN
            IF ((scf_env%iter_delta > sccs_control%eps_scf) .OR. &
                ((qs_env%scf_env%outer_scf%iter_count == 0) .AND. &
                 (qs_env%scf_env%iter_count <= 1))) THEN
               IF (calculate_stress_tensor) THEN
                  ! Request also the calculation of the stress tensor contribution
                  CALL pw_poisson_solve(poisson_env=poisson_env, &
                                        density=rho_tot_gspace, &
                                        ehartree=energy%hartree, &
                                        vhartree=v_hartree_gspace, &
                                        h_stress=h_stress)
               ELSE
                  CALL pw_poisson_solve(poisson_env=poisson_env, &
                                        density=rho_tot_gspace, &
                                        ehartree=energy%hartree, &
                                        vhartree=v_hartree_gspace)
               END IF
               energy%sccs_pol = 0.0_dp
               energy%sccs_cav = 0.0_dp
               energy%sccs_dis = 0.0_dp
               energy%sccs_rep = 0.0_dp
               energy%sccs_sol = 0.0_dp
               energy%sccs_hartree = energy%hartree
               CALL timestop(handle)
               RETURN
            END IF
         END IF
         sccs_control%sccs_activated = .TRUE.
      END IF

      nspin = dft_control%nspins

      ! Manage print output control
      logger => cp_get_default_logger()
      print_level = logger%iter_info%print_level
      print_path = "DFT%PRINT%SCCS"
      should_output = (BTEST(cp_print_key_should_output(logger%iter_info, input, &
                                                        TRIM(print_path)), cp_p_file))
      output_unit = cp_print_key_unit_nr(logger, input, TRIM(print_path), &
                                         extension=".sccs", &
                                         ignore_should_output=should_output, &
                                         log_filename=.FALSE.)

      ! Get rho
      CALL qs_rho_get(rho_struct=rho, &
                      rho_r=rho_pw_r, &
                      rho_r_sccs=rho_pw_r_sccs)

      ! Retrieve the last rho_iter from the previous SCCS cycle if available
      CPASSERT(ASSOCIATED(rho_pw_r_sccs))

      ! Retrieve the total electronic density in r-space
      BLOCK
         TYPE(pw_r3d_rs_type) :: rho_elec
         CALL auxbas_pw_pool%create_pw(rho_elec)

         ! Retrieve grid parameters
         ngpts = rho_elec%pw_grid%ngpts
         dvol = rho_elec%pw_grid%dvol
         cell_volume = rho_elec%pw_grid%vol
         abc(1:3) = REAL(rho_elec%pw_grid%npts(1:3), KIND=dp)*rho_elec%pw_grid%dr(1:3)
         lb(1:3) = rho_elec%pw_grid%bounds_local(1, 1:3)
         ub(1:3) = rho_elec%pw_grid%bounds_local(2, 1:3)

         CALL pw_copy(rho_pw_r(1), rho_elec)
         DO ispin = 2, nspin
            CALL pw_axpy(rho_pw_r(ispin), rho_elec)
         END DO
         tot_rho_elec = pw_integrate_function(rho_elec)

         ! Calculate the dielectric (smoothed) function of rho_elec in r-space
         CALL auxbas_pw_pool%create_pw(eps_elec)
         CALL auxbas_pw_pool%create_pw(deps_elec)
         IF (sccs_control%method_id == sccs_saa_andreussi) THEN
            CALL auxbas_pw_pool%create_pw(s)
            CALL auxbas_pw_pool%create_pw(d_s_rhoel)
            CALL auxbas_pw_pool%create_pw(u)
            CALL auxbas_pw_pool%create_pw(ff)
            CALL auxbas_pw_pool%create_pw(t)
            CALL auxbas_pw_pool%create_pw(dtf)
            f0 = sccs_control%f0
            delta_eta = sccs_control%delta_eta
            alpha_zeta = sccs_control%alpha_zeta
            R_solv = sccs_control%R_solv
            delta_zeta = sccs_control%delta_zeta
         END IF

         ! Relative permittivity or dielectric constant of the solvent (medium)
         epsilon_solvent = sccs_control%epsilon_solvent
         IF (do_nle) THEN
            CALL auxbas_pw_pool%create_pw(nle_ssolv)
            CALL auxbas_pw_pool%create_pw(nle_svdw)
            CALL auxbas_pw_pool%create_pw(nle_sion)
            CALL auxbas_pw_pool%create_pw(nle_sdiel)
            CALL auxbas_pw_pool%create_pw(nle_sdiel_eff)
            CALL auxbas_pw_pool%create_pw(nle_scav)
            CALL auxbas_pw_pool%create_pw(nle_xvdw)
            CALL auxbas_pw_pool%create_pw(nle_xsolv)
            CALL auxbas_pw_pool%create_pw(nle_xion)
            CALL auxbas_pw_pool%create_pw(nle_xdiel)
            CALL auxbas_pw_pool%create_pw(nle_xcav)
            CALL auxbas_pw_pool%create_pw(nle_lambda_diel)
            CALL auxbas_pw_pool%create_pw(nle_lambda_ion)
            CALL auxbas_pw_pool%create_pw(nle_lambda_cav)
            CALL auxbas_pw_pool%create_pw(nle_rho_bound)
            CALL auxbas_pw_pool%create_pw(nle_rho_ion)
            CALL auxbas_pw_pool%create_pw(nle_kernel_b)
            CALL auxbas_pw_pool%create_pw(nle_kernel_cav)
            CALL auxbas_pw_pool%create_pw(nle_kernel_diel)
            CALL auxbas_pw_pool%create_pw(nle_kernel_ion)
            CALL auxbas_pw_pool%create_pw(nle_kernel_solv)
            DO i = 1, 3
               CALL auxbas_pw_pool%create_pw(nle_eloc(i))
            END DO
            CALL pw_zero(nle_lambda_diel)
            CALL pw_zero(nle_rho_ion)
            CALL pw_zero(nle_rho_bound)
            CALL pw_zero(nle_lambda_ion)
            CALL pw_zero(nle_lambda_cav)
            CALL sccs_nle_prepare_gaussian_kernel(nle_kernel_b, &
                                                  MERGE(sccs_control%nle_r_b, sccs_control%nle_a_k, &
                                                        sccs_control%nle_r_b > 0.0_dp), &
                                                  auxbas_pw_pool)
            IF (sccs_control%nle_r_cav > 0.0_dp) &
               CALL sccs_nle_prepare_exp_kernel(nle_kernel_cav, sccs_control%nle_r_cav, &
                                                sccs_control, auxbas_pw_pool)
            IF (sccs_control%nle_r_diel > 0.0_dp) &
               CALL sccs_nle_prepare_exp_kernel(nle_kernel_diel, sccs_control%nle_r_diel, &
                                                sccs_control, auxbas_pw_pool)
            IF (sccs_control%nle_r_ion > 0.0_dp) &
               CALL sccs_nle_prepare_exp_kernel(nle_kernel_ion, sccs_control%nle_r_ion, &
                                                sccs_control, auxbas_pw_pool)
            IF (sccs_control%nle_r_solv > 0.0_dp) &
               CALL sccs_nle_prepare_exp_kernel(nle_kernel_solv, sccs_control%nle_r_solv, &
                                                sccs_control, auxbas_pw_pool)
            CALL sccs_nle_build_shapes(rho_elec, eps_elec, deps_elec, nle_svdw, nle_ssolv, &
                                       nle_sion, nle_sdiel, nle_sdiel_eff, nle_scav, nle_xvdw, &
                                       nle_xsolv, nle_xion, nle_xdiel, nle_xcav, &
                                       nle_kernel_solv, nle_kernel_ion, nle_kernel_diel, &
                                       nle_kernel_cav, sccs_control, auxbas_pw_pool, lb, ub)
         ELSE
            SELECT CASE (sccs_control%method_id)
            CASE (sccs_andreussi)
               CALL andreussi(rho_elec, eps_elec, deps_elec, epsilon_solvent, sccs_control%rho_max, &
                              sccs_control%rho_min)
            CASE (sccs_fattebert_gygi)
               CALL fattebert_gygi(rho_elec, eps_elec, deps_elec, epsilon_solvent, sccs_control%beta, &
                                   sccs_control%rho_zero)
            CASE (sccs_saa_andreussi)
               CALL sa_andreussi(rho_elec, eps_elec, deps_elec, epsilon_solvent, sccs_control%rho_max, &
                                 sccs_control%rho_min, s, d_s_rhoel)
            CASE DEFAULT
               CPABORT("Invalid method specified for SCCS model")
            END SELECT
         END IF

         IF (sccs_control%method_id /= sccs_saa_andreussi) THEN
            ! Optional output of the dielectric function in cube file format
            filename = "DIELECTRIC_FUNCTION"
            cube_path = TRIM(print_path)//"%"//TRIM(filename)
            IF (BTEST(cp_print_key_should_output(logger%iter_info, input, TRIM(cube_path)), &
                      cp_p_file)) THEN
               append_cube = section_get_lval(input, TRIM(cube_path)//"%APPEND")
               my_pos_cube = "REWIND"
               IF (append_cube) my_pos_cube = "APPEND"
               mpi_io = .TRUE.
               cube_unit = cp_print_key_unit_nr(logger, input, TRIM(cube_path), &
                                                extension=".cube", middle_name=TRIM(filename), &
                                                file_position=my_pos_cube, log_filename=.FALSE., &
                                                mpi_io=mpi_io, fout=mpi_filename)
               IF (output_unit > 0) THEN
                  IF (.NOT. mpi_io) THEN
                     INQUIRE (UNIT=cube_unit, NAME=filename)
                  ELSE
                     filename = mpi_filename
                  END IF
                  WRITE (UNIT=output_unit, FMT="(T3,A)") &
                     "SCCS| The dielectric function is written in cube file format to the file:", &
                     "SCCS| "//TRIM(filename)
               END IF
               CALL cp_pw_to_cube(eps_elec, cube_unit, TRIM(filename), particles=particles, &
                                  stride=section_get_ivals(input, TRIM(cube_path)//"%STRIDE"), &
                                  mpi_io=mpi_io)
               CALL cp_print_key_finished_output(cube_unit, logger, input, TRIM(cube_path), mpi_io=mpi_io)
            END IF
         END IF

         IF (do_nle) THEN
            CALL sccs_print_cube(nle_svdw, logger, input, particles, zcharge, output_unit, print_path, &
                                 "NLE_DENSITY_SHAPE", &
                                 "SCCS| The nonlinear electrolyte density cavity shape is written to:")
            CALL sccs_print_cube(nle_ssolv, logger, input, particles, zcharge, output_unit, print_path, &
                                 "NLE_SOLVENT_FUNCTION", &
                                 "SCCS| The nonlinear electrolyte solvent accessibility is written to:")
            CALL sccs_print_cube(nle_sion, logger, input, particles, zcharge, output_unit, print_path, &
                                 "NLE_ION_ACCESS_FUNCTION", &
                                 "SCCS| The nonlinear electrolyte ion accessibility is written to:")
            CALL sccs_print_cube(nle_sdiel, logger, input, particles, zcharge, output_unit, print_path, &
                                 "NLE_DIELECTRIC_SHAPE", &
                                 "SCCS| The nonlinear electrolyte dielectric accessibility is written to:")
            CALL sccs_print_cube(nle_sdiel_eff, logger, input, particles, zcharge, output_unit, print_path, &
                                 "NLE_EFFECTIVE_DIELECTRIC_SHAPE", &
                                 "SCCS| The effective dielectric accessibility is written to:")
            CALL sccs_print_cube(nle_scav, logger, input, particles, zcharge, output_unit, print_path, &
                                 "NLE_CAVITY_FUNCTION", &
                                 "SCCS| The nonlinear electrolyte cavitation accessibility is written to:")
         END IF

         ! Calculate the (quantum) volume and surface of the solute cavity
         cavity_surface = 0.0_dp
         cavity_volume = 0.0_dp

         IF (ABS(epsilon_solvent - 1.0_dp) > epstol) THEN

            BLOCK
               TYPE(pw_r3d_rs_type) :: theta, norm_drho_elec
               TYPE(pw_r3d_rs_type), DIMENSION(3)                        :: drho_elec
               CALL auxbas_pw_pool%create_pw(theta)
               CALL pw_zero(theta)

               ! Calculate the (quantum) volume of the solute cavity
               f = 1.0_dp/(epsilon_solvent - 1.0_dp)
!$OMP    PARALLEL DO DEFAULT(NONE) &
!$OMP                PRIVATE(i,j,k) &
!$OMP                SHARED(do_nle,epsilon_solvent,eps_elec,f,lb,theta,ub,nle_scav)
               DO k = lb(3), ub(3)
                  DO j = lb(2), ub(2)
                     DO i = lb(1), ub(1)
                        IF (do_nle) THEN
                           theta%array(i, j, k) = 1.0_dp - nle_scav%array(i, j, k)
                        ELSE
                           theta%array(i, j, k) = f*(epsilon_solvent - eps_elec%array(i, j, k))
                        END IF
                     END DO
                  END DO
               END DO
!$OMP    END PARALLEL DO
               cavity_volume = pw_integrate_function(theta)

               ! Calculate the derivative of the electronic density in r-space
               ! TODO: Could be retrieved from the QS environment
               DO i = 1, 3
                  CALL auxbas_pw_pool%create_pw(drho_elec(i))
               END DO
               IF (do_nle) THEN
                  CALL derive(nle_scav, drho_elec, sccs_derivative_fft, pw_env, input)
               ELSE
                  CALL derive(rho_elec, drho_elec, sccs_derivative_fft, pw_env, input)
               END IF

               CALL auxbas_pw_pool%create_pw(norm_drho_elec)

               ! Calculate the norm of the relevant gradient in r-space
!$OMP    PARALLEL DO DEFAULT(NONE) &
!$OMP                PRIVATE(i,j,k) &
!$OMP                SHARED(drho_elec,lb,norm_drho_elec,ub)
               DO k = lb(3), ub(3)
                  DO j = lb(2), ub(2)
                     DO i = lb(1), ub(1)
                        norm_drho_elec%array(i, j, k) = SQRT(drho_elec(1)%array(i, j, k)* &
                                                             drho_elec(1)%array(i, j, k) + &
                                                             drho_elec(2)%array(i, j, k)* &
                                                             drho_elec(2)%array(i, j, k) + &
                                                             drho_elec(3)%array(i, j, k)* &
                                                             drho_elec(3)%array(i, j, k))
                     END DO
                  END DO
               END DO
!$OMP    END PARALLEL DO

               ! Optional output of the norm of the density gradient in cube file format
               filename = "DENSITY_GRADIENT"
               cube_path = TRIM(print_path)//"%"//TRIM(filename)
               IF (BTEST(cp_print_key_should_output(logger%iter_info, input, TRIM(cube_path)), &
                         cp_p_file)) THEN
                  append_cube = section_get_lval(input, TRIM(cube_path)//"%APPEND")
                  my_pos_cube = "REWIND"
                  IF (append_cube) my_pos_cube = "APPEND"
                  mpi_io = .TRUE.
                  cube_unit = cp_print_key_unit_nr(logger, input, TRIM(cube_path), &
                                                   extension=".cube", middle_name=TRIM(filename), &
                                                   file_position=my_pos_cube, log_filename=.FALSE., &
                                                   mpi_io=mpi_io, fout=mpi_filename)
                  IF (output_unit > 0) THEN
                     IF (.NOT. mpi_io) THEN
                        INQUIRE (UNIT=cube_unit, NAME=filename)
                     ELSE
                        filename = mpi_filename
                     END IF
                     WRITE (UNIT=output_unit, FMT="(T3,A)") &
                        "SCCS| The norm of the density gradient is written in cube file format to the file:", &
                        "SCCS| "//TRIM(filename)
                  END IF
                  CALL cp_pw_to_cube(norm_drho_elec, cube_unit, "DENSITY_GRADIENT", particles=particles, &
                                     zeff=zcharge, &
                                     stride=section_get_ivals(input, TRIM(cube_path)//"%STRIDE"), &
                                     mpi_io=mpi_io)
                  CALL cp_print_key_finished_output(cube_unit, logger, input, TRIM(cube_path), mpi_io=mpi_io)
               END IF

               ! Calculate the (quantum) surface of the solute cavity
               SELECT CASE (sccs_control%method_id)
               CASE (sccs_andreussi)
                  CALL surface_andreussi(rho_elec, norm_drho_elec, theta, epsilon_solvent, &
                                         sccs_control%rho_max, sccs_control%rho_min, &
                                         sccs_control%delta_rho)
               CASE (sccs_fattebert_gygi)
                  CALL surface_fattebert_gygi(rho_elec, norm_drho_elec, theta, epsilon_solvent, &
                                              sccs_control%beta, sccs_control%rho_zero, &
                                              sccs_control%delta_rho)
               CASE (sccs_saa_andreussi)
               CASE (sccs_nonlinear_electrolyte)
!$OMP             PARALLEL DO DEFAULT(NONE) &
!$OMP                         PRIVATE(i,j,k) &
!$OMP                         SHARED(lb,norm_drho_elec,theta,ub)
                  DO k = lb(3), ub(3)
                     DO j = lb(2), ub(2)
                        DO i = lb(1), ub(1)
                           theta%array(i, j, k) = norm_drho_elec%array(i, j, k)
                        END DO
                     END DO
                  END DO
!$OMP             END PARALLEL DO
                  CALL sccs_nle_cavitation_lambda(nle_lambda_cav, drho_elec, norm_drho_elec, &
                                                  sccs_control, auxbas_pw_pool, pw_env, input, lb, ub)
               CASE DEFAULT
                  CPABORT("Invalid method specified for SCCS model")
               END SELECT
               cavity_surface = pw_integrate_function(theta)

               ! Release storage
               CALL auxbas_pw_pool%give_back_pw(theta)
               CALL auxbas_pw_pool%give_back_pw(norm_drho_elec)
               DO i = 1, 3
                  CALL auxbas_pw_pool%give_back_pw(drho_elec(i))
               END DO
            END BLOCK

         END IF ! epsilon_solvent > 1

         CALL auxbas_pw_pool%give_back_pw(rho_elec)
      END BLOCK

      BLOCK
         TYPE(pw_r3d_rs_type) :: rho_tot, phi_tot, rho_solute, rho_tot_zero
         ! Retrieve the total charge density (core + elec) of the solute in r-space
         CALL auxbas_pw_pool%create_pw(rho_solute)
         CALL pw_zero(rho_solute)
         CALL pw_transfer(rho_tot_gspace, rho_solute)
         tot_rho_solute = pw_integrate_function(rho_solute)
         electrolyte_charge = 0.0_dp
         ion_access_volume = 0.0_dp
         ion_bulk_density = 0.0_dp
         ion_capacity = 0.0_dp
         ion_charge_target = 0.0_dp
         ion_diameter = 0.0_dp
         ion_packing = 0.0_dp
         ion_phi_shift = sccs_control%nle_phi_shift
         kT = 1.0_dp
         solute_charge = REAL(dft_control%charge, KIND=dp)
         IF (do_nle) THEN
            kT = MAX(sccs_control%nle_soltemp, EPSILON(1.0_dp))
            IF (sccs_control%nle_c_molar > 0.0_dp) THEN
               ion_bulk_density = sccs_control%nle_c_molar*1000.0_dp*n_avogadro*a_bohr**3
               ion_access_volume = pw_integrate_function(nle_sion)
               ion_diameter = sccs_control%nle_d_ion
               IF (ion_diameter < 0.0_dp) ion_diameter = 2.0_dp**(5.0_dp/6.0_dp)*sccs_control%nle_r_ion
               IF (ion_diameter > 0.0_dp) THEN
                  ion_packing = 2.0_dp*ion_bulk_density*ion_diameter**3
                  ion_capacity = ABS(sccs_control%nle_zion)*ion_access_volume/ion_diameter**3
                  IF ((ion_packing < 1.0E-8_dp) .OR. (ion_packing > 1.0_dp)) THEN
                     CALL cp_abort(__LOCATION__, &
                                   "Invalid SCCS nonlinear electrolyte C_MOLAR and D_ION/R_ION for nonlinear ionic screening.")
                  END IF
               ELSE
                  ion_capacity = HUGE(1.0_dp)
               END IF
               IF (ABS(solute_charge) > ion_capacity + 1.0E-8_dp) THEN
                  WRITE (UNIT=message, FMT="(A,F12.6,A,F12.6,A)") &
                     "SCCS nonlinear electrolyte capacity is too small for charge ", &
                     ABS(solute_charge), ". Maximum compensating density charge is ", ion_capacity, "."
                  CPABORT(message)
               END IF
            END IF
         END IF

         ! Check total charge
         IF (do_nle) THEN
            IF (dft_control%charge /= 0) THEN
               IF (ion_bulk_density <= 0.0_dp) THEN
                  CPABORT("Charged SCCS METHOD NONLINEAR_ELECTROLYTE calculations require C_MOLAR > 0.")
               END IF
               IF (output_unit > 0) THEN
                  WRITE (UNIT=output_unit, FMT="(T3,A,T56,F25.12)") &
                     "SCCS| Charged solute handled by electrolyte charge", solute_charge, &
                     "SCCS| Maximum electrolyte charge capacity", ion_capacity
               END IF
            END IF
         ELSE IF (ABS(tot_rho_solute) >= 1.0E-6_dp) THEN
            IF ((poisson_env%parameters%solver /= pw_poisson_analytic) .AND. &
                (poisson_env%parameters%solver /= pw_poisson_mt)) THEN
               WRITE (UNIT=message, FMT="(A,SP,F0.6,A)") &
                  "The system (solute) has a non-negligible charge of ", -tot_rho_solute, &
                  ". It is recommended to use non-periodic boundary conditions (PERIODIC none) "// &
                  "combined with an appropriate Poisson solver (POISSON_SOLVER MT or analytic)"
               CPWARN(message)
            END IF
         END IF

         ! Reassign work storage to rho_tot_zero, because rho_elec is no longer needed
         CALL auxbas_pw_pool%create_pw(rho_tot_zero)

         IF (sccs_control%method_id /= sccs_saa_andreussi) THEN

            ! Build the initial (rho_iter = 0) total charge density (solute plus polarisation) in r-space
            ! eps_elec <- ln(eps_elec)
!$OMP PARALLEL DO DEFAULT(NONE) &
!$OMP             PRIVATE(i,j,k) &
!$OMP             SHARED(do_nle,eps_elec,lb,message,output_unit,para_env,ub) &
!$OMP             SHARED(rho_solute,rho_tot_zero)
            DO k = lb(3), ub(3)
               DO j = lb(2), ub(2)
                  DO i = lb(1), ub(1)
                     IF (eps_elec%array(i, j, k) < 1.0_dp) THEN
                        WRITE (UNIT=message, FMT="(A,ES12.3,A,3(I0,A))") &
                           "SCCS| Invalid dielectric function value ", eps_elec%array(i, j, k), &
                           " encountered at grid point (", i, ",", j, ",", k, ")"
                        CPABORT(message)
                     END IF
                     IF (do_nle) THEN
                        ! The nonlinear dielectric response is represented explicitly by its bound charge.
                        rho_tot_zero%array(i, j, k) = rho_solute%array(i, j, k)
                     ELSE
                        rho_tot_zero%array(i, j, k) = rho_solute%array(i, j, k)/eps_elec%array(i, j, k)
                     END IF
                     eps_elec%array(i, j, k) = LOG(eps_elec%array(i, j, k))
                  END DO
               END DO
            END DO
!$OMP END PARALLEL DO

         ELSE

            lxyz(:) = REAL(u%pw_grid%bounds(1, :), dp)
            uxyz(:) = REAL(u%pw_grid%bounds(2, :), dp)
!$OMP PARALLEL DO DEFAULT(NONE) &
!$OMP             PRIVATE(i,j,k,r) &
!$OMP             SHARED(lb,message,output_unit,para_env,ub,dtf,uxyz,lxyz) &
!$OMP             SHARED(u, alpha_zeta, R_solv, delta_zeta, eps_elec)
            DO k = lb(3), ub(3)
               DO j = lb(2), ub(2)
                  DO i = lb(1), ub(1)
                     IF (eps_elec%array(i, j, k) < 1.0_dp) THEN
                        WRITE (UNIT=message, FMT="(A,ES12.3,A,3(I0,A))") &
                           "SCCS| Invalid dielectric function value ", eps_elec%array(i, j, k), &
                           " encountered at grid point (", i, ",", j, ",", k, ")"
                        CPABORT(message)
                     END IF
                     IF (i > -1) r = ((REAL(uxyz(1), dp) - REAL(i, dp) + 0.5_dp)*u%pw_grid%dr(1))**2.0_dp
                     IF (i < 0) r = ((REAL(i, dp) - REAL(lxyz(1), dp) + 0.5_dp)*u%pw_grid%dr(1))**2.0_dp
                     IF (j > -1) r = ((REAL(uxyz(2), dp) - REAL(j, dp) + 0.5_dp)*u%pw_grid%dr(2))**2.0_dp + r
                     IF (j < 0) r = ((REAL(j, dp) - REAL(lxyz(2), dp) + 0.5_dp)*u%pw_grid%dr(2))**2.0_dp + r
                     IF (k > -1) r = ((REAL(uxyz(3), dp) - REAL(k, dp) + 0.5_dp)*u%pw_grid%dr(3))**2.0_dp + r
                     IF (k < 0) r = ((REAL(k, dp) - REAL(lxyz(3), dp) + 0.5_dp)*u%pw_grid%dr(3))**2.0_dp + r
                     r = SQRT(r)
                     u%array(i, j, k) = (1.0_dp/2.0_dp)*erfc((r - alpha_zeta*R_solv)/delta_zeta)
                     dtf%array(i, j, k) = 1.0_dp
                  END DO
               END DO
            END DO
!$OMP END PARALLEL DO

            CALL pw_func_u_convolution(poisson_env=poisson_env, func=dtf, convolution=ff, u=u)
            CALL pw_scale(u, 1.0_dp/ff%array(0, 0, 0))
            CALL pw_func_u_convolution(poisson_env=poisson_env, func=s, convolution=ff, u=u)

!$OMP PARALLEL DO DEFAULT(NONE) &
!$OMP             PRIVATE(i,j,k) &
!$OMP             SHARED(lb,para_env,ub, delta_eta, f0, epsilon_solvent) &
!$OMP             SHARED(ff, t, dtf, eps_elec, rho_tot_zero, rho_solute, s)
            DO k = lb(3), ub(3)
               DO j = lb(2), ub(2)
                  DO i = lb(1), ub(1)
                     t%array(i, j, k) = 0.5_dp*(1.0_dp + erf((ff%array(i, j, k) - f0)/delta_eta))
                     dtf%array(i, j, k) = (1.0_dp/SQRT(pi))*EXP(-(ff%array(i, j, k) - f0)**2.0/(delta_eta**2.0))
                     eps_elec%array(i, j, k) = EXP(LOG(epsilon_solvent)*(1 - (s%array(i, j, k) &
                                                                              + (1 - s%array(i, j, k))*t%array(i, j, k))))
                     rho_tot_zero%array(i, j, k) = rho_solute%array(i, j, k)/eps_elec%array(i, j, k)
!                     eps_elec%array(i, j, k) = LOG(eps_elec%array(i, j, k))
                  END DO
               END DO
            END DO
!$OMP END PARALLEL DO

            ! Optional output of the dielectric function in cube file format (solvent aware algorithm case)
            filename = "DIELECTRIC_FUNCTION"
            cube_path = TRIM(print_path)//"%"//TRIM(filename)
            IF (BTEST(cp_print_key_should_output(logger%iter_info, input, TRIM(cube_path)), &
                      cp_p_file)) THEN
               append_cube = section_get_lval(input, TRIM(cube_path)//"%APPEND")
               my_pos_cube = "REWIND"
               IF (append_cube) my_pos_cube = "APPEND"
               mpi_io = .TRUE.
               cube_unit = cp_print_key_unit_nr(logger, input, TRIM(cube_path), &
                                                extension=".cube", middle_name=TRIM(filename), &
                                                file_position=my_pos_cube, log_filename=.FALSE., &
                                                mpi_io=mpi_io, fout=mpi_filename)
               IF (output_unit > 0) THEN
                  IF (.NOT. mpi_io) THEN
                     INQUIRE (UNIT=cube_unit, NAME=filename)
                  ELSE
                     filename = mpi_filename
                  END IF
                  WRITE (UNIT=output_unit, FMT="(T3,A)") &
                     "SCCS| The dielectric function is written in cube file format to the file:", &
                     "SCCS| "//TRIM(filename)
               END IF
               CALL cp_pw_to_cube(eps_elec, cube_unit, TRIM(filename), particles=particles, &
                                  stride=section_get_ivals(input, TRIM(cube_path)//"%STRIDE"), &
                                  mpi_io=mpi_io)
               CALL cp_print_key_finished_output(cube_unit, logger, input, TRIM(cube_path), mpi_io=mpi_io)
            END IF

!$OMP PARALLEL DO DEFAULT(NONE) &
!$OMP             PRIVATE(i,j,k) &
!$OMP             SHARED(lb,para_env,ub, delta_eta, f0, epsilon_solvent) &
!$OMP             SHARED(ff, t, dtf, eps_elec, rho_tot_zero, rho_solute, s)
            DO k = lb(3), ub(3)
               DO j = lb(2), ub(2)
                  DO i = lb(1), ub(1)
                     eps_elec%array(i, j, k) = LOG(eps_elec%array(i, j, k))
                  END DO
               END DO
            END DO
!$OMP END PARALLEL DO

         END IF

         ! Build the derivative of LOG(eps_elec) for the original SCCS fixed-point methods.
         ! The nonlinear electrolyte Newton solver constructs its response explicitly.
         IF (.NOT. do_nle) THEN
            DO i = 1, 3
               CALL auxbas_pw_pool%create_pw(dln_eps_elec(i))
               CALL pw_zero(dln_eps_elec(i))
            END DO
            CALL derive(eps_elec, dln_eps_elec, sccs_control%derivative_method, pw_env, input)
         END IF
         IF (sccs_control%method_id /= sccs_saa_andreussi) THEN
            CALL auxbas_pw_pool%give_back_pw(eps_elec)
         END IF

         ! Print header for the SCCS cycle
         IF (should_output .AND. (output_unit > 0)) THEN
            IF (print_level > low_print_level) THEN
               WRITE (UNIT=output_unit, FMT="(T3,A,T56,F25.12)") &
                  "SCCS| Total electronic charge density ", -tot_rho_elec, &
                  "SCCS| Total charge density (solute)   ", -tot_rho_solute
               WRITE (UNIT=output_unit, FMT="(T3,A,T56,F25.3)") &
                  "SCCS| Volume of the cell           [bohr^3]", cell_volume, &
                  "SCCS|                              [angstrom^3]", &
                  cp_unit_from_cp2k(cell_volume, "angstrom^3")
               IF (ABS(epsilon_solvent - 1.0_dp) > epstol) THEN
                  WRITE (UNIT=output_unit, FMT="(T3,A,T56,F25.3)") &
                     "SCCS| Volume of the solute cavity  [bohr^3]", cavity_volume, &
                     "SCCS|                              [angstrom^3]", &
                     cp_unit_from_cp2k(cavity_volume, "angstrom^3"), &
                     "SCCS| Surface of the solute cavity [bohr^2]", cavity_surface, &
                     "SCCS|                              [angstrom^2]", &
                     cp_unit_from_cp2k(cavity_surface, "angstrom^2")
               END IF
               WRITE (UNIT=output_unit, FMT="(T3,A)") &
                  "SCCS|", &
                  "SCCS|   Step    Average residual    Maximum residual         E(Hartree) [a.u.]"
            END IF
         END IF

         ! Get storage for the derivative of the total potential (field) in r-space
         DO i = 1, 3
            CALL auxbas_pw_pool%create_pw(dphi_tot(i))
         END DO

         ! Initialise the total charge density in r-space rho_tot with rho_tot_zero + rho_iter_zero
         CALL auxbas_pw_pool%create_pw(rho_tot)
         CALL pw_copy(rho_tot_zero, rho_tot)
         CALL pw_axpy(rho_pw_r_sccs, rho_tot)

         CALL auxbas_pw_pool%create_pw(phi_tot)
         CALL pw_zero(phi_tot)

         IF (sccs_control%method_id == sccs_saa_andreussi) THEN
            CALL auxbas_pw_pool%create_pw(p_e_interface)
            CALL pw_zero(p_e_interface)
         END IF

         ! Main SCCS iteration loop
         iter = 0

         ! Use one Newton update while the electronic density is far from convergence.
         ! Once two consecutive electronic residuals are close to EPS_SCF, keep the
         ! NLE solve strict for the remainder of this SCF cycle.  Counting distinct
         ! SCF iterations avoids treating repeated OT line-search evaluations as
         ! independent evidence of electronic convergence.
         eps_sccs_current = sccs_control%eps_sccs

         IF (do_nle) THEN
            IF ((.NOT. sccs_control%nle_force_strict) .AND. &
                (.NOT. sccs_control%nle_scf_strict) .AND. &
                (scf_env%iter_count > 1) .AND. &
                (scf_env%iter_count /= sccs_control%nle_scf_last_iter)) THEN
               sccs_control%nle_scf_last_iter = scf_env%iter_count
               IF ((scf_env%iter_delta > 0.0_dp) .AND. &
                   (scf_env%iter_delta <= nle_strict_scf_factor*scf_control%eps_scf)) THEN
                  sccs_control%nle_scf_near_count = sccs_control%nle_scf_near_count + 1
               ELSE
                  sccs_control%nle_scf_near_count = 0
               END IF
               IF (sccs_control%nle_scf_near_count >= 2) THEN
                  sccs_control%nle_scf_strict = .TRUE.
                  IF (output_unit > 0) THEN
                     WRITE (UNIT=output_unit, FMT="(T3,A,ES12.4,A,ES12.4)") &
                        "SCCS| Entering strict electrolyte SCF phase at electronic residual ", &
                        scf_env%iter_delta, "; EPS_SCF ", scf_control%eps_scf
                  END IF
               END IF
            END IF

            nle_strict_solve = sccs_control%nle_force_strict .OR. sccs_control%nle_scf_strict
            nle_max_outer = sccs_control%max_iter
            IF (.NOT. nle_strict_solve) nle_max_outer = 1
            CALL sccs_nle_newton_solve(rho_solute, rho_pw_r_sccs, rho_tot, phi_tot, &
                                       nle_rho_bound, nle_rho_ion, nle_lambda_diel, nle_lambda_ion, &
                                       nle_eloc, nle_sdiel, nle_sion, nle_kernel_b, &
                                       sccs_control, poisson_env, auxbas_pw_pool, ion_bulk_density, &
                                       ion_packing, kT, solute_charge, cell_volume, eps_sccs_current, &
                                       nle_max_outer, output_unit,.NOT. nle_strict_solve, &
                                       ion_phi_shift, iter, nle_pcg_iterations, &
                                       nle_residual_rms, nle_converged)
            IF (should_output .AND. (output_unit > 0) .AND. (print_level > low_print_level)) THEN
               WRITE (UNIT=output_unit, FMT="(T3,A,I0,A,I0,A,ES12.4)") &
                  "SCCS| Electrolyte Newton solver used ", iter, " updates and ", &
                  nle_pcg_iterations, " PCG iterations; residual ", nle_residual_rms
            END IF
            IF ((.NOT. nle_converged) .AND. &
                (.NOT. ((.NOT. nle_strict_solve) .AND. (iter == 1)))) THEN
               WRITE (UNIT=message, FMT="(A,I0,A,ES12.4)") &
                  "The SCCS electrolyte Newton solver did not converge in ", iter, &
                  " updates; residual is ", nle_residual_rms
               CPWARN(message)
            END IF
         ELSE
            iter_loop: DO

               ! Increment iteration counter
               iter = iter + 1

               ! Check if the requested maximum number of SCCS iterations is reached
               IF (iter > sccs_control%max_iter) THEN
                  WRITE (UNIT=message, FMT="(A,I0,A)") &
                     "The SCCS iteration cycle did not converge in ", sccs_control%max_iter, " steps"
                  IF (output_unit > 0) THEN
                     WRITE (UNIT=output_unit, FMT="(T3,A,/,T3,A,I0,A)") &
                        "SCCS| Maximum number of SCCS iterations reached", &
                        "SCCS| Iteration cycle did not converge in ", sccs_control%max_iter, " steps"
                  END IF
                  CPWARN(message)
                  EXIT iter_loop
               END IF

               ! Calculate the current total potential and its derivative in r-space.
               IF (do_nle) THEN
                  CALL pw_transfer(rho_tot, rho_tot_gspace)
                  CALL pw_poisson_solve(poisson_env=poisson_env, &
                                        density=rho_tot_gspace, &
                                        ehartree=e_tot, &
                                        vhartree=v_hartree_gspace)
                  CALL pw_transfer(v_hartree_gspace, phi_tot)
                  CALL derive(phi_tot, dphi_tot, sccs_derivative_fft, pw_env, input)
               ELSE
                  CALL pw_poisson_solve(poisson_env=poisson_env, &
                                        density=rho_tot, &
                                        vhartree=phi_tot, &
                                        dvhartree=dphi_tot)
               END IF
               energy%sccs_hartree = 0.5_dp*pw_integral_ab(rho_solute, phi_tot)
               IF (do_nle) THEN
                  CALL sccs_nle_update_dielectric_response(nle_rho_bound, nle_lambda_diel, &
                                                           nle_eloc, phi_tot, nle_sdiel, nle_kernel_b, &
                                                           kT, sccs_control, auxbas_pw_pool, lb, ub)
               END IF
               IF (do_nle .AND. (ion_bulk_density > 0.0_dp)) THEN
                  ion_charge_target = solute_charge
                  CALL sccs_nle_update_ion_density(nle_rho_ion, phi_tot, nle_sion, sccs_control, &
                                                   ion_bulk_density, ion_packing, kT, &
                                                   ion_charge_target, ion_phi_shift, auxbas_pw_pool, lb, ub)
               END IF

               ! Update total charge density (solute plus polarisation) in r-space
               ! based on the iterated polarisation charge density
               f = 1.0_dp/fourpi
               rho_delta_avg = 0.0_dp
               rho_delta_max = 0.0_dp
!$OMP    PARALLEL DO DEFAULT(NONE) &
!$OMP                PRIVATE(i,j,k,rho_delta,rho_iter_new) &
!$OMP                SHARED(dln_eps_elec,do_nle,dphi_tot,f,lb,rho_pw_r_sccs,rho_tot) &
!$OMP                SHARED(rho_tot_zero,sccs_control,ub,nle_rho_bound,nle_rho_ion) &
!$OMP                REDUCTION(+:rho_delta_avg) &
!$OMP                REDUCTION(MAX:rho_delta_max)
               DO k = lb(3), ub(3)
                  DO j = lb(2), ub(2)
                     DO i = lb(1), ub(1)
                        IF (do_nle) THEN
                           rho_iter_new = nle_rho_bound%array(i, j, k) + nle_rho_ion%array(i, j, k)
                        ELSE
                           rho_iter_new = (dln_eps_elec(1)%array(i, j, k)*dphi_tot(1)%array(i, j, k) + &
                                           dln_eps_elec(2)%array(i, j, k)*dphi_tot(2)%array(i, j, k) + &
                                           dln_eps_elec(3)%array(i, j, k)*dphi_tot(3)%array(i, j, k))*f
                        END IF
                        rho_delta = ABS(rho_iter_new - rho_pw_r_sccs%array(i, j, k))
                        rho_iter_new = rho_pw_r_sccs%array(i, j, k) + &
                                       sccs_control%mixing*(rho_iter_new - rho_pw_r_sccs%array(i, j, k))
                        ! Keep the 2026.2 SCCS convergence measure for the existing methods. The
                        ! nonlinear electrolyte reference implementation measures the unmixed update.
                        IF (.NOT. do_nle) rho_delta = ABS(rho_iter_new - rho_pw_r_sccs%array(i, j, k))
                        rho_delta_max = MAX(rho_delta, rho_delta_max)
                        rho_delta_avg = rho_delta_avg + rho_delta
                        rho_tot%array(i, j, k) = rho_tot_zero%array(i, j, k) + rho_iter_new
                        rho_pw_r_sccs%array(i, j, k) = rho_iter_new
                     END DO
                  END DO
               END DO
!$OMP    END PARALLEL DO

               IF (do_nle) THEN
                  response_charge = pw_integrate_function(rho_pw_r_sccs)
                  charge_correction = (solute_charge - response_charge)/cell_volume
!$OMP       PARALLEL DO DEFAULT(NONE) PRIVATE(i,j,k) &
!$OMP                   SHARED(charge_correction,lb,rho_pw_r_sccs,rho_tot,ub)
                  DO k = lb(3), ub(3)
                     DO j = lb(2), ub(2)
                        DO i = lb(1), ub(1)
                           rho_pw_r_sccs%array(i, j, k) = rho_pw_r_sccs%array(i, j, k) + charge_correction
                           rho_tot%array(i, j, k) = rho_tot%array(i, j, k) + charge_correction
                        END DO
                     END DO
                  END DO
!$OMP       END PARALLEL DO
               END IF

               CALL para_env%sum(rho_delta_avg)
               rho_delta_avg = rho_delta_avg/REAL(ngpts, KIND=dp)
               CALL para_env%max(rho_delta_max)

               IF (should_output .AND. (output_unit > 0)) THEN
                  IF (print_level > low_print_level) THEN
                     IF ((ABS(rho_delta_avg) < 1.0E-8_dp) .OR. &
                         (ABS(rho_delta_avg) >= 1.0E5_dp)) THEN
                        WRITE (UNIT=output_unit, FMT="(T3,A,I6,4X,ES16.4,4X,ES16.4,1X,F25.12)") &
                           "SCCS| ", iter, rho_delta_avg, rho_delta_max, energy%sccs_hartree
                     ELSE
                        WRITE (UNIT=output_unit, FMT="(T3,A,I6,4X,F16.8,4X,F16.8,1X,F25.12)") &
                           "SCCS| ", iter, rho_delta_avg, rho_delta_max, energy%sccs_hartree
                     END IF
                  END IF
               END IF

               ! Check if the SCCS iteration cycle is converged to the requested tolerance
               IF (rho_delta_max <= eps_sccs_current) THEN
                  IF (should_output .AND. (output_unit > 0)) THEN
                     WRITE (UNIT=output_unit, FMT="(T3,A,I0,A)") &
                        "SCCS| Iteration cycle converged in ", iter, " steps"
                  END IF
                  EXIT iter_loop
               END IF

            END DO iter_loop
         END IF

         IF (sccs_control%method_id == sccs_saa_andreussi) THEN
            CALL auxbas_pw_pool%create_pw(tem_convolution)
            CALL auxbas_pw_pool%create_pw(tem_func)

!$OMP PARALLEL DO DEFAULT(NONE) &
!$OMP                PRIVATE(i,j,k,dphi2) &
!$OMP                SHARED(lb,ub, tem_func, s, dtf) &
!$OMP                SHARED(dphi_tot, p_e_interface, eps_elec, epsilon_solvent)
            DO k = lb(3), ub(3)
               DO j = lb(2), ub(2)
                  DO i = lb(1), ub(1)
                     p_e_interface%array(i, j, k) = 1.0_dp/(4.0_dp*twopi)*(EXP(eps_elec%array(i, j, k))*LOG(epsilon_solvent))
                     dphi2 = dphi_tot(1)%array(i, j, k)*dphi_tot(1)%array(i, j, k) + &
                             dphi_tot(2)%array(i, j, k)*dphi_tot(2)%array(i, j, k) + &
                             dphi_tot(3)%array(i, j, k)*dphi_tot(3)%array(i, j, k)
                     p_e_interface%array(i, j, k) = p_e_interface%array(i, j, k)*dphi2
                     tem_func%array(i, j, k) = (1.0_dp - s%array(i, j, k))*p_e_interface%array(i, j, k)*dtf%array(i, j, k)
                  END DO
               END DO
            END DO
!$OMP END PARALLEL DO

            CALL pw_func_u_convolution(poisson_env=poisson_env, func=tem_func, convolution=tem_convolution, u=u)

!$OMP PARALLEL DO DEFAULT(NONE) &
!$OMP                PRIVATE(i,j,k) &
!$OMP                SHARED(lb,ub) &
!$OMP                SHARED(p_e_interface, t, tem_convolution, phi_tot, d_s_rhoel)
            DO k = lb(3), ub(3)
               DO j = lb(2), ub(2)
                  DO i = lb(1), ub(1)
                     p_e_interface%array(i, j, k) = p_e_interface%array(i, j, k)*(1.0_dp - t%array(i, j, k)) &
                                                    + tem_convolution%array(i, j, k)
                  END DO
               END DO
            END DO
!$OMP END PARALLEL DO

         END IF

         ! Release work storage which is no longer needed
         CALL auxbas_pw_pool%give_back_pw(rho_tot_zero)
         IF (.NOT. do_nle) THEN
            DO i = 1, 3
               CALL auxbas_pw_pool%give_back_pw(dln_eps_elec(i))
            END DO
         END IF

         ! Optional output of the total charge density in cube file format
         filename = "TOTAL_CHARGE_DENSITY"
         cube_path = TRIM(print_path)//"%"//TRIM(filename)
         IF (BTEST(cp_print_key_should_output(logger%iter_info, input, TRIM(cube_path)), cp_p_file)) THEN
            append_cube = section_get_lval(input, TRIM(cube_path)//"%APPEND")
            my_pos_cube = "REWIND"
            IF (append_cube) my_pos_cube = "APPEND"
            mpi_io = .TRUE.
            cube_unit = cp_print_key_unit_nr(logger, input, TRIM(cube_path), &
                                             extension=".cube", middle_name=TRIM(filename), &
                                             file_position=my_pos_cube, log_filename=.FALSE., &
                                             mpi_io=mpi_io, fout=mpi_filename)
            IF (output_unit > 0) THEN
               IF (.NOT. mpi_io) THEN
                  INQUIRE (UNIT=cube_unit, NAME=filename)
               ELSE
                  filename = mpi_filename
               END IF
               WRITE (UNIT=output_unit, FMT="(T3,A)") &
                  "SCCS| The total SCCS charge density is written in cube file format to the file:", &
                  "SCCS| "//TRIM(filename)
            END IF
            CALL cp_pw_to_cube(rho_tot, cube_unit, "TOTAL_CHARGE_DENSITY", particles=particles, &
                               zeff=zcharge, &
                               stride=section_get_ivals(input, TRIM(cube_path)//"%STRIDE"), mpi_io=mpi_io)
            CALL cp_print_key_finished_output(cube_unit, logger, input, TRIM(cube_path), mpi_io=mpi_io)
         END IF

         ! Calculate the total SCCS Hartree energy, potential, and its
         ! derivatives of the solute and the implicit solvent
         CALL pw_transfer(rho_tot, rho_tot_gspace)
         IF (calculate_stress_tensor) THEN
            ! Request also the calculation of the stress tensor contribution
            CALL pw_poisson_solve(poisson_env=poisson_env, &
                                  density=rho_tot_gspace, &
                                  ehartree=e_tot, &
                                  vhartree=v_hartree_gspace, &
                                  dvhartree=dphi_tot, &
                                  h_stress=h_stress)
         ELSE IF (do_nle) THEN
            CALL pw_poisson_solve(poisson_env=poisson_env, &
                                  density=rho_tot_gspace, &
                                  ehartree=e_tot, &
                                  vhartree=v_hartree_gspace)
         ELSE
            CALL pw_poisson_solve(poisson_env=poisson_env, &
                                  density=rho_tot_gspace, &
                                  ehartree=e_tot, &
                                  vhartree=v_hartree_gspace, &
                                  dvhartree=dphi_tot)
         END IF
         CALL pw_transfer(v_hartree_gspace, phi_tot)
         IF (do_nle) CALL derive(phi_tot, dphi_tot, sccs_derivative_fft, pw_env, input)

         ! The SCCS loop evaluates the nonlinear response before mixing it into
         ! rho_tot.  The Poisson solve above uses the mixed response, so refresh
         ! the nonlinear fields once at that final potential.  Otherwise the
         ! printed response fields and the free-energy densities lag one SCCS
         ! iteration behind the total potential, which is visible when a loose
         ! SCCS tolerance is used.
         IF (do_nle) THEN
            CALL sccs_nle_update_dielectric_response(nle_rho_bound, nle_lambda_diel, &
                                                     nle_eloc, phi_tot, nle_sdiel, nle_kernel_b, &
                                                     kT, sccs_control, auxbas_pw_pool, lb, ub)
            IF (ion_bulk_density > 0.0_dp) THEN
               ion_charge_target = solute_charge
               CALL sccs_nle_update_ion_density(nle_rho_ion, phi_tot, nle_sion, sccs_control, &
                                                ion_bulk_density, ion_packing, kT, &
                                                ion_charge_target, ion_phi_shift, auxbas_pw_pool, lb, ub)
            END IF
         END IF
         energy%sccs_hartree = 0.5_dp*pw_integral_ab(rho_solute, phi_tot)

         ! nonlinear electrolyte evaluates the ionic grand-potential density at the total
         ! electrostatic potential, including its charge-neutrality shift.
         IF (do_nle .AND. (ion_bulk_density > 0.0_dp)) THEN
!$OMP       PARALLEL DO DEFAULT(NONE) &
!$OMP                   PRIVATE(beta_phi,i,ion_denom,ion_free_density,j,k) &
!$OMP                   SHARED(ion_bulk_density,ion_packing,ion_phi_shift,kT,lb,phi_tot) &
!$OMP                   SHARED(sccs_control,ub,nle_lambda_ion)
            DO k = lb(3), ub(3)
               DO j = lb(2), ub(2)
                  DO i = lb(1), ub(1)
                     beta_phi = sccs_control%nle_zion*(phi_tot%array(i, j, k) + ion_phi_shift)/kT
                     beta_phi = MAX(-40.0_dp, MIN(40.0_dp, beta_phi))
                     IF (ion_packing > 0.0_dp) THEN
                        ion_denom = 1.0_dp + ion_packing*(COSH(beta_phi) - 1.0_dp)
                        ion_free_density = -2.0_dp*kT*ion_bulk_density*LOG(ion_denom)/ion_packing
                     ELSE
                        ion_free_density = -2.0_dp*kT*ion_bulk_density*(COSH(beta_phi) - 1.0_dp)
                     END IF
                     nle_lambda_ion%array(i, j, k) = ion_free_density
                  END DO
               END DO
            END DO
!$OMP       END PARALLEL DO
         END IF

         ! Calculate the Hartree energy and potential of the solute only
         BLOCK
            TYPE(pw_r3d_rs_type) :: phi_solute
            CALL auxbas_pw_pool%create_pw(phi_solute)
            CALL pw_zero(phi_solute)
            ! Keep v_hartree_gspace as the total solute-plus-response potential
            ! returned to the KS matrix.  The solute-only reference potential is
            ! needed only for the polarisation energy and must use separate storage.
            CALL pw_poisson_solve(poisson_env=poisson_env, &
                                  density=rho_solute, &
                                  ehartree=energy%hartree, &
                                  vhartree=phi_solute)

            ! Calculate the polarisation potential (store it in phi_tot)
            ! phi_pol = phi_tot - phi_solute
            CALL pw_axpy(phi_solute, phi_tot, alpha=-1.0_dp)
            CALL auxbas_pw_pool%give_back_pw(phi_solute)
         END BLOCK

         ! Optional output of the SCCS polarisation potential in cube file format
         filename = "POLARISATION_POTENTIAL"
         cube_path = TRIM(print_path)//"%"//TRIM(filename)
         IF (BTEST(cp_print_key_should_output(logger%iter_info, input, TRIM(cube_path)), &
                   cp_p_file)) THEN
            append_cube = section_get_lval(input, TRIM(cube_path)//"%APPEND")
            my_pos_cube = "REWIND"
            IF (append_cube) my_pos_cube = "APPEND"
            mpi_io = .TRUE.
            cube_unit = cp_print_key_unit_nr(logger, input, TRIM(cube_path), &
                                             extension=".cube", middle_name=TRIM(filename), &
                                             file_position=my_pos_cube, log_filename=.FALSE., &
                                             mpi_io=mpi_io, fout=mpi_filename)
            IF (output_unit > 0) THEN
               IF (.NOT. mpi_io) THEN
                  INQUIRE (UNIT=cube_unit, NAME=filename)
               ELSE
                  filename = mpi_filename
               END IF
               WRITE (UNIT=output_unit, FMT="(T3,A)") &
                  "SCCS| The SCCS polarisation potential is written in cube file format to the file:", &
                  "SCCS| "//TRIM(filename)
            END IF
            ! Match the potential used by the nonlinear ionic functional.
            phi_tot%array = phi_tot%array + ion_phi_shift
            CALL cp_pw_to_cube(phi_tot, cube_unit, "POLARISATION_POTENTIAL (INCLUDING G=0 SHIFT)", particles=particles, &
                               zeff=zcharge, &
                               stride=section_get_ivals(input, TRIM(cube_path)//"%STRIDE"), mpi_io=mpi_io)
            phi_tot%array = phi_tot%array - ion_phi_shift
            CALL cp_print_key_finished_output(cube_unit, logger, input, TRIM(cube_path), mpi_io=mpi_io)
         END IF

         ! Calculate the polarisation charge (store it in rho_tot)
         ! rho_pol = rho_tot - rho_solute
         CALL pw_axpy(rho_solute, rho_tot, alpha=-1.0_dp)
         polarisation_charge = pw_integrate_function(rho_tot)
         IF (do_nle) electrolyte_charge = pw_integrate_function(nle_rho_ion)

         ! Optional output of the SCCS polarisation charge in cube file format
         filename = "POLARISATION_CHARGE_DENSITY"
         cube_path = TRIM(print_path)//"%"//TRIM(filename)
         IF (BTEST(cp_print_key_should_output(logger%iter_info, input, TRIM(cube_path)), &
                   cp_p_file)) THEN
            append_cube = section_get_lval(input, TRIM(cube_path)//"%APPEND")
            my_pos_cube = "REWIND"
            IF (append_cube) my_pos_cube = "APPEND"
            mpi_io = .TRUE.
            cube_unit = cp_print_key_unit_nr(logger, input, TRIM(cube_path), &
                                             extension=".cube", middle_name=TRIM(filename), &
                                             file_position=my_pos_cube, log_filename=.FALSE., &
                                             mpi_io=mpi_io, fout=mpi_filename)
            IF (output_unit > 0) THEN
               IF (.NOT. mpi_io) THEN
                  INQUIRE (UNIT=cube_unit, NAME=filename)
               ELSE
                  filename = mpi_filename
               END IF
               WRITE (UNIT=output_unit, FMT="(T3,A)") &
                  "SCCS| The SCCS polarisation charge density is written in cube file format to the file:", &
                  "SCCS| "//TRIM(filename)
            END IF
            CALL cp_pw_to_cube(rho_tot, cube_unit, "POLARISATION_CHARGE_DENSITY", particles=particles, &
                               zeff=zcharge, &
                               stride=section_get_ivals(input, TRIM(cube_path)//"%STRIDE"), mpi_io=mpi_io)
            CALL cp_print_key_finished_output(cube_unit, logger, input, TRIM(cube_path), mpi_io=mpi_io)
         END IF
         IF (do_nle) THEN
            CALL sccs_print_cube(nle_rho_bound, logger, input, particles, zcharge, output_unit, print_path, &
                                 "NLE_BOUND_CHARGE_DENSITY", &
                                 "SCCS| The nonlinear electrolyte bound charge density is written to:")
            CALL sccs_print_cube(nle_rho_ion, logger, input, particles, zcharge, output_unit, print_path, &
                                 "NLE_ION_CHARGE_DENSITY", &
                                 "SCCS| The nonlinear electrolyte mobile ionic charge density is written to:")
         END IF

         ! Calculate SCCS polarisation energy.  The nonlinear electrolyte
         ! functional is not the linear-response half interaction energy.
         e_nle_shift = 0.0_dp
         e_nle_field = 0.0_dp
         e_nle_diel = 0.0_dp
         e_nle_ion = 0.0_dp
         IF (do_nle) THEN
            ! solute_charge is the neutralising response-density target.  The
            ! solute charge in the Poisson-density convention has the opposite sign.
            e_nle_shift = -solute_charge*ion_phi_shift
            e_nle_field = -0.5_dp*pw_integral_ab(rho_tot, phi_tot)
            e_nle_diel = pw_integral_ab(nle_sdiel, nle_lambda_diel)
            IF (ion_bulk_density > 0.0_dp) THEN
               e_nle_ion = pw_integral_ab(nle_sion, nle_lambda_ion)
            END IF
            energy%sccs_pol = e_nle_shift + e_nle_field + e_nle_diel + e_nle_ion
         ELSE
            energy%sccs_pol = 0.5_dp*pw_integral_ab(rho_solute, phi_tot)
         END IF
         CALL auxbas_pw_pool%give_back_pw(rho_solute)
         CALL auxbas_pw_pool%give_back_pw(phi_tot)
         CALL auxbas_pw_pool%give_back_pw(rho_tot)
         IF (sccs_control%method_id == sccs_saa_andreussi) THEN
            CALL auxbas_pw_pool%give_back_pw(s)
            CALL auxbas_pw_pool%give_back_pw(u)
            CALL auxbas_pw_pool%give_back_pw(ff)
            CALL auxbas_pw_pool%give_back_pw(t)
            CALL auxbas_pw_pool%give_back_pw(dtf)
            CALL auxbas_pw_pool%give_back_pw(tem_convolution)
            CALL auxbas_pw_pool%give_back_pw(tem_func)
            CALL auxbas_pw_pool%give_back_pw(eps_elec)
         END IF
      END BLOCK

      ! Preserve the electrolyte reference potential for post-SCF output.
      sccs_control%nle_phi_shift = ion_phi_shift

      ! Calculate additional solvation terms
      energy%sccs_cav = sccs_control%gamma_solvent*cavity_surface
      energy%sccs_dis = sccs_control%beta_solvent*cavity_volume
      energy%sccs_rep = sccs_control%alpha_solvent*cavity_surface
      ! Calculate solvation free energy: \delta G^el + (alpha + gamma)*S + beta*V
      energy%sccs_sol = energy%sccs_pol + energy%sccs_rep + energy%sccs_cav + energy%sccs_dis

      IF (should_output .AND. (output_unit > 0)) THEN
         WRITE (UNIT=output_unit, FMT="(T3,A)") &
            "SCCS|"
         WRITE (UNIT=output_unit, FMT="(T3,A,T56,F25.12)") &
            "SCCS| Polarisation charge", polarisation_charge
         IF (do_nle) THEN
            WRITE (UNIT=output_unit, FMT="(T3,A,T56,F25.12)") &
               "SCCS| Ionic response charge", electrolyte_charge, &
               "SCCS| NLE G=0 potential shift [Hartree]", sccs_control%nle_phi_shift
            WRITE (UNIT=output_unit, FMT="(T3,A,T56,F25.12)") &
               "SCCS| NLE G=0 energy [Hartree]", e_nle_shift, &
               "SCCS| NLE field energy [Hartree]", e_nle_field, &
               "SCCS| NLE dielectric free energy [Hartree]", e_nle_diel, &
               "SCCS| NLE ionic free energy [Hartree]", e_nle_ion
         END IF
         !MK   "SCCS| Total interaction energy [a.u.]", e_tot
         WRITE (UNIT=output_unit, FMT="(T3,A)") &
            "SCCS|"
         CALL print_sccs_results(energy, sccs_control, output_unit)
      END IF

      ! Calculate SCCS contribution to the Kohn-Sham potential
      IF (do_nle) THEN
         CALL sccs_nle_add_dielectric_potential(v_sccs, nle_lambda_diel, nle_xvdw, nle_xsolv, &
                                                nle_xdiel, nle_sdiel, nle_svdw, sccs_control, &
                                                nle_kernel_diel, nle_kernel_solv, auxbas_pw_pool, &
                                                dvol, lb, ub)
         IF (ion_bulk_density > 0.0_dp) THEN
            CALL sccs_nle_add_ion_potential(v_sccs, nle_lambda_ion, nle_xvdw, nle_xion, &
                                            sccs_control, nle_kernel_ion, auxbas_pw_pool, &
                                            dvol, lb, ub)
         END IF
         CALL sccs_nle_add_cavity_potential(v_sccs, nle_lambda_cav, nle_xvdw, nle_xsolv, &
                                            nle_xcav, sccs_control, nle_kernel_cav, nle_kernel_solv, &
                                            auxbas_pw_pool, dvol, lb, ub)

         ! Remove the real-space integration weight before printing the physical potential.
         CALL pw_copy(v_sccs, deps_elec)
         CALL pw_scale(deps_elec, 1.0_dp/dvol)
         CALL sccs_print_cube(deps_elec, logger, input, particles, zcharge, output_unit, print_path, &
                              "NLE_SHAPE_POTENTIAL", &
                              "SCCS| The nonlinear electrolyte cavity shape potential is written to:")
      ELSE IF (sccs_control%method_id /= sccs_saa_andreussi) THEN
         f = -0.5_dp*dvol/fourpi
!$OMP PARALLEL DO DEFAULT(NONE) &
!$OMP             PRIVATE(dphi2,i,j,k) &
!$OMP             SHARED(f,deps_elec,dphi_tot) &
!$OMP             SHARED(lb,ub,v_sccs)
         DO k = lb(3), ub(3)
            DO j = lb(2), ub(2)
               DO i = lb(1), ub(1)
                  dphi2 = dphi_tot(1)%array(i, j, k)*dphi_tot(1)%array(i, j, k) + &
                          dphi_tot(2)%array(i, j, k)*dphi_tot(2)%array(i, j, k) + &
                          dphi_tot(3)%array(i, j, k)*dphi_tot(3)%array(i, j, k)
                  v_sccs%array(i, j, k) = v_sccs%array(i, j, k) + f*deps_elec%array(i, j, k)*dphi2
               END DO
            END DO
         END DO
!$OMP END PARALLEL DO
      ELSE
!$OMP PARALLEL DO DEFAULT(NONE) &
!$OMP                PRIVATE(i,j,k) &
!$OMP                SHARED(lb,ub, dvol) &
!$OMP                SHARED(p_E_interface, v_sccs, d_s_rhoel)
         DO k = lb(3), ub(3)
            DO j = lb(2), ub(2)
               DO i = lb(1), ub(1)
                  v_sccs%array(i, j, k) = dvol*d_s_rhoel%array(i, j, k)*p_E_interface%array(i, j, k)
               END DO
            END DO
         END DO
!$OMP END PARALLEL DO
         CALL auxbas_pw_pool%give_back_pw(d_s_rhoel)
         CALL auxbas_pw_pool%give_back_pw(p_E_interface)
      END IF

      CALL auxbas_pw_pool%give_back_pw(deps_elec)
      IF (do_nle) THEN
         CALL auxbas_pw_pool%give_back_pw(nle_svdw)
         CALL auxbas_pw_pool%give_back_pw(nle_ssolv)
         CALL auxbas_pw_pool%give_back_pw(nle_sion)
         CALL auxbas_pw_pool%give_back_pw(nle_sdiel)
         CALL auxbas_pw_pool%give_back_pw(nle_sdiel_eff)
         CALL auxbas_pw_pool%give_back_pw(nle_scav)
         CALL auxbas_pw_pool%give_back_pw(nle_xvdw)
         CALL auxbas_pw_pool%give_back_pw(nle_xsolv)
         CALL auxbas_pw_pool%give_back_pw(nle_xion)
         CALL auxbas_pw_pool%give_back_pw(nle_xdiel)
         CALL auxbas_pw_pool%give_back_pw(nle_xcav)
         CALL auxbas_pw_pool%give_back_pw(nle_lambda_diel)
         CALL auxbas_pw_pool%give_back_pw(nle_lambda_ion)
         CALL auxbas_pw_pool%give_back_pw(nle_lambda_cav)
         CALL auxbas_pw_pool%give_back_pw(nle_rho_bound)
         CALL auxbas_pw_pool%give_back_pw(nle_rho_ion)
         CALL auxbas_pw_pool%give_back_pw(nle_kernel_b)
         CALL auxbas_pw_pool%give_back_pw(nle_kernel_cav)
         CALL auxbas_pw_pool%give_back_pw(nle_kernel_diel)
         CALL auxbas_pw_pool%give_back_pw(nle_kernel_ion)
         CALL auxbas_pw_pool%give_back_pw(nle_kernel_solv)
         DO i = 1, 3
            CALL auxbas_pw_pool%give_back_pw(nle_eloc(i))
         END DO
      END IF
      DO i = 1, 3
         CALL auxbas_pw_pool%give_back_pw(dphi_tot(i))
      END DO
      DEALLOCATE (zcharge)

      ! Release the SCCS printout environment
      CALL cp_print_key_finished_output(output_unit, logger, input, TRIM(print_path), &
                                        ignore_should_output=should_output)

      CALL timestop(handle)

   END SUBROUTINE sccs

! **************************************************************************************************
!> \brief Print a real-space cube field controlled by a SCCS print key.
!> \param field ...
!> \param logger ...
!> \param input ...
!> \param particles ...
!> \param zcharge ...
!> \param output_unit ...
!> \param print_path ...
!> \param field_name ...
!> \param log_message ...
! **************************************************************************************************
   SUBROUTINE sccs_print_cube(field, logger, input, particles, zcharge, output_unit, print_path, &
                              field_name, log_message)

      TYPE(pw_r3d_rs_type), INTENT(IN)                   :: field
      TYPE(cp_logger_type), POINTER                      :: logger
      TYPE(section_vals_type), POINTER                   :: input
      TYPE(particle_list_type), POINTER                  :: particles
      REAL(KIND=dp), DIMENSION(:), INTENT(IN)            :: zcharge
      INTEGER, INTENT(IN)                                :: output_unit
      CHARACTER(LEN=*), INTENT(IN)                       :: print_path, field_name, log_message

      CHARACTER(LEN=default_path_length)                 :: cube_path, filename, mpi_filename, &
                                                            my_pos_cube
      INTEGER                                            :: cube_unit
      LOGICAL                                            :: append_cube, mpi_io

      filename = field_name
      cube_path = TRIM(print_path)//"%"//TRIM(field_name)
      IF (BTEST(cp_print_key_should_output(logger%iter_info, input, TRIM(cube_path)), &
                cp_p_file)) THEN
         append_cube = section_get_lval(input, TRIM(cube_path)//"%APPEND")
         my_pos_cube = "REWIND"
         IF (append_cube) my_pos_cube = "APPEND"
         mpi_io = .TRUE.
         cube_unit = cp_print_key_unit_nr(logger, input, TRIM(cube_path), &
                                          extension=".cube", middle_name=TRIM(field_name), &
                                          file_position=my_pos_cube, log_filename=.FALSE., &
                                          mpi_io=mpi_io, fout=mpi_filename)
         IF (output_unit > 0) THEN
            IF (.NOT. mpi_io) THEN
               INQUIRE (UNIT=cube_unit, NAME=filename)
            ELSE
               filename = mpi_filename
            END IF
            WRITE (UNIT=output_unit, FMT="(T3,A)") &
               TRIM(log_message), &
               "SCCS| "//TRIM(filename)
         END IF
         CALL cp_pw_to_cube(field, cube_unit, TRIM(field_name), particles=particles, zeff=zcharge, &
                            stride=section_get_ivals(input, TRIM(cube_path)//"%STRIDE"), mpi_io=mpi_io)
         CALL cp_print_key_finished_output(cube_unit, logger, input, TRIM(cube_path), mpi_io=mpi_io)
      END IF

   END SUBROUTINE sccs_print_cube

! **************************************************************************************************
!> \brief Solve the nonlinear electrolyte equation with Newton updates and a vacuum-preconditioned
!>        conjugate-gradient solution of each linearized equation.
! **************************************************************************************************
   SUBROUTINE sccs_nle_newton_solve(rho_solute, rho_response, rho_tot, phi, rho_bound, rho_ion, &
                                    lambda_diel, lambda_ion, eloc, sdiel, sion, kernel_g, &
                                    sccs_control, poisson_env, pw_pool, ion_bulk_density, &
                                    ion_packing, kT, solute_charge, cell_volume, tolerance, &
                                    max_outer, output_unit, adaptive_tolerance, phi_shift, nouter, npcgtot, &
                                    residual_rms, converged)

      TYPE(pw_r3d_rs_type), INTENT(IN)                   :: rho_solute
      TYPE(pw_r3d_rs_type), INTENT(INOUT)                :: rho_response, rho_tot, phi, rho_bound, &
                                                            rho_ion, lambda_diel, lambda_ion
      TYPE(pw_r3d_rs_type), DIMENSION(3), INTENT(INOUT)  :: eloc
      TYPE(pw_r3d_rs_type), INTENT(IN)                   :: sdiel, sion
      TYPE(pw_c1d_gs_type), INTENT(IN)                   :: kernel_g
      TYPE(sccs_control_type), POINTER                   :: sccs_control
      TYPE(pw_poisson_type), POINTER                     :: poisson_env
      TYPE(pw_pool_type), POINTER                        :: pw_pool
      REAL(KIND=dp), INTENT(IN)                          :: ion_bulk_density, ion_packing, kT, &
                                                            solute_charge, cell_volume, tolerance
      INTEGER, INTENT(IN)                                :: max_outer, output_unit
      LOGICAL, INTENT(IN)                                :: adaptive_tolerance
      REAL(KIND=dp), INTENT(INOUT)                       :: phi_shift
      INTEGER, INTENT(OUT)                               :: nouter, npcgtot
      REAL(KIND=dp), INTENT(OUT)                         :: residual_rms
      LOGICAL, INTENT(OUT)                               :: converged

      INTEGER                                            :: i, ipcgit, j, k, line_iter
      INTEGER, DIMENSION(3)                              :: lb, ub
      LOGICAL                                            :: pcg_converged
      REAL(KIND=dp)                                      :: alpha, charge_error, free_energy, &
                                                            free_energy_old, slope, &
                                                            target_tolerance, tol_linear
      TYPE(pw_r3d_rs_type)                               :: delta_response, dphi, kappa, phi_base, &
                                                            phi_solute, residual, response_base, &
                                                            response_nonlinear, work
      TYPE(pw_r3d_rs_type), DIMENSION(6)                 :: chi

      lb(:) = rho_response%pw_grid%bounds_local(1, :)
      ub(:) = rho_response%pw_grid%bounds_local(2, :)
      CALL pw_pool%create_pw(delta_response)
      CALL pw_pool%create_pw(dphi)
      CALL pw_pool%create_pw(kappa)
      CALL pw_pool%create_pw(phi_base)
      CALL pw_pool%create_pw(phi_solute)
      CALL pw_pool%create_pw(residual)
      CALL pw_pool%create_pw(response_base)
      CALL pw_pool%create_pw(response_nonlinear)
      CALL pw_pool%create_pw(work)
      DO i = 1, 6
         CALL pw_pool%create_pw(chi(i))
      END DO

      ! The uniform response component enforces charge neutrality.  It has no
      ! periodic electrostatic potential, but it fixes the ionic G=0 equation.
      charge_error = (solute_charge - pw_integrate_function(rho_response))/cell_volume
!$OMP PARALLEL DO DEFAULT(NONE) PRIVATE(i,j,k) SHARED(charge_error,lb,rho_response,ub)
      DO k = lb(3), ub(3)
         DO j = lb(2), ub(2)
            DO i = lb(1), ub(1)
               rho_response%array(i, j, k) = rho_response%array(i, j, k) + charge_error
            END DO
         END DO
      END DO
!$OMP END PARALLEL DO

      CALL pw_copy(rho_solute, rho_tot)
      CALL pw_axpy(rho_response, rho_tot)
      CALL sccs_nle_apply_preconditioner(rho_tot, phi, poisson_env, pw_pool)
      CALL sccs_nle_apply_preconditioner(rho_solute, phi_solute, poisson_env, pw_pool)
      CALL sccs_nle_evaluate(phi, rho_bound, rho_ion, lambda_diel, lambda_ion, eloc, &
                             sdiel, sion, kernel_g, sccs_control, pw_pool, ion_bulk_density, &
                             ion_packing, kT, solute_charge, phi_shift, lb, ub)

      nouter = 0
      npcgtot = 0
      converged = .FALSE.
      residual_rms = HUGE(1.0_dp)
      target_tolerance = tolerance
      DO
         CALL pw_copy(rho_bound, response_nonlinear)
         CALL pw_axpy(rho_ion, response_nonlinear)
         CALL pw_copy(response_nonlinear, residual)
         CALL pw_axpy(rho_response, residual, alpha=-1.0_dp)
         charge_error = pw_integrate_function(residual)/cell_volume
!$OMP    PARALLEL DO DEFAULT(NONE) PRIVATE(i,j,k) SHARED(charge_error,lb,residual,ub)
         DO k = lb(3), ub(3)
            DO j = lb(2), ub(2)
               DO i = lb(1), ub(1)
                  residual%array(i, j, k) = residual%array(i, j, k) - charge_error
               END DO
            END DO
         END DO
!$OMP    END PARALLEL DO
         CALL sccs_nle_apply_preconditioner(residual, work, poisson_env, pw_pool)
         residual_rms = SQRT(MAX(pw_integral_ab(work, work)/cell_volume, 0.0_dp))
         IF (adaptive_tolerance .AND. (nouter == 0)) THEN
            target_tolerance = MAX(tolerance, 0.1_dp*residual_rms)
         END IF
         IF (sccs_control%nle_force_strict .AND. (output_unit > 0)) THEN
            WRITE (UNIT=output_unit, FMT="(T3,A,I0,A,ES12.4,A,ES12.4)") &
               "SCCS| Electrolyte Newton update ", nouter, " residual ", residual_rms, &
               " target ", target_tolerance
         END IF
         IF (residual_rms <= target_tolerance) THEN
            converged = .TRUE.
            EXIT
         END IF
         IF (nouter >= max_outer) EXIT

         CALL sccs_nle_build_linear_response(chi, kappa, eloc, phi, sdiel, sion, sccs_control, &
                                             ion_bulk_density, ion_packing, kT, phi_shift, lb, ub)
         tol_linear = MAX(residual_rms/10.0_dp, 10.0_dp*EPSILON(1.0_dp))
         CALL sccs_nle_pcg(residual, dphi, chi, kappa, kernel_g, poisson_env, pw_pool, &
                           cell_volume, tol_linear, ipcgit, pcg_converged)
         npcgtot = npcgtot + ipcgit
         IF (.NOT. pcg_converged) THEN
            IF (sccs_control%nle_force_strict .AND. (output_unit > 0)) THEN
               WRITE (UNIT=output_unit, FMT="(T3,A,I0,A)") &
                  "SCCS| Electrolyte Newton linear solve stopped after ", ipcgit, " PCG iterations"
            END IF
            EXIT
         END IF

         free_energy_old = sccs_nle_free_energy(rho_response, phi, phi_solute, lambda_diel, &
                                                lambda_ion, sdiel, sion, solute_charge, &
                                                phi_shift, work)
         slope = pw_integral_ab(residual, dphi)
         CALL sccs_nle_apply_vacuum_operator(dphi, delta_response, pw_pool)
         CALL pw_copy(phi, phi_base)
         CALL pw_copy(rho_response, response_base)
         alpha = 1.0_dp
         DO line_iter = 0, 10
            CALL pw_copy(phi_base, phi)
            CALL pw_axpy(dphi, phi, alpha=alpha)
            CALL pw_copy(response_base, rho_response)
            CALL pw_axpy(delta_response, rho_response, alpha=alpha)
            CALL pw_copy(rho_solute, rho_tot)
            CALL pw_axpy(rho_response, rho_tot)
            CALL sccs_nle_evaluate(phi, rho_bound, rho_ion, lambda_diel, lambda_ion, eloc, &
                                   sdiel, sion, kernel_g, sccs_control, pw_pool, ion_bulk_density, &
                                   ion_packing, kT, solute_charge, phi_shift, lb, ub)
            free_energy = sccs_nle_free_energy(rho_response, phi, phi_solute, lambda_diel, &
                                               lambda_ion, sdiel, sion, solute_charge, &
                                               phi_shift, work)
            IF (free_energy >= free_energy_old + 0.1_dp*alpha*slope) EXIT
            IF (line_iter < 10) alpha = 0.5_dp*alpha
         END DO
         nouter = nouter + 1
         IF (sccs_control%nle_force_strict .AND. (output_unit > 0)) THEN
            WRITE (UNIT=output_unit, FMT="(T3,A,I0,A,I0,A,ES12.4,A,ES16.8)") &
               "SCCS| Electrolyte Newton accepted update ", nouter, " with ", ipcgit, &
               " PCG iterations; alpha ", alpha, "; free energy ", free_energy
         END IF
      END DO

      CALL pw_copy(rho_solute, rho_tot)
      CALL pw_axpy(rho_response, rho_tot)
      DO i = 1, 6
         CALL pw_pool%give_back_pw(chi(i))
      END DO
      CALL pw_pool%give_back_pw(work)
      CALL pw_pool%give_back_pw(response_nonlinear)
      CALL pw_pool%give_back_pw(response_base)
      CALL pw_pool%give_back_pw(residual)
      CALL pw_pool%give_back_pw(phi_solute)
      CALL pw_pool%give_back_pw(phi_base)
      CALL pw_pool%give_back_pw(kappa)
      CALL pw_pool%give_back_pw(dphi)
      CALL pw_pool%give_back_pw(delta_response)

   END SUBROUTINE sccs_nle_newton_solve

! **************************************************************************************************
!> \brief Evaluate nonlinear dielectric and ionic response fields at a fixed total potential.
! **************************************************************************************************
   SUBROUTINE sccs_nle_evaluate(phi, rho_bound, rho_ion, lambda_diel, lambda_ion, eloc, &
                                sdiel, sion, kernel_g, sccs_control, pw_pool, ion_bulk_density, &
                                ion_packing, kT, solute_charge, phi_shift, lb, ub)

      TYPE(pw_r3d_rs_type), INTENT(IN)                   :: phi
      TYPE(pw_r3d_rs_type), INTENT(INOUT)                :: rho_bound, rho_ion, lambda_diel, &
                                                            lambda_ion
      TYPE(pw_r3d_rs_type), DIMENSION(3), INTENT(INOUT)  :: eloc
      TYPE(pw_r3d_rs_type), INTENT(IN)                   :: sdiel, sion
      TYPE(pw_c1d_gs_type), INTENT(IN)                   :: kernel_g
      TYPE(sccs_control_type), POINTER                   :: sccs_control
      TYPE(pw_pool_type), POINTER                        :: pw_pool
      REAL(KIND=dp), INTENT(IN)                          :: ion_bulk_density, ion_packing, kT, &
                                                            solute_charge
      REAL(KIND=dp), INTENT(INOUT)                       :: phi_shift
      INTEGER, DIMENSION(3), INTENT(IN)                  :: lb, ub

      CALL sccs_nle_update_dielectric_response(rho_bound, lambda_diel, eloc, phi, sdiel, &
                                               kernel_g, kT, sccs_control, pw_pool, lb, ub)
      IF (ion_bulk_density > 0.0_dp) THEN
         CALL sccs_nle_update_ion_density(rho_ion, phi, sion, sccs_control, ion_bulk_density, &
                                          ion_packing, kT, solute_charge, phi_shift, pw_pool, lb, ub)
         CALL sccs_nle_update_ion_free_energy(lambda_ion, phi, sccs_control, ion_bulk_density, &
                                              ion_packing, kT, phi_shift, lb, ub)
      ELSE
         phi_shift = 0.0_dp
         CALL pw_zero(rho_ion)
         CALL pw_zero(lambda_ion)
      END IF

   END SUBROUTINE sccs_nle_evaluate

! **************************************************************************************************
!> \brief Build the analytic local dielectric susceptibility and ionic curvature.
! **************************************************************************************************
   SUBROUTINE sccs_nle_build_linear_response(chi, kappa, eloc, phi, sdiel, sion, sccs_control, &
                                             ion_bulk_density, ion_packing, kT, phi_shift, lb, ub)

      TYPE(pw_r3d_rs_type), DIMENSION(6), INTENT(INOUT)  :: chi
      TYPE(pw_r3d_rs_type), INTENT(INOUT)                :: kappa
      TYPE(pw_r3d_rs_type), DIMENSION(3), INTENT(IN)     :: eloc
      TYPE(pw_r3d_rs_type), INTENT(IN)                   :: phi, sdiel, sion
      TYPE(sccs_control_type), POINTER                   :: sccs_control
      REAL(KIND=dp), INTENT(IN)                          :: ion_bulk_density, ion_packing, kT, &
                                                            phi_shift
      INTEGER, DIMENSION(3), INTENT(IN)                  :: lb, ub

      INTEGER                                            :: i, j, k
      REAL(KIND=dp) :: alpha0_rot, alpha_par, alpha_perp, alpha_pol, beta_phi, chi_delta, chi_par, &
         chi_perp, coshm1, denom, e2, e_inv2, g_par, g_perp, invalpha_sic, x

      alpha0_rot = fourpi*sccs_control%nle_p_mol**2/(3.0_dp*kT)
      alpha_pol = alpha0_rot/(sccs_control%nle_eb_k - sccs_control%nle_epsilon_inf)* &
                  (sccs_control%nle_epsilon_inf - 1.0_dp)
      invalpha_sic = ((sccs_control%nle_eb_k - sccs_control%nle_epsilon_inf)/alpha0_rot - &
                      sccs_control%nle_n_mol)/(sccs_control%nle_eb_k - 1.0_dp)
!$OMP PARALLEL DO DEFAULT(NONE) &
!$OMP             PRIVATE(alpha_par,alpha_perp,beta_phi,chi_delta,chi_par,chi_perp,coshm1) &
!$OMP             PRIVATE(denom,e2,e_inv2,g_par,g_perp,i,j,k,x) &
!$OMP             SHARED(alpha0_rot,alpha_pol,chi,eloc,invalpha_sic,ion_bulk_density) &
!$OMP             SHARED(ion_packing,kT,kappa,lb,phi,phi_shift,sccs_control,sdiel,sion,ub)
      DO k = lb(3), ub(3)
         DO j = lb(2), ub(2)
            DO i = lb(1), ub(1)
               e2 = eloc(1)%array(i, j, k)**2 + eloc(2)%array(i, j, k)**2 + &
                    eloc(3)%array(i, j, k)**2
               x = sccs_control%nle_p_mol*SQRT(e2)/kT
               IF (x < 2.0E-4_dp) THEN
                  g_par = 1.0_dp
                  g_perp = 1.0_dp
                  e_inv2 = 0.0_dp
               ELSE IF (x > 100.0_dp) THEN
                  g_par = 3.0_dp/(x*x)
                  g_perp = 3.0_dp*(1.0_dp - 1.0_dp/x)/x
                  e_inv2 = 1.0_dp/e2
               ELSE
                  g_par = 3.0_dp*(1.0_dp/(x*x) - 1.0_dp/(SINH(x)**2))
                  g_perp = sccs_nle_grot(x)
                  e_inv2 = 1.0_dp/e2
               END IF
               alpha_par = alpha_pol + alpha0_rot*g_par
               alpha_perp = alpha_pol + alpha0_rot*g_perp
               chi_par = sccs_control%nle_n_mol*sdiel%array(i, j, k)/ &
                         (1.0_dp/alpha_par - invalpha_sic)
               chi_perp = sccs_control%nle_n_mol*sdiel%array(i, j, k)/ &
                          (1.0_dp/alpha_perp - invalpha_sic)
               chi_delta = chi_par - chi_perp
               chi(1)%array(i, j, k) = chi_perp + chi_delta*eloc(1)%array(i, j, k)**2*e_inv2
               chi(2)%array(i, j, k) = chi_perp + chi_delta*eloc(2)%array(i, j, k)**2*e_inv2
               chi(3)%array(i, j, k) = chi_perp + chi_delta*eloc(3)%array(i, j, k)**2*e_inv2
               chi(4)%array(i, j, k) = chi_delta*eloc(1)%array(i, j, k)*eloc(2)%array(i, j, k)*e_inv2
               chi(5)%array(i, j, k) = chi_delta*eloc(1)%array(i, j, k)*eloc(3)%array(i, j, k)*e_inv2
               chi(6)%array(i, j, k) = chi_delta*eloc(2)%array(i, j, k)*eloc(3)%array(i, j, k)*e_inv2

               kappa%array(i, j, k) = 0.0_dp
               beta_phi = sccs_control%nle_zion*(phi%array(i, j, k) + phi_shift)/kT
               IF ((ion_bulk_density > 0.0_dp) .AND. (ABS(beta_phi) < 40.0_dp)) THEN
                  coshm1 = COSH(beta_phi) - 1.0_dp
                  denom = 1.0_dp + ion_packing*coshm1
                  kappa%array(i, j, k) = 2.0_dp*sccs_control%nle_zion**2*ion_bulk_density* &
                                         sion%array(i, j, k)/kT*(1.0_dp + (1.0_dp - ion_packing)*coshm1)/(denom*denom)
               END IF
            END DO
         END DO
      END DO
!$OMP END PARALLEL DO

   END SUBROUTINE sccs_nle_build_linear_response

! **************************************************************************************************
!> \brief Apply the linearized nonlinear-electrolyte operator.
! **************************************************************************************************
   SUBROUTINE sccs_nle_apply_linear_operator(field_g, result_g, chi, kappa, kernel_g, pw_pool)

      TYPE(pw_c1d_gs_type), INTENT(IN)                   :: field_g
      TYPE(pw_c1d_gs_type), INTENT(INOUT)                :: result_g
      TYPE(pw_r3d_rs_type), DIMENSION(6), INTENT(IN)     :: chi
      TYPE(pw_r3d_rs_type), INTENT(IN)                   :: kappa
      TYPE(pw_c1d_gs_type), INTENT(IN)                   :: kernel_g
      TYPE(pw_pool_type), POINTER                        :: pw_pool

      INTEGER                                            :: i, idir, ig, j, k
      INTEGER, DIMENSION(3)                              :: derivative, lb, ub
      REAL(KIND=dp)                                      :: kappa_integral, weighted_average
      TYPE(pw_c1d_gs_type)                               :: work_g
      TYPE(pw_r3d_rs_type)                               :: field, work
      TYPE(pw_r3d_rs_type), DIMENSION(3)                 :: flux, gradient

      CALL pw_pool%create_pw(field)
      CALL pw_pool%create_pw(work)
      CALL pw_pool%create_pw(work_g)
      CALL pw_transfer(field_g, field)
      lb(:) = field%pw_grid%bounds_local(1, :)
      ub(:) = field%pw_grid%bounds_local(2, :)

      ! Vacuum part, kept in reciprocal space.
      CALL pw_copy(field_g, result_g)
      CALL pw_laplace(result_g)
      CALL pw_scale(result_g, -1.0_dp/fourpi)

      DO idir = 1, 3
         CALL pw_pool%create_pw(flux(idir))
         CALL pw_pool%create_pw(gradient(idir))
         CALL pw_copy(field_g, work_g)
         derivative(:) = 0
         derivative(idir) = 1
         CALL pw_derive(work_g, derivative)
!$OMP    PARALLEL DO DEFAULT(NONE) PRIVATE(ig) SHARED(kernel_g,work_g)
         DO ig = 1, SIZE(work_g%array)
            work_g%array(ig) = work_g%array(ig)*kernel_g%array(ig)
         END DO
!$OMP    END PARALLEL DO
         CALL pw_transfer(work_g, gradient(idir))
      END DO
!$OMP PARALLEL DO DEFAULT(NONE) PRIVATE(i,j,k) SHARED(chi,flux,gradient,lb,ub)
      DO k = lb(3), ub(3)
         DO j = lb(2), ub(2)
            DO i = lb(1), ub(1)
               flux(1)%array(i, j, k) = (chi(1)%array(i, j, k)*gradient(1)%array(i, j, k) + &
                                         chi(4)%array(i, j, k)*gradient(2)%array(i, j, k) + &
                                         chi(5)%array(i, j, k)*gradient(3)%array(i, j, k))/fourpi
               flux(2)%array(i, j, k) = (chi(4)%array(i, j, k)*gradient(1)%array(i, j, k) + &
                                         chi(2)%array(i, j, k)*gradient(2)%array(i, j, k) + &
                                         chi(6)%array(i, j, k)*gradient(3)%array(i, j, k))/fourpi
               flux(3)%array(i, j, k) = (chi(5)%array(i, j, k)*gradient(1)%array(i, j, k) + &
                                         chi(6)%array(i, j, k)*gradient(2)%array(i, j, k) + &
                                         chi(3)%array(i, j, k)*gradient(3)%array(i, j, k))/fourpi
            END DO
         END DO
      END DO
!$OMP END PARALLEL DO
      DO idir = 1, 3
         CALL pw_transfer(flux(idir), work_g)
!$OMP    PARALLEL DO DEFAULT(NONE) PRIVATE(ig) SHARED(kernel_g,work_g)
         DO ig = 1, SIZE(work_g%array)
            work_g%array(ig) = work_g%array(ig)*kernel_g%array(ig)
         END DO
!$OMP    END PARALLEL DO
         derivative(:) = 0
         derivative(idir) = 1
         CALL pw_derive(work_g, derivative)
         CALL pw_axpy(work_g, result_g, alpha=-1.0_dp)
      END DO

      kappa_integral = pw_integrate_function(kappa)
      weighted_average = 0.0_dp
      IF (kappa_integral > EPSILON(1.0_dp)) THEN
         weighted_average = pw_integral_ab(kappa, field)/kappa_integral
      END IF
!$OMP PARALLEL DO DEFAULT(NONE) PRIVATE(i,j,k) &
!$OMP             SHARED(field,kappa,lb,ub,weighted_average,work)
      DO k = lb(3), ub(3)
         DO j = lb(2), ub(2)
            DO i = lb(1), ub(1)
               work%array(i, j, k) = kappa%array(i, j, k)* &
                                     (field%array(i, j, k) - weighted_average)
            END DO
         END DO
      END DO
!$OMP END PARALLEL DO
      CALL pw_transfer(work, work_g)
      CALL pw_axpy(work_g, result_g)
      DO idir = 1, 3
         CALL pw_pool%give_back_pw(gradient(idir))
         CALL pw_pool%give_back_pw(flux(idir))
      END DO
      CALL pw_pool%give_back_pw(work_g)
      CALL pw_pool%give_back_pw(work)
      CALL pw_pool%give_back_pw(field)

   END SUBROUTINE sccs_nle_apply_linear_operator

! **************************************************************************************************
!> \brief Apply the periodic vacuum Poisson operator -laplacian/(4*pi).
! **************************************************************************************************
   SUBROUTINE sccs_nle_apply_vacuum_operator(phi, density, pw_pool)

      TYPE(pw_r3d_rs_type), INTENT(IN)                   :: phi
      TYPE(pw_r3d_rs_type), INTENT(INOUT)                :: density
      TYPE(pw_pool_type), POINTER                        :: pw_pool

      TYPE(pw_c1d_gs_type)                               :: work_g

      CALL pw_pool%create_pw(work_g)
      CALL pw_transfer(phi, work_g)
      CALL pw_laplace(work_g)
      CALL pw_transfer(work_g, density)
      CALL pw_scale(density, -1.0_dp/fourpi)
      CALL pw_pool%give_back_pw(work_g)

   END SUBROUTINE sccs_nle_apply_vacuum_operator

! **************************************************************************************************
!> \brief Apply the periodic vacuum inverse Poisson operator.
! **************************************************************************************************
   SUBROUTINE sccs_nle_apply_preconditioner(density, potential, poisson_env, pw_pool)

      TYPE(pw_r3d_rs_type), INTENT(IN)                   :: density
      TYPE(pw_r3d_rs_type), INTENT(INOUT)                :: potential
      TYPE(pw_poisson_type), POINTER                     :: poisson_env
      TYPE(pw_pool_type), POINTER                        :: pw_pool

      TYPE(pw_c1d_gs_type)                               :: density_g, potential_g

      CALL pw_pool%create_pw(density_g)
      CALL pw_pool%create_pw(potential_g)
      CALL pw_transfer(density, density_g)
      CALL pw_poisson_solve(poisson_env=poisson_env, density=density_g, vhartree=potential_g)
      CALL pw_transfer(potential_g, potential)
      CALL pw_pool%give_back_pw(potential_g)
      CALL pw_pool%give_back_pw(density_g)

   END SUBROUTINE sccs_nle_apply_preconditioner

! **************************************************************************************************
!> \brief Solve one linearized Newton equation by preconditioned conjugate gradients.
! **************************************************************************************************
   SUBROUTINE sccs_nle_pcg(rhs, solution, chi, kappa, kernel_g, poisson_env, pw_pool, &
                           cell_volume, tolerance, niter, converged)

      TYPE(pw_r3d_rs_type), INTENT(IN)                   :: rhs
      TYPE(pw_r3d_rs_type), INTENT(INOUT)                :: solution
      TYPE(pw_r3d_rs_type), DIMENSION(6), INTENT(IN)     :: chi
      TYPE(pw_r3d_rs_type), INTENT(IN)                   :: kappa
      TYPE(pw_c1d_gs_type), INTENT(IN)                   :: kernel_g
      TYPE(pw_poisson_type), POINTER                     :: poisson_env
      TYPE(pw_pool_type), POINTER                        :: pw_pool
      REAL(KIND=dp), INTENT(IN)                          :: cell_volume, tolerance
      INTEGER, INTENT(OUT)                               :: niter
      LOGICAL, INTENT(OUT)                               :: converged

      REAL(KIND=dp)                                      :: alpha, beta, pap, rms, rz, rz_new
      TYPE(pw_c1d_gs_type)                               :: ap, p, residual, solution_g, z

      CALL pw_pool%create_pw(ap)
      CALL pw_pool%create_pw(p)
      CALL pw_pool%create_pw(residual)
      CALL pw_pool%create_pw(solution_g)
      CALL pw_pool%create_pw(z)
      CALL pw_zero(solution_g)
      CALL pw_transfer(rhs, residual)
      CALL pw_poisson_solve(poisson_env=poisson_env, density=residual, vhartree=z)
      CALL pw_copy(z, p)
      rz = pw_integral_ab(residual, z)
      converged = .FALSE.
      niter = 0
      DO niter = 1, 200
         CALL sccs_nle_apply_linear_operator(p, ap, chi, kappa, kernel_g, pw_pool)
         pap = pw_integral_ab(p, ap)
         IF (pap <= EPSILON(1.0_dp)) EXIT
         alpha = rz/pap
         CALL pw_axpy(p, solution_g, alpha=alpha)
         CALL pw_axpy(ap, residual, alpha=-alpha)
         CALL pw_poisson_solve(poisson_env=poisson_env, density=residual, vhartree=z)
         rms = SQRT(MAX(pw_integral_ab(z, z)/cell_volume, 0.0_dp))
         IF ((rms <= tolerance) .AND. (niter >= 4)) THEN
            converged = .TRUE.
            EXIT
         END IF
         rz_new = pw_integral_ab(residual, z)
         beta = rz_new/rz
         CALL pw_axpy(z, p, alpha=1.0_dp, beta=beta)
         rz = rz_new
      END DO
      CALL pw_transfer(solution_g, solution)
      CALL pw_pool%give_back_pw(z)
      CALL pw_pool%give_back_pw(solution_g)
      CALL pw_pool%give_back_pw(residual)
      CALL pw_pool%give_back_pw(p)
      CALL pw_pool%give_back_pw(ap)

   END SUBROUTINE sccs_nle_pcg

! **************************************************************************************************
!> \brief Electrolyte free-energy functional used by the Newton line search.
! **************************************************************************************************
   REAL(KIND=dp) FUNCTION sccs_nle_free_energy(rho_response, phi, phi_solute, lambda_diel, &
                                               lambda_ion, sdiel, sion, solute_charge, &
                                               phi_shift, work) RESULT(free_energy)

      TYPE(pw_r3d_rs_type), INTENT(IN)                   :: rho_response, phi, phi_solute, &
                                                            lambda_diel, lambda_ion, sdiel, sion
      REAL(KIND=dp), INTENT(IN)                          :: solute_charge, phi_shift
      TYPE(pw_r3d_rs_type), INTENT(INOUT)                :: work

      CALL pw_copy(phi, work)
      CALL pw_axpy(phi_solute, work, alpha=-1.0_dp)
      free_energy = -solute_charge*phi_shift - 0.5_dp*pw_integral_ab(rho_response, work) + &
                    pw_integral_ab(sdiel, lambda_diel) + pw_integral_ab(sion, lambda_ion)

   END FUNCTION sccs_nle_free_energy

! **************************************************************************************************
!> \brief Update the mobile-ion grand-potential density.
! **************************************************************************************************
   SUBROUTINE sccs_nle_update_ion_free_energy(lambda_ion, phi, sccs_control, ion_bulk_density, &
                                              ion_packing, kT, phi_shift, lb, ub)

      TYPE(pw_r3d_rs_type), INTENT(INOUT)                :: lambda_ion
      TYPE(pw_r3d_rs_type), INTENT(IN)                   :: phi
      TYPE(sccs_control_type), POINTER                   :: sccs_control
      REAL(KIND=dp), INTENT(IN)                          :: ion_bulk_density, ion_packing, kT, &
                                                            phi_shift
      INTEGER, DIMENSION(3), INTENT(IN)                  :: lb, ub

      INTEGER                                            :: i, j, k
      REAL(KIND=dp)                                      :: beta_phi, ion_denom

!$OMP PARALLEL DO DEFAULT(NONE) PRIVATE(beta_phi,i,ion_denom,j,k) &
!$OMP             SHARED(ion_bulk_density,ion_packing,kT,lambda_ion,lb,phi,phi_shift,sccs_control,ub)
      DO k = lb(3), ub(3)
         DO j = lb(2), ub(2)
            DO i = lb(1), ub(1)
               beta_phi = sccs_control%nle_zion*(phi%array(i, j, k) + phi_shift)/kT
               beta_phi = MAX(-40.0_dp, MIN(40.0_dp, beta_phi))
               IF (ion_packing > 0.0_dp) THEN
                  ion_denom = 1.0_dp + ion_packing*(COSH(beta_phi) - 1.0_dp)
                  lambda_ion%array(i, j, k) = -2.0_dp*kT*ion_bulk_density*LOG(ion_denom)/ion_packing
               ELSE
                  lambda_ion%array(i, j, k) = -2.0_dp*kT*ion_bulk_density*(COSH(beta_phi) - 1.0_dp)
               END IF
            END DO
         END DO
      END DO
!$OMP END PARALLEL DO

   END SUBROUTINE sccs_nle_update_ion_free_energy

! **************************************************************************************************
!> \brief Update nonlinear electrolyte mobile ion density with a neutralising G=0 potential shift.
!> \param rho_ion ...
!> \param phi ...
!> \param sion ...
!> \param sccs_control ...
!> \param ion_bulk_density ...
!> \param ion_packing ...
!> \param kT ...
!> \param target_charge ...
!> \param phi_shift ...
!> \param lb ...
!> \param ub ...
! **************************************************************************************************
   SUBROUTINE sccs_nle_update_ion_density(rho_ion, phi, sion, sccs_control, ion_bulk_density, &
                                          ion_packing, kT, target_charge, phi_shift, pw_pool, lb, ub)

      TYPE(pw_r3d_rs_type), INTENT(INOUT)                :: rho_ion
      TYPE(pw_r3d_rs_type), INTENT(IN)                   :: phi, sion
      TYPE(sccs_control_type), POINTER                   :: sccs_control
      REAL(KIND=dp), INTENT(IN)                          :: ion_bulk_density, ion_packing, kT, &
                                                            target_charge
      REAL(KIND=dp), INTENT(INOUT)                       :: phi_shift
      TYPE(pw_pool_type), POINTER                        :: pw_pool
      INTEGER, DIMENSION(3), INTENT(IN)                  :: lb, ub

      INTEGER                                            :: bracket_iter, iter
      LOGICAL                                            :: shift_converged
      REAL(KIND=dp)                                      :: charge_derivative, charge_high, &
                                                            charge_low, charge_mid, high, low, &
                                                            mid, shift_step
      TYPE(pw_r3d_rs_type)                               :: rho_derivative

      CALL pw_pool%create_pw(rho_derivative)
      shift_converged = .FALSE.
      DO iter = 1, 12
         CALL sccs_nle_set_ion_density(rho_ion, rho_derivative, phi, sion, sccs_control, &
                                       ion_bulk_density, ion_packing, kT, phi_shift, charge_mid, &
                                       charge_derivative, lb, ub)
         IF (ABS(charge_mid - target_charge) < 1.0E-10_dp) THEN
            shift_converged = .TRUE.
            EXIT
         END IF
         IF (ABS(charge_derivative) <= EPSILON(1.0_dp)) EXIT
         shift_step = -(charge_mid - target_charge)/charge_derivative
         shift_step = MAX(-0.5_dp, MIN(0.5_dp, shift_step))
         phi_shift = phi_shift + shift_step
      END DO
      IF (shift_converged) THEN
         CALL pw_pool%give_back_pw(rho_derivative)
         RETURN
      END IF

      ! Safeguarded fallback for potentials close to ionic saturation.
      low = phi_shift - 1.0_dp
      high = phi_shift + 1.0_dp
      CALL sccs_nle_set_ion_density(rho_ion, rho_derivative, phi, sion, sccs_control, &
                                    ion_bulk_density, ion_packing, kT, low, charge_low, &
                                    charge_derivative, lb, ub)
      CALL sccs_nle_set_ion_density(rho_ion, rho_derivative, phi, sion, sccs_control, &
                                    ion_bulk_density, ion_packing, kT, high, charge_high, &
                                    charge_derivative, lb, ub)
      DO bracket_iter = 1, 12
         IF ((target_charge >= MIN(charge_low, charge_high) - 1.0E-8_dp) .AND. &
             (target_charge <= MAX(charge_low, charge_high) + 1.0E-8_dp)) EXIT
         low = low - 2.0_dp**bracket_iter
         high = high + 2.0_dp**bracket_iter
         CALL sccs_nle_set_ion_density(rho_ion, rho_derivative, phi, sion, sccs_control, &
                                       ion_bulk_density, ion_packing, kT, low, charge_low, &
                                       charge_derivative, lb, ub)
         CALL sccs_nle_set_ion_density(rho_ion, rho_derivative, phi, sion, sccs_control, &
                                       ion_bulk_density, ion_packing, kT, high, charge_high, &
                                       charge_derivative, lb, ub)
      END DO
      IF ((target_charge < MIN(charge_low, charge_high) - 1.0E-8_dp) .OR. &
          (target_charge > MAX(charge_low, charge_high) + 1.0E-8_dp)) THEN
         CPABORT("SCCS nonlinear electrolyte cannot neutralise the requested charge.")
      END IF

      DO iter = 1, 60
         mid = 0.5_dp*(low + high)
         CALL sccs_nle_set_ion_density(rho_ion, rho_derivative, phi, sion, sccs_control, &
                                       ion_bulk_density, ion_packing, kT, mid, charge_mid, &
                                       charge_derivative, lb, ub)
         IF (ABS(charge_mid - target_charge) < 1.0E-10_dp) EXIT
         IF (charge_high > charge_low) THEN
            IF (charge_mid < target_charge) THEN
               low = mid
            ELSE
               high = mid
            END IF
         ELSE
            IF (charge_mid > target_charge) THEN
               low = mid
            ELSE
               high = mid
            END IF
         END IF
      END DO

      phi_shift = 0.5_dp*(low + high)
      CALL sccs_nle_set_ion_density(rho_ion, rho_derivative, phi, sion, sccs_control, &
                                    ion_bulk_density, ion_packing, kT, phi_shift, charge_mid, &
                                    charge_derivative, lb, ub)
      CALL pw_pool%give_back_pw(rho_derivative)

   END SUBROUTINE sccs_nle_update_ion_density

! **************************************************************************************************
!> \brief Fill nonlinear electrolyte mobile ion density for a fixed potential shift.
!> \param rho_ion ...
!> \param phi ...
!> \param sion ...
!> \param sccs_control ...
!> \param ion_bulk_density ...
!> \param ion_packing ...
!> \param kT ...
!> \param phi_shift ...
!> \param charge ...
!> \param lb ...
!> \param ub ...
! **************************************************************************************************
   SUBROUTINE sccs_nle_set_ion_density(rho_ion, rho_derivative, phi, sion, sccs_control, &
                                       ion_bulk_density, ion_packing, kT, phi_shift, charge, &
                                       charge_derivative, lb, ub)

      TYPE(pw_r3d_rs_type), INTENT(INOUT)                :: rho_ion, rho_derivative
      TYPE(pw_r3d_rs_type), INTENT(IN)                   :: phi, sion
      TYPE(sccs_control_type), POINTER                   :: sccs_control
      REAL(KIND=dp), INTENT(IN)                          :: ion_bulk_density, ion_packing, kT, &
                                                            phi_shift
      REAL(KIND=dp), INTENT(OUT)                         :: charge, charge_derivative
      INTEGER, DIMENSION(3), INTENT(IN)                  :: lb, ub

      INTEGER                                            :: i, j, k
      REAL(KIND=dp)                                      :: beta_phi, coshm1, ion_denom

!$OMP PARALLEL DO DEFAULT(NONE) &
!$OMP             PRIVATE(beta_phi,coshm1,i,ion_denom,j,k) &
!$OMP             SHARED(ion_bulk_density,ion_packing,kT,lb,phi,phi_shift,rho_ion,sccs_control) &
!$OMP             SHARED(rho_derivative,sion,ub)
      DO k = lb(3), ub(3)
         DO j = lb(2), ub(2)
            DO i = lb(1), ub(1)
               beta_phi = MAX(-40.0_dp, MIN(40.0_dp, &
                                            sccs_control%nle_zion*(phi%array(i, j, k) + phi_shift)/kT))
               coshm1 = COSH(beta_phi) - 1.0_dp
               ion_denom = 1.0_dp + ion_packing*coshm1
               rho_ion%array(i, j, k) = -2.0_dp*sccs_control%nle_zion*ion_bulk_density* &
                                        sion%array(i, j, k)*SINH(beta_phi)/ion_denom
               rho_derivative%array(i, j, k) = -2.0_dp*sccs_control%nle_zion**2* &
                                               ion_bulk_density*sion%array(i, j, k)/kT* &
                                               (1.0_dp + (1.0_dp - ion_packing)*coshm1)/(ion_denom*ion_denom)
            END DO
         END DO
      END DO
!$OMP END PARALLEL DO

      charge = pw_integrate_function(rho_ion)
      charge_derivative = pw_integrate_function(rho_derivative)

   END SUBROUTINE sccs_nle_set_ion_density

! **************************************************************************************************
!> \brief Build the nonlinear, nonlocal nonlinear electrolyte dielectric response for a fixed potential.
! **************************************************************************************************
   SUBROUTINE sccs_nle_update_dielectric_response(rho_bound, lambda_diel, eloc, phi, sdiel, &
                                                  kernel_g, kT, sccs_control, pw_pool, lb, ub)

      TYPE(pw_r3d_rs_type), INTENT(INOUT)                :: rho_bound, lambda_diel
      TYPE(pw_r3d_rs_type), DIMENSION(3), INTENT(INOUT)  :: eloc
      TYPE(pw_r3d_rs_type), INTENT(IN)                   :: phi, sdiel
      TYPE(pw_c1d_gs_type), INTENT(IN)                   :: kernel_g
      REAL(KIND=dp), INTENT(IN)                          :: kT
      TYPE(sccs_control_type), POINTER                   :: sccs_control
      TYPE(pw_pool_type), POINTER                        :: pw_pool
      INTEGER, DIMENSION(3), INTENT(IN)                  :: lb, ub

      INTEGER                                            :: i, idir, iter, j, k
      REAL(KIND=dp)                                      :: alpha0_rot, alpha_pol, coeff, e0, emag, &
                                                            f_high, f_loc, f_low, f_mid, g_rot, &
                                                            high, invalpha_sic, lambda, low, &
                                                            nwork, p_beta, x, x0
      TYPE(pw_r3d_rs_type), DIMENSION(3)                 :: polarisation

      alpha0_rot = fourpi*sccs_control%nle_p_mol**2/(3.0_dp*kT)
      alpha_pol = alpha0_rot/(sccs_control%nle_eb_k - sccs_control%nle_epsilon_inf)* &
                  (sccs_control%nle_epsilon_inf - 1.0_dp)
      invalpha_sic = ((sccs_control%nle_eb_k - sccs_control%nle_epsilon_inf)/alpha0_rot - &
                      sccs_control%nle_n_mol)/(sccs_control%nle_eb_k - 1.0_dp)
      IF ((alpha0_rot <= 0.0_dp) .OR. (invalpha_sic < 0.0_dp)) THEN
         CPABORT("Invalid nonlinear electrolyte nonlinear dielectric parameters.")
      END IF
      f_low = 1.0_dp/(1.0_dp - alpha_pol*invalpha_sic)
      f_high = 1.0_dp/(1.0_dp - (alpha_pol + alpha0_rot)*invalpha_sic)
      IF ((f_low <= 0.0_dp) .OR. (f_high < f_low)) THEN
         CPABORT("Invalid nonlinear electrolyte local-field correction interval.")
      END IF
      p_beta = sccs_control%nle_p_mol/kT
      CALL sccs_nle_smoothed_gradient(phi, eloc, kernel_g, pw_pool)
      DO idir = 1, 3
         CALL pw_scale(eloc(idir), -1.0_dp)
         CALL pw_pool%create_pw(polarisation(idir))
      END DO

!$OMP PARALLEL DO DEFAULT(NONE) &
!$OMP             PRIVATE(coeff,e0,emag,f_loc,f_mid,g_rot,high,i,idir,iter,j,k,lambda,low,nwork,x,x0) &
!$OMP             SHARED(alpha0_rot,alpha_pol,eloc,f_high,f_low,invalpha_sic,kT,lambda_diel,lb) &
!$OMP             SHARED(p_beta,polarisation,sccs_control,sdiel,ub)
      DO k = lb(3), ub(3)
         DO j = lb(2), ub(2)
            DO i = lb(1), ub(1)
               e0 = SQRT(eloc(1)%array(i, j, k)**2 + eloc(2)%array(i, j, k)**2 + &
                         eloc(3)%array(i, j, k)**2)
               x0 = p_beta*e0
               low = f_low
               high = f_high
               IF (x0 > 0.0_dp) THEN
                  DO iter = 1, 40
                     f_mid = 0.5_dp*(low + high)
                     g_rot = sccs_nle_grot(f_mid*x0)
                     IF (1.0_dp/(1.0_dp - (g_rot*alpha0_rot + alpha_pol)*invalpha_sic) - &
                         f_mid > 0.0_dp) THEN
                        low = f_mid
                     ELSE
                        high = f_mid
                     END IF
                  END DO
               END IF
               f_loc = 0.5_dp*(low + high)
               DO idir = 1, 3
                  eloc(idir)%array(i, j, k) = f_loc*eloc(idir)%array(i, j, k)
               END DO
               emag = f_loc*e0
               x = p_beta*emag
               g_rot = sccs_nle_grot(x)
               nwork = (alpha0_rot*g_rot + alpha_pol)/fourpi
               coeff = sccs_control%nle_n_mol*sdiel%array(i, j, k)*nwork
               DO idir = 1, 3
                  polarisation(idir)%array(i, j, k) = coeff*eloc(idir)%array(i, j, k)
               END DO
               IF (x < 2.0E-4_dp) THEN
                  lambda = -kT*x*x/6.0_dp
               ELSE IF (x > 100.0_dp) THEN
                  lambda = -kT*(x - LOG(2.0_dp*x))
               ELSE
                  lambda = -kT*LOG(SINH(x)/x)
               END IF
               lambda = lambda - 0.5_dp*alpha_pol/fourpi*emag**2 + &
                        0.5_dp*fourpi*invalpha_sic*(nwork*emag)**2
               lambda_diel%array(i, j, k) = sccs_control%nle_n_mol*lambda
            END DO
         END DO
      END DO
!$OMP END PARALLEL DO

      CALL sccs_nle_gaussian_divergence(polarisation, rho_bound, kernel_g, pw_pool)
      DO idir = 1, 3
         CALL pw_pool%give_back_pw(polarisation(idir))
      END DO

   END SUBROUTINE sccs_nle_update_dielectric_response

! **************************************************************************************************
!> \brief nonlinear electrolyte rotational response factor.
! **************************************************************************************************
   PURE REAL(KIND=dp) FUNCTION sccs_nle_grot(x) RESULT(g_rot)

      REAL(KIND=dp), INTENT(IN)                          :: x

      IF (x < 2.0E-4_dp) THEN
         g_rot = 1.0_dp
      ELSE IF (x > 100.0_dp) THEN
         g_rot = 3.0_dp*(1.0_dp - 1.0_dp/x)/x
      ELSE
         g_rot = 3.0_dp*(1.0_dp/TANH(x) - 1.0_dp/x)/x
      END IF

   END FUNCTION sccs_nle_grot

! **************************************************************************************************
!> \brief Prepare a normalized Gaussian kernel in reciprocal space.
! **************************************************************************************************
   SUBROUTINE sccs_nle_prepare_gaussian_kernel(kernel_g, sigma, pw_pool)

      TYPE(pw_c1d_gs_type), INTENT(INOUT)                :: kernel_g
      REAL(KIND=dp), INTENT(IN)                          :: sigma
      TYPE(pw_pool_type), POINTER                        :: pw_pool

      REAL(KIND=dp)                                      :: factor
      TYPE(pw_r3d_rs_type)                               :: kernel_rs

      CALL pw_pool%create_pw(kernel_rs)
      CALL sccs_nle_fill_gaussian_kernel(kernel_rs, sigma)
      CALL pw_transfer(kernel_rs, kernel_g)
      factor = REAL(kernel_rs%pw_grid%ngpts, KIND=dp)*kernel_rs%pw_grid%dvol
      CALL pw_scale(kernel_g, factor)
      CALL pw_pool%give_back_pw(kernel_rs)

   END SUBROUTINE sccs_nle_prepare_gaussian_kernel

! **************************************************************************************************
!> \brief Compute the Gaussian-smoothed gradient of a scalar field.
! **************************************************************************************************
   SUBROUTINE sccs_nle_smoothed_gradient(field, gradient, kernel_g, pw_pool)

      TYPE(pw_r3d_rs_type), INTENT(IN)                   :: field
      TYPE(pw_r3d_rs_type), DIMENSION(3), INTENT(INOUT)  :: gradient
      TYPE(pw_c1d_gs_type), INTENT(IN)                   :: kernel_g
      TYPE(pw_pool_type), POINTER                        :: pw_pool

      INTEGER                                            :: idir, ig
      INTEGER, DIMENSION(3)                              :: derivative
      TYPE(pw_c1d_gs_type)                               :: field_g, work_g

      CALL pw_pool%create_pw(field_g)
      CALL pw_pool%create_pw(work_g)
      CALL pw_transfer(field, field_g)
      DO idir = 1, 3
         CALL pw_copy(field_g, work_g)
         derivative(:) = 0
         derivative(idir) = 1
         CALL pw_derive(work_g, derivative)
!$OMP    PARALLEL DO DEFAULT(NONE) PRIVATE(ig) SHARED(kernel_g,work_g)
         DO ig = 1, SIZE(work_g%array)
            work_g%array(ig) = work_g%array(ig)*kernel_g%array(ig)
         END DO
!$OMP    END PARALLEL DO
         CALL pw_transfer(work_g, gradient(idir))
      END DO
      CALL pw_pool%give_back_pw(work_g)
      CALL pw_pool%give_back_pw(field_g)

   END SUBROUTINE sccs_nle_smoothed_gradient

! **************************************************************************************************
!> \brief Compute minus the divergence of a Gaussian-convolved vector field.
! **************************************************************************************************
   SUBROUTINE sccs_nle_gaussian_divergence(field, divergence, kernel_g, pw_pool)

      TYPE(pw_r3d_rs_type), DIMENSION(3), INTENT(IN)     :: field
      TYPE(pw_r3d_rs_type), INTENT(INOUT)                :: divergence
      TYPE(pw_c1d_gs_type), INTENT(IN)                   :: kernel_g
      TYPE(pw_pool_type), POINTER                        :: pw_pool

      INTEGER                                            :: idir, ig
      INTEGER, DIMENSION(3)                              :: derivative
      TYPE(pw_c1d_gs_type)                               :: divergence_g, field_g

      CALL pw_pool%create_pw(field_g)
      CALL pw_pool%create_pw(divergence_g)
      CALL pw_zero(divergence_g)
      DO idir = 1, 3
         CALL pw_transfer(field(idir), field_g)
!$OMP    PARALLEL DO DEFAULT(NONE) PRIVATE(ig) SHARED(field_g,kernel_g)
         DO ig = 1, SIZE(field_g%array)
            field_g%array(ig) = field_g%array(ig)*kernel_g%array(ig)
         END DO
!$OMP    END PARALLEL DO
         derivative(:) = 0
         derivative(idir) = 1
         CALL pw_derive(field_g, derivative)
         CALL pw_axpy(field_g, divergence_g, alpha=-1.0_dp)
      END DO
      CALL pw_transfer(divergence_g, divergence)
      CALL pw_pool%give_back_pw(divergence_g)
      CALL pw_pool%give_back_pw(field_g)

   END SUBROUTINE sccs_nle_gaussian_divergence

! **************************************************************************************************
!> \brief Fill and discretely normalize the nonlinear electrolyte Gaussian kernel.
! **************************************************************************************************
   SUBROUTINE sccs_nle_fill_gaussian_kernel(kernel, sigma)

      TYPE(pw_r3d_rs_type), INTENT(INOUT)                :: kernel
      REAL(KIND=dp), INTENT(IN)                          :: sigma

      INTEGER                                            :: i, ii, j, jj, k, kk
      INTEGER, DIMENSION(2, 3)                           :: bounds
      INTEGER, DIMENSION(3)                              :: lb, npts, ub
      REAL(KIND=dp)                                      :: norm, r2
      REAL(KIND=dp), DIMENSION(3)                        :: rvec

      npts(:) = kernel%pw_grid%npts(:)
      bounds(:, :) = kernel%pw_grid%bounds(:, :)
      lb(:) = kernel%pw_grid%bounds_local(1, :)
      ub(:) = kernel%pw_grid%bounds_local(2, :)
!$OMP PARALLEL DO DEFAULT(NONE) PRIVATE(i,ii,j,jj,k,kk,r2,rvec) &
!$OMP             SHARED(bounds,kernel,lb,npts,sigma,ub)
      DO k = lb(3), ub(3)
         kk = MODULO(k - bounds(1, 3), npts(3))
         IF ((kk > 0) .AND. (2*kk >= npts(3))) kk = kk - npts(3)
         DO j = lb(2), ub(2)
            jj = MODULO(j - bounds(1, 2), npts(2))
            IF ((jj > 0) .AND. (2*jj >= npts(2))) jj = jj - npts(2)
            DO i = lb(1), ub(1)
               ii = MODULO(i - bounds(1, 1), npts(1))
               IF ((ii > 0) .AND. (2*ii >= npts(1))) ii = ii - npts(1)
               rvec(:) = REAL(ii, KIND=dp)*kernel%pw_grid%dh(:, 1) + &
                         REAL(jj, KIND=dp)*kernel%pw_grid%dh(:, 2) + &
                         REAL(kk, KIND=dp)*kernel%pw_grid%dh(:, 3)
               r2 = SUM(rvec(:)**2)
               kernel%array(i, j, k) = EXP(-0.5_dp*r2/sigma**2)/(sigma*SQRT(twopi))**3
            END DO
         END DO
      END DO
!$OMP END PARALLEL DO
      norm = pw_integrate_function(kernel)
      IF (norm <= 0.0_dp) CPABORT("Invalid nonlinear electrolyte Gaussian kernel normalization.")
      CALL pw_scale(kernel, 1.0_dp/norm)

   END SUBROUTINE sccs_nle_fill_gaussian_kernel

! **************************************************************************************************
!> \brief Build nonlinear electrolyte style cavity and electrolyte shape functions.
!> \param rho_elec ...
!> \param eps_elec ...
!> \param deps_elec ...
!> \param svdw ...
!> \param ssolv ...
!> \param sion ...
!> \param sdiel ...
!> \param sdiel_eff ...
!> \param scav ...
!> \param xvdw ...
!> \param xsolv ...
!> \param xion ...
!> \param xdiel ...
!> \param xcav ...
!> \param sccs_control ...
!> \param pw_pool ...
!> \param lb ...
!> \param ub ...
! **************************************************************************************************
   SUBROUTINE sccs_nle_build_shapes(rho_elec, eps_elec, deps_elec, svdw, ssolv, sion, sdiel, &
                                    sdiel_eff, scav, xvdw, xsolv, xion, xdiel, xcav, &
                                    kernel_solv, kernel_ion, kernel_diel, kernel_cav, &
                                    sccs_control, pw_pool, lb, ub)

      TYPE(pw_r3d_rs_type), INTENT(IN)                   :: rho_elec
      TYPE(pw_r3d_rs_type), INTENT(INOUT)                :: eps_elec, deps_elec, svdw, ssolv, sion, &
                                                            sdiel, sdiel_eff, scav, xvdw, xsolv, &
                                                            xion, xdiel, xcav
      TYPE(pw_c1d_gs_type), INTENT(IN)                   :: kernel_solv, kernel_ion, kernel_diel, &
                                                            kernel_cav
      TYPE(sccs_control_type), POINTER                   :: sccs_control
      TYPE(pw_pool_type), POINTER                        :: pw_pool
      INTEGER, DIMENSION(3), INTENT(IN)                  :: lb, ub

      INTEGER                                            :: i, j, k
      REAL(KIND=dp)                                      :: rho_scaled, x

!$OMP PARALLEL DO DEFAULT(NONE) &
!$OMP             PRIVATE(i,j,k,rho_scaled,x) &
!$OMP             SHARED(deps_elec,lb,rho_elec,sccs_control,ssolv,svdw,ub,xvdw)
      DO k = lb(3), ub(3)
         DO j = lb(2), ub(2)
            DO i = lb(1), ub(1)
               rho_scaled = MAX(rho_elec%array(i, j, k)/sccs_control%nle_nc_k, nle_n_min)
               x = LOG(rho_scaled)
               xvdw%array(i, j, k) = x
               svdw%array(i, j, k) = sccs_nle_shape_value(x, sccs_control%nle_sigma_k)
               ssolv%array(i, j, k) = svdw%array(i, j, k)
               deps_elec%array(i, j, k) = 0.0_dp
            END DO
         END DO
      END DO
!$OMP END PARALLEL DO

      IF (sccs_control%nle_r_ion > 0.0_dp) THEN
         CALL pw_copy(ssolv, sion)
!$OMP    PARALLEL DO DEFAULT(NONE) PRIVATE(i,j,k) SHARED(lb,sion,ub)
         DO k = lb(3), ub(3)
            DO j = lb(2), ub(2)
               DO i = lb(1), ub(1)
                  sion%array(i, j, k) = 1.0_dp - sion%array(i, j, k)
               END DO
            END DO
         END DO
!$OMP    END PARALLEL DO
         CALL sccs_nle_apply_kernel(sion, kernel_ion, pw_pool)
         CALL sccs_nle_shape_from_field(sion, xion, sion, sccs_control, lb, ub)
      ELSE
         CALL pw_copy(ssolv, sion)
         CALL pw_copy(xvdw, xion)
      END IF

      IF (sccs_control%nle_r_solv > 0.0_dp) THEN
!$OMP    PARALLEL DO DEFAULT(NONE) PRIVATE(i,j,k) SHARED(lb,ssolv,ub)
         DO k = lb(3), ub(3)
            DO j = lb(2), ub(2)
               DO i = lb(1), ub(1)
                  ssolv%array(i, j, k) = 1.0_dp - ssolv%array(i, j, k)
               END DO
            END DO
         END DO
!$OMP    END PARALLEL DO
         CALL sccs_nle_apply_kernel(ssolv, kernel_solv, pw_pool)
         CALL sccs_nle_shape_from_field(ssolv, xsolv, ssolv, sccs_control, lb, ub)
      ELSE
         CALL pw_copy(xvdw, xsolv)
      END IF

      IF (sccs_control%nle_r_cav > 0.0_dp) THEN
         CALL pw_copy(ssolv, scav)
         CALL sccs_nle_apply_kernel(scav, kernel_cav, pw_pool)
         CALL sccs_nle_shape_from_field(scav, xcav, scav, sccs_control, lb, ub)
!$OMP    PARALLEL DO DEFAULT(NONE) PRIVATE(i,j,k) SHARED(lb,scav,ub)
         DO k = lb(3), ub(3)
            DO j = lb(2), ub(2)
               DO i = lb(1), ub(1)
                  scav%array(i, j, k) = 1.0_dp - scav%array(i, j, k)
               END DO
            END DO
         END DO
!$OMP    END PARALLEL DO
      ELSE
         CALL pw_copy(ssolv, scav)
         CALL pw_copy(xsolv, xcav)
      END IF

      IF (sccs_control%nle_r_diel > 0.0_dp) THEN
         CALL pw_copy(ssolv, sdiel)
         CALL sccs_nle_apply_kernel(sdiel, kernel_diel, pw_pool)
         CALL sccs_nle_shape_from_field(sdiel, xdiel, sdiel, sccs_control, lb, ub)
!$OMP    PARALLEL DO DEFAULT(NONE) PRIVATE(i,j,k) SHARED(lb,sdiel,ub)
         DO k = lb(3), ub(3)
            DO j = lb(2), ub(2)
               DO i = lb(1), ub(1)
                  sdiel%array(i, j, k) = 1.0_dp - sdiel%array(i, j, k)
               END DO
            END DO
         END DO
!$OMP    END PARALLEL DO
      ELSE
         CALL pw_copy(ssolv, sdiel)
         CALL pw_copy(xsolv, xdiel)
      END IF

!$OMP PARALLEL DO DEFAULT(NONE) PRIVATE(i,j,k) &
!$OMP             SHARED(eps_elec,lb,sccs_control,sdiel,sdiel_eff,ub)
      DO k = lb(3), ub(3)
         DO j = lb(2), ub(2)
            DO i = lb(1), ub(1)
               sdiel_eff%array(i, j, k) = sdiel%array(i, j, k)
               eps_elec%array(i, j, k) = 1.0_dp + &
                                         (sccs_control%nle_eb_k - 1.0_dp)*sdiel_eff%array(i, j, k)
            END DO
         END DO
      END DO
!$OMP END PARALLEL DO

   END SUBROUTINE sccs_nle_build_shapes

! **************************************************************************************************
!> \brief Add the dielectric shape derivative to the nonlinear electrolyte Kohn-Sham potential.
!> \param v_sccs ...
!> \param lambda_diel ...
!> \param xvdw ...
!> \param xsolv ...
!> \param xdiel ...
!> \param sdiel ...
!> \param svdw ...
!> \param sccs_control ...
!> \param pw_pool ...
!> \param dvol ...
!> \param lb ...
!> \param ub ...
! **************************************************************************************************
   SUBROUTINE sccs_nle_add_dielectric_potential(v_sccs, lambda_diel, xvdw, xsolv, xdiel, &
                                                sdiel, svdw, sccs_control, kernel_diel, kernel_solv, &
                                                pw_pool, dvol, lb, ub)

      TYPE(pw_r3d_rs_type), INTENT(INOUT)                :: v_sccs
      TYPE(pw_r3d_rs_type), INTENT(IN)                   :: lambda_diel, xvdw, xsolv, xdiel, sdiel, &
                                                            svdw
      TYPE(sccs_control_type), POINTER                   :: sccs_control
      TYPE(pw_c1d_gs_type), INTENT(IN)                   :: kernel_diel, kernel_solv
      TYPE(pw_pool_type), POINTER                        :: pw_pool
      REAL(KIND=dp), INTENT(IN)                          :: dvol
      INTEGER, DIMENSION(3), INTENT(IN)                  :: lb, ub

      TYPE(pw_r3d_rs_type)                               :: work

      CALL pw_pool%create_pw(work)
      CALL pw_copy(lambda_diel, work)
      IF (sccs_control%nle_r_diel > 0.0_dp) THEN
         CALL sccs_nle_shape_potential(work, xdiel, sccs_control, lb, ub)
         CALL sccs_nle_apply_kernel(work, kernel_diel, pw_pool, transpose=.TRUE.)
         CALL pw_scale(work, -1.0_dp)
      END IF
      CALL sccs_nle_backprop_ssolv(work, xsolv, sccs_control, kernel_solv, pw_pool, lb, ub)
      CALL sccs_nle_backprop_base(v_sccs, work, xvdw, sccs_control, dvol, lb, ub)
      CALL pw_pool%give_back_pw(work)

   END SUBROUTINE sccs_nle_add_dielectric_potential

! **************************************************************************************************
!> \brief Add the ionic shape derivative to the nonlinear electrolyte Kohn-Sham potential.
!> \param v_sccs ...
!> \param lambda_ion ...
!> \param xvdw ...
!> \param xion ...
!> \param sccs_control ...
!> \param pw_pool ...
!> \param dvol ...
!> \param lb ...
!> \param ub ...
! **************************************************************************************************
   SUBROUTINE sccs_nle_add_ion_potential(v_sccs, lambda_ion, xvdw, xion, sccs_control, &
                                         kernel_ion, pw_pool, dvol, lb, ub)

      TYPE(pw_r3d_rs_type), INTENT(INOUT)                :: v_sccs
      TYPE(pw_r3d_rs_type), INTENT(IN)                   :: lambda_ion, xvdw, xion
      TYPE(sccs_control_type), POINTER                   :: sccs_control
      TYPE(pw_c1d_gs_type), INTENT(IN)                   :: kernel_ion
      TYPE(pw_pool_type), POINTER                        :: pw_pool
      REAL(KIND=dp), INTENT(IN)                          :: dvol
      INTEGER, DIMENSION(3), INTENT(IN)                  :: lb, ub

      TYPE(pw_r3d_rs_type)                               :: work

      CALL pw_pool%create_pw(work)
      CALL pw_copy(lambda_ion, work)
      IF (sccs_control%nle_r_ion > 0.0_dp) THEN
         CALL sccs_nle_shape_potential(work, xion, sccs_control, lb, ub)
         CALL sccs_nle_apply_kernel(work, kernel_ion, pw_pool, transpose=.TRUE.)
         CALL pw_scale(work, -1.0_dp)
      END IF
      CALL sccs_nle_backprop_base(v_sccs, work, xvdw, sccs_control, dvol, lb, ub)
      CALL pw_pool%give_back_pw(work)

   END SUBROUTINE sccs_nle_add_ion_potential

! **************************************************************************************************
!> \brief Add the cavitation shape derivative to the nonlinear electrolyte Kohn-Sham potential.
!> \param v_sccs ...
!> \param lambda_cav ...
!> \param xvdw ...
!> \param xsolv ...
!> \param xcav ...
!> \param sccs_control ...
!> \param pw_pool ...
!> \param dvol ...
!> \param lb ...
!> \param ub ...
! **************************************************************************************************
   SUBROUTINE sccs_nle_add_cavity_potential(v_sccs, lambda_cav, xvdw, xsolv, xcav, &
                                            sccs_control, kernel_cav, kernel_solv, pw_pool, &
                                            dvol, lb, ub)

      TYPE(pw_r3d_rs_type), INTENT(INOUT)                :: v_sccs
      TYPE(pw_r3d_rs_type), INTENT(IN)                   :: lambda_cav, xvdw, xsolv, xcav
      TYPE(sccs_control_type), POINTER                   :: sccs_control
      TYPE(pw_c1d_gs_type), INTENT(IN)                   :: kernel_cav, kernel_solv
      TYPE(pw_pool_type), POINTER                        :: pw_pool
      REAL(KIND=dp), INTENT(IN)                          :: dvol
      INTEGER, DIMENSION(3), INTENT(IN)                  :: lb, ub

      TYPE(pw_r3d_rs_type)                               :: work

      CALL pw_pool%create_pw(work)
      CALL pw_copy(lambda_cav, work)
      IF (sccs_control%nle_r_cav > 0.0_dp) THEN
         CALL sccs_nle_shape_potential(work, xcav, sccs_control, lb, ub)
         CALL sccs_nle_apply_kernel(work, kernel_cav, pw_pool, transpose=.TRUE.)
         CALL pw_scale(work, -1.0_dp)
      END IF
      CALL sccs_nle_backprop_ssolv(work, xsolv, sccs_control, kernel_solv, pw_pool, lb, ub)
      CALL sccs_nle_backprop_base(v_sccs, work, xvdw, sccs_control, dvol, lb, ub)
      CALL pw_pool%give_back_pw(work)

   END SUBROUTINE sccs_nle_add_cavity_potential

! **************************************************************************************************
!> \brief Back-propagate a response field from S_solv to S_vdW.
!> \param response ...
!> \param xsolv ...
!> \param sccs_control ...
!> \param pw_pool ...
!> \param lb ...
!> \param ub ...
   ! **************************************************************************************************
   SUBROUTINE sccs_nle_backprop_ssolv(response, xsolv, sccs_control, kernel_solv, pw_pool, lb, ub)

      TYPE(pw_r3d_rs_type), INTENT(INOUT)                :: response
      TYPE(pw_r3d_rs_type), INTENT(IN)                   :: xsolv
      TYPE(sccs_control_type), POINTER                   :: sccs_control
      TYPE(pw_c1d_gs_type), INTENT(IN)                   :: kernel_solv
      TYPE(pw_pool_type), POINTER                        :: pw_pool
      INTEGER, DIMENSION(3), INTENT(IN)                  :: lb, ub

      IF (sccs_control%nle_r_solv > 0.0_dp) THEN
         CALL sccs_nle_shape_potential(response, xsolv, sccs_control, lb, ub)
         CALL sccs_nle_apply_kernel(response, kernel_solv, pw_pool, transpose=.TRUE.)
         CALL pw_scale(response, -1.0_dp)
      END IF

   END SUBROUTINE sccs_nle_backprop_ssolv

! **************************************************************************************************
!> \brief Back-propagate a response field from S_vdW to electron density.
!> \param v_sccs ...
!> \param response ...
!> \param xvdw ...
!> \param sccs_control ...
!> \param dvol ...
!> \param lb ...
!> \param ub ...
   ! **************************************************************************************************
   SUBROUTINE sccs_nle_backprop_base(v_sccs, response, xvdw, sccs_control, dvol, lb, ub)

      TYPE(pw_r3d_rs_type), INTENT(INOUT)                :: v_sccs, response
      TYPE(pw_r3d_rs_type), INTENT(IN)                   :: xvdw
      TYPE(sccs_control_type), POINTER                   :: sccs_control
      REAL(KIND=dp), INTENT(IN)                          :: dvol
      INTEGER, DIMENSION(3), INTENT(IN)                  :: lb, ub

      INTEGER                                            :: i, j, k

      CALL sccs_nle_shape_potential(response, xvdw, sccs_control, lb, ub)
      CALL pw_scale(response, 1.0_dp/sccs_control%nle_nc_k)

!$OMP PARALLEL DO DEFAULT(NONE) PRIVATE(i,j,k) SHARED(dvol,lb,response,ub,v_sccs)
      DO k = lb(3), ub(3)
         DO j = lb(2), ub(2)
            DO i = lb(1), ub(1)
               v_sccs%array(i, j, k) = v_sccs%array(i, j, k) + dvol*response%array(i, j, k)
            END DO
         END DO
      END DO
!$OMP END PARALLEL DO

   END SUBROUTINE sccs_nle_backprop_base

! **************************************************************************************************
!> \brief Build the cavitation surface/volume derivative with respect to S_cav.
!> \param lambda_cav ...
!> \param grad_scav ...
!> \param norm_grad_scav ...
!> \param sccs_control ...
!> \param pw_pool ...
!> \param pw_env ...
!> \param input ...
!> \param lb ...
!> \param ub ...
! **************************************************************************************************
   SUBROUTINE sccs_nle_cavitation_lambda(lambda_cav, grad_scav, norm_grad_scav, sccs_control, &
                                         pw_pool, pw_env, input, lb, ub)

      TYPE(pw_r3d_rs_type), INTENT(INOUT)                :: lambda_cav
      TYPE(pw_r3d_rs_type), DIMENSION(3), INTENT(IN)     :: grad_scav
      TYPE(pw_r3d_rs_type), INTENT(IN)                   :: norm_grad_scav
      TYPE(sccs_control_type), POINTER                   :: sccs_control
      TYPE(pw_pool_type), POINTER                        :: pw_pool
      TYPE(pw_env_type), POINTER                         :: pw_env
      TYPE(section_vals_type), POINTER                   :: input
      INTEGER, DIMENSION(3), INTENT(IN)                  :: lb, ub

      INTEGER                                            :: i, idir, j, k
      REAL(KIND=dp)                                      :: coeff_surface, denom
      TYPE(pw_r3d_rs_type)                               :: unit_grad
      TYPE(pw_r3d_rs_type), DIMENSION(3)                 :: dunit_grad

      coeff_surface = sccs_control%gamma_solvent + sccs_control%alpha_solvent
      CALL pw_zero(lambda_cav)
      DO idir = 1, 3
         CALL pw_pool%create_pw(dunit_grad(idir))
      END DO
      CALL pw_pool%create_pw(unit_grad)

      IF (ABS(coeff_surface) > 0.0_dp) THEN
         DO idir = 1, 3
!$OMP       PARALLEL DO DEFAULT(NONE) PRIVATE(denom,i,j,k) &
!$OMP                   SHARED(grad_scav,idir,lb,norm_grad_scav,ub,unit_grad)
            DO k = lb(3), ub(3)
               DO j = lb(2), ub(2)
                  DO i = lb(1), ub(1)
                     denom = norm_grad_scav%array(i, j, k)
                     IF (denom > 1.0E-12_dp) THEN
                        unit_grad%array(i, j, k) = grad_scav(idir)%array(i, j, k)/denom
                     ELSE
                        unit_grad%array(i, j, k) = 0.0_dp
                     END IF
                  END DO
               END DO
            END DO
!$OMP       END PARALLEL DO
            CALL derive(unit_grad, dunit_grad, sccs_derivative_fft, pw_env, input)
!$OMP       PARALLEL DO DEFAULT(NONE) PRIVATE(i,j,k) &
!$OMP                   SHARED(coeff_surface,dunit_grad,idir,lambda_cav,lb,ub)
            DO k = lb(3), ub(3)
               DO j = lb(2), ub(2)
                  DO i = lb(1), ub(1)
                     lambda_cav%array(i, j, k) = lambda_cav%array(i, j, k) - &
                                                 coeff_surface*dunit_grad(idir)%array(i, j, k)
                  END DO
               END DO
            END DO
!$OMP       END PARALLEL DO
         END DO
      END IF

      IF (ABS(sccs_control%beta_solvent) > 0.0_dp) THEN
!$OMP    PARALLEL DO DEFAULT(NONE) PRIVATE(i,j,k) SHARED(lambda_cav,lb,sccs_control,ub)
         DO k = lb(3), ub(3)
            DO j = lb(2), ub(2)
               DO i = lb(1), ub(1)
                  lambda_cav%array(i, j, k) = lambda_cav%array(i, j, k) - sccs_control%beta_solvent
               END DO
            END DO
         END DO
!$OMP    END PARALLEL DO
      END IF

      CALL pw_pool%give_back_pw(unit_grad)
      DO idir = 1, 3
         CALL pw_pool%give_back_pw(dunit_grad(idir))
      END DO

   END SUBROUTINE sccs_nle_cavitation_lambda

! **************************************************************************************************
!> \brief Convert a nonlinear electrolyte shape response dE/dS to dE/dexp(X).
!> \param response ...
!> \param xfield ...
!> \param sccs_control ...
!> \param lb ...
!> \param ub ...
! **************************************************************************************************
   SUBROUTINE sccs_nle_shape_potential(response, xfield, sccs_control, lb, ub)

      TYPE(pw_r3d_rs_type), INTENT(INOUT)                :: response
      TYPE(pw_r3d_rs_type), INTENT(IN)                   :: xfield
      TYPE(sccs_control_type), POINTER                   :: sccs_control
      INTEGER, DIMENSION(3), INTENT(IN)                  :: lb, ub

      INTEGER                                            :: i, j, k
      REAL(KIND=dp)                                      :: x

!$OMP PARALLEL DO DEFAULT(NONE) PRIVATE(i,j,k,x) SHARED(lb,response,sccs_control,ub,xfield)
      DO k = lb(3), ub(3)
         DO j = lb(2), ub(2)
            DO i = lb(1), ub(1)
               x = xfield%array(i, j, k)
               response%array(i, j, k) = response%array(i, j, k)* &
                                         sccs_nle_shape_prime(x, sccs_control%nle_sigma_k)*EXP(-x)
            END DO
         END DO
      END DO
!$OMP END PARALLEL DO

   END SUBROUTINE sccs_nle_shape_potential

! **************************************************************************************************
!> \brief Build X=log(field) and the nonlinear electrolyte erfc shape from a convolution field.
!> \param field ...
!> \param xfield ...
!> \param shape ...
!> \param sccs_control ...
!> \param lb ...
!> \param ub ...
! **************************************************************************************************
   SUBROUTINE sccs_nle_shape_from_field(field, xfield, shape, sccs_control, lb, ub)

      TYPE(pw_r3d_rs_type), INTENT(IN)                   :: field
      TYPE(pw_r3d_rs_type), INTENT(INOUT)                :: xfield, shape
      TYPE(sccs_control_type), POINTER                   :: sccs_control
      INTEGER, DIMENSION(3), INTENT(IN)                  :: lb, ub

      INTEGER                                            :: i, j, k
      REAL(KIND=dp)                                      :: x

!$OMP PARALLEL DO DEFAULT(NONE) PRIVATE(i,j,k,x) SHARED(field,lb,sccs_control,shape,ub,xfield)
      DO k = lb(3), ub(3)
         DO j = lb(2), ub(2)
            DO i = lb(1), ub(1)
               x = LOG(MAX(field%array(i, j, k), nle_n_min))
               xfield%array(i, j, k) = x
               shape%array(i, j, k) = sccs_nle_shape_value(x, sccs_control%nle_sigma_k)
            END DO
         END DO
      END DO
!$OMP END PARALLEL DO

   END SUBROUTINE sccs_nle_shape_from_field

! **************************************************************************************************
!> \brief Prepare a nonlinear electrolyte F_exp convolution kernel in reciprocal space.
!> \param kernel_g ...
!> \param radius ...
!> \param sccs_control ...
!> \param pw_pool ...
! **************************************************************************************************
   SUBROUTINE sccs_nle_prepare_exp_kernel(kernel_g, radius, sccs_control, pw_pool)

      TYPE(pw_c1d_gs_type), INTENT(INOUT)                :: kernel_g
      REAL(KIND=dp), INTENT(IN)                          :: radius
      TYPE(sccs_control_type), POINTER                   :: sccs_control
      TYPE(pw_pool_type), POINTER                        :: pw_pool

      REAL(KIND=dp)                                      :: factor, sigma
      TYPE(pw_r3d_rs_type)                               :: kernel_rs

      sigma = sccs_control%nle_a_k/sccs_control%nle_sigma_k
      CALL pw_pool%create_pw(kernel_rs)
      CALL sccs_nle_fill_exp_kernel(kernel_rs, radius, sigma)
      CALL pw_transfer(kernel_rs, kernel_g)
      factor = REAL(kernel_rs%pw_grid%ngpts, KIND=dp)*kernel_rs%pw_grid%dvol
      CALL pw_scale(kernel_g, factor)
      CALL pw_pool%give_back_pw(kernel_rs)

   END SUBROUTINE sccs_nle_prepare_exp_kernel

! **************************************************************************************************
!> \brief Apply a prepared nonlinear electrolyte F_exp kernel to a real-space field.
!> \param field ...
!> \param kernel_g ...
!> \param pw_pool ...
!> \param transpose ...
! **************************************************************************************************
   SUBROUTINE sccs_nle_apply_kernel(field, kernel_g, pw_pool, transpose)

      TYPE(pw_r3d_rs_type), INTENT(INOUT)                :: field
      TYPE(pw_c1d_gs_type), INTENT(IN)                   :: kernel_g
      TYPE(pw_pool_type), POINTER                        :: pw_pool
      LOGICAL, INTENT(IN), OPTIONAL                      :: transpose

      INTEGER                                            :: ig
      LOGICAL                                            :: do_transpose
      TYPE(pw_c1d_gs_type)                               :: field_g

      do_transpose = .FALSE.
      IF (PRESENT(transpose)) do_transpose = transpose

      CALL pw_pool%create_pw(field_g)
      CALL pw_transfer(field, field_g)

      IF (do_transpose) THEN
!$OMP    PARALLEL DO DEFAULT(NONE) PRIVATE(ig) SHARED(field_g,kernel_g)
         DO ig = 1, SIZE(field_g%array)
            field_g%array(ig) = field_g%array(ig)*CONJG(kernel_g%array(ig))
         END DO
!$OMP    END PARALLEL DO
      ELSE
!$OMP    PARALLEL DO DEFAULT(NONE) PRIVATE(ig) SHARED(field_g,kernel_g)
         DO ig = 1, SIZE(field_g%array)
            field_g%array(ig) = field_g%array(ig)*kernel_g%array(ig)
         END DO
!$OMP    END PARALLEL DO
      END IF

      CALL pw_transfer(field_g, field)

      CALL pw_pool%give_back_pw(field_g)

   END SUBROUTINE sccs_nle_apply_kernel

! **************************************************************************************************
!> \brief Fill the real-space nonlinear electrolyte F_exp kernel on the periodic grid.
!> \param kernel ...
!> \param radius ...
!> \param sigma ...
! **************************************************************************************************
   SUBROUTINE sccs_nle_fill_exp_kernel(kernel, radius, sigma)

      TYPE(pw_r3d_rs_type), INTENT(INOUT)                :: kernel
      REAL(KIND=dp), INTENT(IN)                          :: radius, sigma

      INTEGER                                            :: i, ii, j, jj, k, kk
      INTEGER, DIMENSION(2, 3)                           :: bounds
      INTEGER, DIMENSION(3)                              :: lb, npts, ub
      REAL(KIND=dp)                                      :: arg, denom, r
      REAL(KIND=dp), DIMENSION(3)                        :: rvec

      npts(:) = kernel%pw_grid%npts(:)
      bounds(:, :) = kernel%pw_grid%bounds(:, :)
      lb(:) = kernel%pw_grid%bounds_local(1, :)
      ub(:) = kernel%pw_grid%bounds_local(2, :)

!$OMP PARALLEL DO DEFAULT(NONE) PRIVATE(arg,denom,i,ii,j,jj,k,kk,r,rvec) &
!$OMP             SHARED(bounds,kernel,lb,npts,radius,sigma,ub)
      DO k = lb(3), ub(3)
         kk = MODULO(k - bounds(1, 3), npts(3))
         IF ((kk > 0) .AND. (2*kk >= npts(3))) kk = kk - npts(3)
         DO j = lb(2), ub(2)
            jj = MODULO(j - bounds(1, 2), npts(2))
            IF ((jj > 0) .AND. (2*jj >= npts(2))) jj = jj - npts(2)
            DO i = lb(1), ub(1)
               ii = MODULO(i - bounds(1, 1), npts(1))
               IF ((ii > 0) .AND. (2*ii >= npts(1))) ii = ii - npts(1)
               rvec(:) = REAL(ii, KIND=dp)*kernel%pw_grid%dh(:, 1) + &
                         REAL(jj, KIND=dp)*kernel%pw_grid%dh(:, 2) + &
                         REAL(kk, KIND=dp)*kernel%pw_grid%dh(:, 3)
               r = SQRT(SUM(rvec(:)**2))
               arg = MAX(-80.0_dp, MIN(80.0_dp, (r - radius)/sigma))
               denom = 1.0_dp/nle_exp_cutoff + EXP(arg)
               kernel%array(i, j, k) = 1.0_dp/denom/(fourpi*sigma**3)* &
                                       4.0_dp/(2.0_dp + radius/sigma)
            END DO
         END DO
      END DO
!$OMP END PARALLEL DO

   END SUBROUTINE sccs_nle_fill_exp_kernel

! **************************************************************************************************
!> \brief nonlinear electrolyte erfc shape value.
! **************************************************************************************************
   PURE FUNCTION sccs_nle_shape_value(x, sigma) RESULT(value)

      REAL(KIND=dp), INTENT(IN)                          :: x, sigma
      REAL(KIND=dp)                                      :: value

      value = 0.5_dp*ERFC(x/(SQRT(2.0_dp)*sigma))
      value = MAX(0.0_dp, MIN(1.0_dp, value))

   END FUNCTION sccs_nle_shape_value

! **************************************************************************************************
!> \brief nonlinear electrolyte erfc shape derivative with respect to X.
! **************************************************************************************************
   PURE FUNCTION sccs_nle_shape_prime(x, sigma) RESULT(value)

      REAL(KIND=dp), INTENT(IN)                          :: x, sigma
      REAL(KIND=dp)                                      :: value

      value = -EXP(-0.5_dp*(x/sigma)**2)/(SQRT(twopi)*sigma)

   END FUNCTION sccs_nle_shape_prime

! **************************************************************************************************
!> \brief      Calculate the smoothed dielectric function of Andreussi et al.
!> \param rho_elec ...
!> \param eps_elec ...
!> \param deps_elec ...
!> \param epsilon_solvent ...
!> \param rho_max ...
!> \param rho_min ...
!> \par History:
!>      - Creation (16.10.2013,MK)
!>      - Finite difference of isosurfaces implemented (21.12.2013,MK)
!> \author     Matthias Krack (MK)
!> \version    1.1
! **************************************************************************************************
   SUBROUTINE andreussi(rho_elec, eps_elec, deps_elec, epsilon_solvent, rho_max, &
                        rho_min)

      TYPE(pw_r3d_rs_type), INTENT(IN)                   :: rho_elec, eps_elec, deps_elec
      REAL(KIND=dp), INTENT(IN)                          :: epsilon_solvent, rho_max, rho_min

      CHARACTER(LEN=*), PARAMETER                        :: routineN = 'andreussi'
      REAL(KIND=dp), PARAMETER                           :: rhotol = 1.0E-12_dp

      INTEGER                                            :: handle, i, j, k
      INTEGER, DIMENSION(3)                              :: lb, ub
      REAL(KIND=dp)                                      :: diff, dq, dt, f, ln_rho_max, ln_rho_min, &
                                                            q, rho, t, x, y

      CALL timeset(routineN, handle)

      f = LOG(epsilon_solvent)/twopi
      diff = rho_max - rho_min
      IF (diff < SQRT(rhotol)) CPABORT("SCCS: Difference between rho(min) and rho(max) is too small")
      IF (rho_min >= rhotol) THEN
         ln_rho_max = LOG(rho_max)
         ln_rho_min = LOG(rho_min)
         q = twopi/(ln_rho_max - ln_rho_min)
         dq = -f*q
      END IF

      lb(1:3) = rho_elec%pw_grid%bounds_local(1, 1:3)
      ub(1:3) = rho_elec%pw_grid%bounds_local(2, 1:3)

      ! Calculate the dielectric function and its derivative
!$OMP PARALLEL DO DEFAULT(NONE) &
!$OMP             PRIVATE(dt,i,j,k,rho,t,x,y) &
!$OMP             SHARED(deps_elec,dq,eps_elec,epsilon_solvent,f,lb,ub) &
!$OMP             SHARED(ln_rho_max,rho_elec,q,rho_max,rho_min)
      DO k = lb(3), ub(3)
         DO j = lb(2), ub(2)
            DO i = lb(1), ub(1)
               rho = rho_elec%array(i, j, k)
               IF (rho < rho_min) THEN
                  eps_elec%array(i, j, k) = epsilon_solvent
                  deps_elec%array(i, j, k) = 0.0_dp
               ELSE IF (rho <= rho_max) THEN
                  x = LOG(rho)
                  y = q*(ln_rho_max - x)
                  t = f*(y - SIN(y))
                  eps_elec%array(i, j, k) = EXP(t)
                  dt = dq*(1.0_dp - COS(y))
                  deps_elec%array(i, j, k) = eps_elec%array(i, j, k)*dt/rho
               ELSE
                  eps_elec%array(i, j, k) = 1.0_dp
                  deps_elec%array(i, j, k) = 0.0_dp
               END IF
            END DO
         END DO
      END DO
!$OMP END PARALLEL DO

      CALL timestop(handle)

   END SUBROUTINE andreussi

! **************************************************************************************************
!> \brief      Calculate the smoothed dielectric function of Fattebert and Gygi
!> \param rho_elec ...
!> \param eps_elec ...
!> \param deps_elec ...
!> \param epsilon_solvent ...
!> \param beta ...
!> \param rho_zero ...
!> \par History:
!>      - Creation (15.10.2013,MK)
!> \author     Matthias Krack (MK)
!> \version    1.0
! **************************************************************************************************
   SUBROUTINE fattebert_gygi(rho_elec, eps_elec, deps_elec, epsilon_solvent, beta, &
                             rho_zero)

      TYPE(pw_r3d_rs_type), INTENT(IN)                   :: rho_elec, eps_elec, deps_elec
      REAL(KIND=dp), INTENT(IN)                          :: epsilon_solvent, beta, rho_zero

      CHARACTER(LEN=*), PARAMETER                        :: routineN = 'fattebert_gygi'
      REAL(KIND=dp), PARAMETER                           :: rhotol = 1.0E-12_dp

      INTEGER                                            :: handle, i, j, k
      INTEGER, DIMENSION(3)                              :: lb, ub
      REAL(KIND=dp)                                      :: df, f, p, q, rho, s, t, twobeta

      CALL timeset(routineN, handle)

      df = (1.0_dp - epsilon_solvent)/rho_zero
      f = 0.5_dp*(epsilon_solvent - 1.0_dp)
      q = 1.0_dp/rho_zero
      twobeta = 2.0_dp*beta

      lb(1:3) = rho_elec%pw_grid%bounds_local(1, 1:3)
      ub(1:3) = rho_elec%pw_grid%bounds_local(2, 1:3)

      ! Calculate the smoothed dielectric function and its derivative
!$OMP PARALLEL DO DEFAULT(NONE) &
!$OMP             PRIVATE(i,j,k,p,rho,s,t) &
!$OMP             SHARED(df,deps_elec,eps_elec,epsilon_solvent,f,lb,ub) &
!$OMP             SHARED(q,rho_elec,twobeta)
      DO k = lb(3), ub(3)
         DO j = lb(2), ub(2)
            DO i = lb(1), ub(1)
               rho = rho_elec%array(i, j, k)
               IF (rho < rhotol) THEN
                  eps_elec%array(i, j, k) = epsilon_solvent
                  deps_elec%array(i, j, k) = 0.0_dp
               ELSE
                  s = rho*q
                  p = s**twobeta
                  t = 1.0_dp/(1.0_dp + p)
                  eps_elec%array(i, j, k) = 1.0_dp + f*(1.0_dp + (1.0_dp - p)*t)
                  deps_elec%array(i, j, k) = df*twobeta*t*t*p/s
               END IF
            END DO
         END DO
      END DO
!$OMP END PARALLEL DO

      CALL timestop(handle)

   END SUBROUTINE fattebert_gygi

! **************************************************************************************************
!> \brief      Calculate the smoothed dielectric function of solvent aware algorithm
!> \param rho_elec ...
!> \param eps_elec ...
!> \param deps_elec ...
!> \param epsilon_solvent ...
!> \param rho_max ...
!> \param rho_min ...
!> \param s ...
!> \param d_s_rhoel ...
!> \author     Ziwei Chai
!> \version    1.0
! **************************************************************************************************
   SUBROUTINE sa_andreussi(rho_elec, eps_elec, deps_elec, epsilon_solvent, rho_max, &
                           rho_min, s, d_s_rhoel)

      TYPE(pw_r3d_rs_type), INTENT(IN)                   :: rho_elec, eps_elec, deps_elec
      REAL(KIND=dp), INTENT(IN)                          :: epsilon_solvent, rho_max, rho_min
      TYPE(pw_r3d_rs_type), INTENT(IN)                   :: s, d_s_rhoel

      CHARACTER(LEN=*), PARAMETER                        :: routineN = 'sa_andreussi'
      REAL(KIND=dp), PARAMETER                           :: rhotol = 1.0E-12_dp

      INTEGER                                            :: handle, i, j, k
      INTEGER, DIMENSION(3)                              :: lb, ub
      REAL(KIND=dp)                                      :: diff, dq, dt, f, ln_rho_max, ln_rho_min, &
                                                            q, rho, t, x, y

      CALL timeset(routineN, handle)

      f = LOG(epsilon_solvent)/twopi
      diff = rho_max - rho_min
      IF (diff < SQRT(rhotol)) CPABORT("SCCS: Difference between rho(min) and rho(max) is too small")
      IF (rho_min >= rhotol) THEN
         ln_rho_max = LOG(rho_max)
         ln_rho_min = LOG(rho_min)
         q = twopi/(ln_rho_max - ln_rho_min)
         dq = -f*q
      END IF

      lb(1:3) = rho_elec%pw_grid%bounds_local(1, 1:3)
      ub(1:3) = rho_elec%pw_grid%bounds_local(2, 1:3)

      ! Calculate the dielectric function and its derivative
!$OMP PARALLEL DO DEFAULT(NONE) &
!$OMP             PRIVATE(dt,i,j,k,rho,t,x,y) &
!$OMP             SHARED(deps_elec,dq,eps_elec,epsilon_solvent,f,lb,ub) &
!$OMP             SHARED(ln_rho_max,rho_elec,q,rho_max,rho_min, s, d_s_rhoel)
      DO k = lb(3), ub(3)
         DO j = lb(2), ub(2)
            DO i = lb(1), ub(1)
               rho = rho_elec%array(i, j, k)
               IF (rho < rho_min) THEN
                  eps_elec%array(i, j, k) = epsilon_solvent
                  s%array(i, j, k) = 0.0_dp
                  deps_elec%array(i, j, k) = 0.0_dp
                  d_s_rhoel%array(i, j, k) = 0.0_dp
               ELSE IF (rho <= rho_max) THEN
                  x = LOG(rho)
                  y = q*(ln_rho_max - x)
                  t = f*(y - SIN(y))
                  s%array(i, j, k) = 1 - t/(f*twopi)
                  eps_elec%array(i, j, k) = EXP(t)
                  dt = dq*(1.0_dp - COS(y))
                  deps_elec%array(i, j, k) = eps_elec%array(i, j, k)*dt/rho
                  d_s_rhoel%array(i, j, k) = deps_elec%array(i, j, k)/(-eps_elec%array(i, j, k)*LOG(epsilon_solvent))
               ELSE
                  eps_elec%array(i, j, k) = 1.0_dp
                  deps_elec%array(i, j, k) = 0.0_dp
                  s%array(i, j, k) = 1.0_dp
                  d_s_rhoel%array(i, j, k) = 0.0_dp
               END IF
            END DO
         END DO
      END DO
!$OMP END PARALLEL DO

      CALL timestop(handle)

   END SUBROUTINE sa_andreussi

! **************************************************************************************************
!> \brief      Build the numerical derivative of a function on realspace grid
!> \param f ...
!> \param df ...
!> \param method ...
!> \param pw_env ...
!> \param input ...
!> \par History:
!>      - Creation (15.11.2013,MK)
!> \author     Matthias Krack (MK)
!> \version    1.0
! **************************************************************************************************
   SUBROUTINE derive(f, df, method, pw_env, input)

      TYPE(pw_r3d_rs_type), INTENT(IN)                   :: f
      TYPE(pw_r3d_rs_type), DIMENSION(3), INTENT(INOUT)  :: df
      INTEGER, INTENT(IN)                                :: method
      TYPE(pw_env_type), POINTER                         :: pw_env
      TYPE(section_vals_type), POINTER                   :: input

      CHARACTER(LEN=*), PARAMETER                        :: routineN = 'derive'

      INTEGER                                            :: border_points, handle, i
      INTEGER, DIMENSION(3)                              :: lb, n, ub
      TYPE(pw_c1d_gs_type), DIMENSION(2)                 :: work_g1d
      TYPE(pw_pool_type), POINTER                        :: auxbas_pw_pool
      TYPE(realspace_grid_desc_type), POINTER            :: rs_desc
      TYPE(realspace_grid_input_type)                    :: input_settings
      TYPE(realspace_grid_type), POINTER                 :: rs_grid
      TYPE(section_vals_type), POINTER                   :: rs_grid_section

      CALL timeset(routineN, handle)

      CPASSERT(ASSOCIATED(pw_env))

      ! Perform method specific setup
      SELECT CASE (method)
      CASE (sccs_derivative_cd3, sccs_derivative_cd5, sccs_derivative_cd7)
         NULLIFY (rs_desc)
         rs_grid_section => section_vals_get_subs_vals(input, "DFT%MGRID%RS_GRID")
         SELECT CASE (method)
         CASE (sccs_derivative_cd3)
            border_points = 1
         CASE (sccs_derivative_cd5)
            border_points = 2
         CASE (sccs_derivative_cd7)
            border_points = 3
         END SELECT
         CALL init_input_type(input_settings, 2*border_points + 1, rs_grid_section, &
                              1, [-1, -1, -1])
         CALL rs_grid_create_descriptor(rs_desc, f%pw_grid, input_settings, &
                                        border_points=border_points)
         ALLOCATE (rs_grid)
         CALL rs_grid_create(rs_grid, rs_desc)
!MK      CALL rs_grid_print(rs_grid, 6)
      CASE (sccs_derivative_fft)
         lb(1:3) = f%pw_grid%bounds_local(1, 1:3)
         ub(1:3) = f%pw_grid%bounds_local(2, 1:3)
         NULLIFY (auxbas_pw_pool)
         CALL pw_env_get(pw_env, auxbas_pw_pool=auxbas_pw_pool)
         ! Get work storage for the 1d grids in g-space (derivative calculation)
         DO i = 1, SIZE(work_g1d)
            CALL auxbas_pw_pool%create_pw(work_g1d(i))
         END DO
      END SELECT

      ! Calculate the derivatives
      SELECT CASE (method)
      CASE (sccs_derivative_cd3)
         CALL derive_fdm_cd3(f, df, rs_grid)
      CASE (sccs_derivative_cd5)
         CALL derive_fdm_cd5(f, df, rs_grid)
      CASE (sccs_derivative_cd7)
         CALL derive_fdm_cd7(f, df, rs_grid)
      CASE (sccs_derivative_fft)
         ! FFT
         CALL pw_transfer(f, work_g1d(1))
         DO i = 1, 3
            n(:) = 0
            n(i) = 1
            CALL pw_copy(work_g1d(1), work_g1d(2))
            CALL pw_derive(work_g1d(2), n(:))
            CALL pw_transfer(work_g1d(2), df(i))
         END DO
      CASE DEFAULT
         CPABORT("Invalid derivative method for SCCS specified")
      END SELECT

      ! Perform method specific cleanup
      SELECT CASE (method)
      CASE (sccs_derivative_cd3, sccs_derivative_cd5, sccs_derivative_cd7)
         CALL rs_grid_release(rs_grid)
         DEALLOCATE (rs_grid)
         CALL rs_grid_release_descriptor(rs_desc)
      CASE (sccs_derivative_fft)
         DO i = 1, SIZE(work_g1d)
            CALL auxbas_pw_pool%give_back_pw(work_g1d(i))
         END DO
      END SELECT

      CALL timestop(handle)

   END SUBROUTINE derive

! **************************************************************************************************
!> \brief      Calculate the finite difference between two isosurfaces of the
!>             electronic density. The smoothed dielectric function of
!>             Andreussi et al. is used as switching function eventually
!>             defining the quantum volume and surface of the cavity.
!> \param rho_elec ...
!> \param norm_drho_elec ...
!> \param dtheta ...
!> \param epsilon_solvent ...
!> \param rho_max ...
!> \param rho_min ...
!> \param delta_rho ...
!> \par History:
!>      - Creation (21.12.2013,MK)
!> \author     Matthias Krack (MK)
!> \version    1.0
! **************************************************************************************************
   SUBROUTINE surface_andreussi(rho_elec, norm_drho_elec, dtheta, &
                                epsilon_solvent, rho_max, rho_min, delta_rho)

      TYPE(pw_r3d_rs_type), INTENT(IN)                   :: rho_elec, norm_drho_elec, dtheta
      REAL(KIND=dp), INTENT(IN)                          :: epsilon_solvent, rho_max, rho_min, &
                                                            delta_rho

      CHARACTER(LEN=*), PARAMETER                        :: routineN = 'surface_andreussi'
      REAL(KIND=dp), PARAMETER                           :: rhotol = 1.0E-12_dp

      INTEGER                                            :: handle, i, j, k, l
      INTEGER, DIMENSION(3)                              :: lb, ub
      REAL(KIND=dp)                                      :: diff, e, eps_elec, f, ln_rho_max, &
                                                            ln_rho_min, q, rho, t, x, y
      REAL(KIND=dp), DIMENSION(2)                        :: theta

      CALL timeset(routineN, handle)

      e = epsilon_solvent - 1.0_dp
      f = LOG(epsilon_solvent)/twopi
      diff = rho_max - rho_min
      IF (diff < SQRT(rhotol)) CPABORT("SCCS: Difference between rho(min) and rho(max) is too small")
      IF (rho_min >= rhotol) THEN
         ln_rho_max = LOG(rho_max)
         ln_rho_min = LOG(rho_min)
         q = twopi/(ln_rho_max - ln_rho_min)
      END IF

      lb(1:3) = rho_elec%pw_grid%bounds_local(1, 1:3)
      ub(1:3) = rho_elec%pw_grid%bounds_local(2, 1:3)

      ! Calculate finite difference between two isosurfaces
!$OMP PARALLEL DO DEFAULT(NONE) &
!$OMP             PRIVATE(eps_elec,i,j,k,l,rho,t,theta,x,y) &
!$OMP             SHARED(delta_rho,dtheta,e,epsilon_solvent,f,lb) &
!$OMP             SHARED(ln_rho_max,norm_drho_elec,rho_elec,q,rho_max,rho_min,ub)
      DO k = lb(3), ub(3)
         DO j = lb(2), ub(2)
            DO i = lb(1), ub(1)
               DO l = 1, 2
                  rho = rho_elec%array(i, j, k) + (REAL(l, KIND=dp) - 1.5_dp)*delta_rho
                  IF (rho < rho_min) THEN
                     eps_elec = epsilon_solvent
                  ELSE IF (rho <= rho_max) THEN
                     x = LOG(rho)
                     y = q*(ln_rho_max - x)
                     t = f*(y - SIN(y))
                     eps_elec = EXP(t)
                  ELSE
                     eps_elec = 1.0_dp
                  END IF
                  theta(l) = (epsilon_solvent - eps_elec)/e
               END DO
               dtheta%array(i, j, k) = (theta(2) - theta(1))*norm_drho_elec%array(i, j, k)/delta_rho
            END DO
         END DO
      END DO
!$OMP END PARALLEL DO

      CALL timestop(handle)

   END SUBROUTINE surface_andreussi

! **************************************************************************************************
!> \brief      Calculate the finite difference between two isosurfaces of the
!>             the electronic density. The smoothed dielectric function of
!>             Fattebert and Gygi is used as switching function eventually
!>             defining the quantum volume and surface of the cavity.
!> \param rho_elec ...
!> \param norm_drho_elec ...
!> \param dtheta ...
!> \param epsilon_solvent ...
!> \param beta ...
!> \param rho_zero ...
!> \param delta_rho ...
!> \par History:
!>      - Creation (21.12.2013,MK)
!> \author     Matthias Krack (MK)
!> \version    1.0
! **************************************************************************************************
   SUBROUTINE surface_fattebert_gygi(rho_elec, norm_drho_elec, dtheta, &
                                     epsilon_solvent, beta, rho_zero, delta_rho)

      TYPE(pw_r3d_rs_type), INTENT(IN)                   :: rho_elec, norm_drho_elec, dtheta
      REAL(KIND=dp), INTENT(IN)                          :: epsilon_solvent, beta, rho_zero, &
                                                            delta_rho

      CHARACTER(LEN=*), PARAMETER :: routineN = 'surface_fattebert_gygi'
      REAL(KIND=dp), PARAMETER                           :: rhotol = 1.0E-12_dp

      INTEGER                                            :: handle, i, j, k, l
      INTEGER, DIMENSION(3)                              :: lb, ub
      REAL(KIND=dp)                                      :: e, eps_elec, f, p, q, rho, s, t, twobeta
      REAL(KIND=dp), DIMENSION(2)                        :: theta

      CALL timeset(routineN, handle)

      e = epsilon_solvent - 1.0_dp
      f = 0.5_dp*e
      q = 1.0_dp/rho_zero
      twobeta = 2.0_dp*beta

      lb(1:3) = rho_elec%pw_grid%bounds_local(1, 1:3)
      ub(1:3) = rho_elec%pw_grid%bounds_local(2, 1:3)

      ! Calculate finite difference between two isosurfaces
!$OMP PARALLEL DO DEFAULT(NONE) &
!$OMP             PRIVATE(eps_elec,i,j,k,l,p,rho,s,t,theta) &
!$OMP             SHARED(delta_rho,dtheta,e,epsilon_solvent,f,lb) &
!$OMP             SHARED(norm_drho_elec,q,rho_elec,twobeta,ub)
      DO k = lb(3), ub(3)
         DO j = lb(2), ub(2)
            DO i = lb(1), ub(1)
               DO l = 1, 2
                  rho = rho_elec%array(i, j, k) + (REAL(l, KIND=dp) - 1.5_dp)*delta_rho
                  IF (rho < rhotol) THEN
                     eps_elec = epsilon_solvent
                  ELSE
                     s = rho*q
                     p = s**twobeta
                     t = 1.0_dp/(1.0_dp + p)
                     eps_elec = 1.0_dp + f*(1.0_dp + (1.0_dp - p)*t)
                  END IF
                  theta(l) = (epsilon_solvent - eps_elec)/e
               END DO
               dtheta%array(i, j, k) = (theta(2) - theta(1))*norm_drho_elec%array(i, j, k)/delta_rho
            END DO
         END DO
      END DO
!$OMP END PARALLEL DO

      CALL timestop(handle)

   END SUBROUTINE surface_fattebert_gygi

! **************************************************************************************************
!> \brief      Print SCCS results
!> \param energy ...
!> \param sccs_control ...
!> \param output_unit ...
!> \par History:
!>      - Creation (11.11.2022,MK)
!> \author     Matthias Krack (MK)
!> \version    1.0
! **************************************************************************************************
   SUBROUTINE print_sccs_results(energy, sccs_control, output_unit)

      TYPE(qs_energy_type), POINTER                      :: energy
      TYPE(sccs_control_type), POINTER                   :: sccs_control
      INTEGER, INTENT(IN)                                :: output_unit

      IF (output_unit > 0) THEN
         CPASSERT(ASSOCIATED(energy))
         CPASSERT(ASSOCIATED(sccs_control))
         WRITE (UNIT=output_unit, FMT="(T3,A,T56,F25.14)") &
            "SCCS| Hartree energy of solute and solvent [Hartree]", energy%sccs_hartree, &
            "SCCS| Hartree energy of the solute only    [Hartree]", energy%hartree
         WRITE (UNIT=output_unit, FMT="(T3,A,T56,F25.14,/,T3,A,T61,F20.3)") &
            "SCCS| Polarisation energy                  [Hartree]", energy%sccs_pol, &
            "SCCS|                                      [kcal/mol]", &
            cp_unit_from_cp2k(energy%sccs_pol, "kcalmol"), &
            "SCCS| Cavitation energy                    [Hartree]", energy%sccs_cav, &
            "SCCS|                                      [kcal/mol]", &
            cp_unit_from_cp2k(energy%sccs_cav, "kcalmol"), &
            "SCCS| Dispersion free energy               [Hartree]", energy%sccs_dis, &
            "SCCS|                                      [kcal/mol]", &
            cp_unit_from_cp2k(energy%sccs_dis, "kcalmol"), &
            "SCCS| Repulsion free energy                [Hartree]", energy%sccs_rep, &
            "SCCS|                                      [kcal/mol]", &
            cp_unit_from_cp2k(energy%sccs_rep, "kcalmol"), &
            "SCCS| Solvation free energy                [Hartree]", energy%sccs_sol, &
            "SCCS|                                      [kcal/mol]", &
            cp_unit_from_cp2k(energy%sccs_sol, "kcalmol")
      END IF

   END SUBROUTINE print_sccs_results

END MODULE qs_sccs
