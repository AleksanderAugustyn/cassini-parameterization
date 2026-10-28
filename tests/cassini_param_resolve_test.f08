!> Node sets and the resolve: c, centre of mass, poles, flags and statuses,
!! the volume against an independent Simpson rule, the t-substituted
!! quadrature through the pinch, and the quadrupole.
program cassini_param_resolve_test
    use precision_utilities_mod, only: ik, rk
    use mathematical_and_physical_constants_mod, only: PI_C
    use cassini_parameterization_workers_mod, only: tables_t, resolve_t, build_tables_s, &
            quadrature_nodes_s, cosine_nodes_s, resolve_s, quadrupole_f, series_on_tables_s, &
            lemniscate_nodes_s, CASSINI_N_GL_HALF, CASSINI_N_R_CHECK, CASSINI_MAX_ORDER, &
            CASSINI_ERROR_R_NONPOSITIVE, CASSINI_ERROR_TOROIDAL
    use cassini_test_fixtures_mod, only: SPHERE, BARE_05, BARE_095, OBLATE, GS, SADDLE1, SD_MIN, &
            SADDLE2, ASYM, DIMPLED, FOLD_BOX, TWO_NECKS, STAR_FAIL, PINCHED, PINCHED_ASYM, &
            SEPARATED, NEAR_PINCH, R_NEG, TOROIDAL, TOROIDAL_ONE
    use cassini_test_reference_mod, only: reference_volume_s, reference_quadrupole_f, &
            reference_gl_volume_s
    use test_utils_mod, only: assert_true, assert_int_eq, assert_abs_close, assert_close, &
            assert_bits_eq, is_zero_f, test_summary
    implicit none

    real(kind = rk), parameter :: FOUR_THIRDS_PI = 4.0_rk * PI_C / 3.0_rk
    type(tables_t) :: quad_tab, check_tab
    real(kind = rk) :: quad_x(2 * CASSINI_N_GL_HALF), quad_w(2 * CASSINI_N_GL_HALF)
    real(kind = rk) :: check_x(CASSINI_N_R_CHECK)

    call quadrature_nodes_s(quad_x, quad_w)
    call build_tables_s(quad_x, CASSINI_MAX_ORDER, quad_tab)
    call cosine_nodes_s(CASSINI_N_R_CHECK, check_x)
    call build_tables_s(check_x, CASSINI_MAX_ORDER, check_tab)

    call test_nodes()
    call test_sphere_and_bare()
    call test_236u()
    call test_volume(SPHERE, 'SPHERE')
    call test_volume(BARE_05, 'BARE_05')
    call test_volume(BARE_095, 'BARE_095')
    call test_volume(OBLATE, 'OBLATE')
    call test_volume(GS, 'GS')
    call test_volume(SADDLE1, 'SADDLE1')
    call test_volume(SD_MIN, 'SD_MIN')
    call test_volume(SADDLE2, 'SADDLE2')
    call test_volume(ASYM, 'ASYM')
    call test_volume(DIMPLED, 'DIMPLED')
    call test_volume(FOLD_BOX, 'FOLD_BOX')
    call test_volume(TWO_NECKS, 'TWO_NECKS')
    call test_volume(STAR_FAIL, 'STAR_FAIL')
    call test_volume(NEAR_PINCH, 'NEAR_PINCH')
    call test_volume(SEPARATED, 'SEPARATED')
    call test_volume(PINCHED, 'PINCHED')
    call test_volume(PINCHED_ASYM, 'PINCHED_ASYM')
    call test_centre_of_mass(SADDLE2, 'SADDLE2')
    call test_centre_of_mass(ASYM, 'ASYM')
    call test_centre_of_mass(FOLD_BOX, 'FOLD_BOX')
    call test_centre_of_mass(PINCHED_ASYM, 'PINCHED_ASYM')
    call test_symmetric(SPHERE, 'SPHERE')
    call test_symmetric(BARE_05, 'BARE_05')
    call test_symmetric(OBLATE, 'OBLATE')
    call test_symmetric(GS, 'GS')
    call test_symmetric(SADDLE1, 'SADDLE1')
    call test_symmetric(SD_MIN, 'SD_MIN')
    call test_symmetric(DIMPLED, 'DIMPLED')
    call test_symmetric(TWO_NECKS, 'TWO_NECKS')
    call test_symmetric(PINCHED, 'PINCHED')
    call test_symmetric(SEPARATED, 'SEPARATED')
    call test_flags()
    call test_status()
    call test_pinch_quadrature()
    call test_quadrupole()
    call test_summary()

contains

    subroutine resolve_fixture_s(params, res, status)
        real(kind = rk), intent(in) :: params(:)
        type(resolve_t), intent(out) :: res
        integer(kind = ik), intent(out) :: status
        call resolve_s(params, quad_tab, quad_w, check_tab, res, status)
    end subroutine resolve_fixture_s

    !> V11: the cosine grid is symmetric with exact ends and middle; each
    !! half of the t-substituted weights integrates 1 to 1.
    subroutine test_nodes()
        real(kind = rk) :: x(65)
        integer(kind = ik) :: j
        logical :: mirrored
        call cosine_nodes_s(65_ik, x)
        call assert_true(is_zero_f(x(33)), 'V11 cosine middle node is exactly 0')
        call assert_bits_eq(x(1), 1.0_rk, 'V11 x(1) = 1')
        call assert_bits_eq(x(65), -1.0_rk, 'V11 x(65) = -1')
        mirrored = .true.
        do j = 1_ik, 65_ik
            if (abs(x(j) + x(66_ik - j)) > 0.0_rk) mirrored = .false.
        end do
        call assert_true(mirrored, 'V11 x(j) = -x(66 - j) bits')
        ! 2e-15: the exact sum of the weights is 1 + 7.5e-16 and a 64-term double
        ! sum adds ~1e-15; fortran-foundations' own weights miss by 2.4e-14
        call assert_abs_close(sum(quad_w(1:CASSINI_N_GL_HALF)), 1.0_rk, 2.0e-15_rk, 'V11 north weights')
        call assert_abs_close(sum(quad_w(CASSINI_N_GL_HALF + 1:)), 1.0_rk, 2.0e-15_rk, 'V11 south weights')
        call assert_true(all(quad_x(1:CASSINI_N_GL_HALF) > 0.0_rk) .and. &
                all(quad_x(CASSINI_N_GL_HALF + 1:) < 0.0_rk), 'V11 north then south')
    end subroutine test_nodes

    !> V1, V2. The sphere's c and poles hold to 2e-15: a 512-node half summed
    !! in order (Debug) lands 5 ulps from 1.
    subroutine test_sphere_and_bare()
        type(resolve_t) :: res
        integer(kind = ik) :: status
        call resolve_fixture_s(SPHERE, res, status)
        call assert_int_eq(status, 0_ik, 'V1 sphere status')
        call assert_abs_close(res%eps, 0.0_rk, 1.0e-15_rk, 'V1 sphere eps')
        call assert_abs_close(res%c, 1.0_rk, 2.0e-15_rk, 'V1 sphere c')
        call assert_abs_close(res%z_shift, 0.0_rk, 1.0e-15_rk, 'V1 sphere z_shift')
        call assert_abs_close(res%r_north, 1.0_rk, 2.0e-15_rk, 'V1 sphere r_north')
        call assert_abs_close(res%r_south, 1.0_rk, 2.0e-15_rk, 'V1 sphere r_south')
        call resolve_fixture_s(BARE_05, res, status)
        call assert_abs_close(res%c, 0.887236_rk, 1.0e-5_rk, 'V2 BARE_05 c')
        call assert_abs_close(res%r_north, 1.380404_rk, 1.0e-5_rk, 'V2 BARE_05 r_north')
        call assert_abs_close(res%r_north, sqrt(1.5_rk) / res%c, 1.0e-14_rk, 'V2 BARE_05 r_north = sqrt(1.5)/c')
    end subroutine test_sphere_and_bare

    !> V3: the four 236U shapes of J fig. 6.
    subroutine test_236u()
        type(resolve_t) :: res
        integer(kind = ik) :: status
        call resolve_fixture_s(GS, res, status)
        call assert_abs_close(res%c, 0.957813_rk, 1.0e-5_rk, 'V3 GS c')
        call assert_abs_close(res%r_north, 1.227770_rk, 1.0e-5_rk, 'V3 GS r_north')
        call assert_abs_close(res%r_south, 1.227770_rk, 1.0e-5_rk, 'V3 GS r_south')
        call resolve_fixture_s(SADDLE1, res, status)
        call assert_abs_close(res%c, 0.924768_rk, 1.0e-5_rk, 'V3 SADDLE1 c')
        call assert_abs_close(res%r_north, 1.191997_rk, 1.0e-5_rk, 'V3 SADDLE1 r_north')
        call resolve_fixture_s(SD_MIN, res, status)
        call assert_abs_close(res%c, 0.880261_rk, 1.0e-5_rk, 'V3 SD_MIN c')
        call assert_abs_close(res%r_north, 1.427043_rk, 1.0e-5_rk, 'V3 SD_MIN r_north')
        call resolve_fixture_s(SADDLE2, res, status)
        call assert_abs_close(res%c, 0.858229_rk, 1.0e-5_rk, 'V3 SADDLE2 c')
        call assert_abs_close(res%z_shift, 0.316153_rk, 1.0e-5_rk, 'V3 SADDLE2 z_shift')
        call assert_abs_close(res%r_north, 1.443267_rk, 1.0e-5_rk, 'V3 SADDLE2 r_north')
        call assert_abs_close(res%r_south, 1.662467_rk, 1.0e-5_rk, 'V3 SADDLE2 r_south')
        call assert_abs_close(res%eps, 0.654431_rk, 1.0e-5_rk, 'V3 SADDLE2 eps')
    end subroutine test_236u

    !> V4: pi int rho^2 dz on the resolve's own nodes is 4 pi/3 c^3, and the
    !! independent Simpson rule agrees.
    subroutine test_volume(params, name)
        real(kind = rk), intent(in) :: params(:)
        character(len = *), intent(in) :: name
        type(resolve_t) :: res
        integer(kind = ik) :: status
        real(kind = rk) :: r(quad_tab%n_nodes), dr(quad_tab%n_nodes), rho(quad_tab%n_nodes)
        real(kind = rk) :: z(quad_tab%n_nodes), drho(quad_tab%n_nodes), dz(quad_tab%n_nodes)
        real(kind = rk) :: vbar, vol, zcm
        call resolve_fixture_s(params, res, status)
        call assert_int_eq(status, 0_ik, 'V4 status ' // name)
        call series_on_tables_s(params, res%n_active, quad_tab, r, dr)
        call lemniscate_nodes_s(quad_tab%x, r, dr, res%eps, rho, z, drho, dz)
        vbar = PI_C * sum(quad_w * rho**2 * dz)
        call assert_close(vbar / (FOUR_THIRDS_PI * res%c**3), 1.0_rk, 1.0e-13_rk, 'V4 own nodes ' // name)
        call reference_volume_s(params, res, vol, zcm)
        call assert_close(vol / (FOUR_THIRDS_PI * res%c**3), 1.0_rk, 1.0e-7_rk, 'V4 Simpson ' // name)
    end subroutine test_volume

    !> V5: the raw centre of mass against the Simpson rule.
    subroutine test_centre_of_mass(params, name)
        real(kind = rk), intent(in) :: params(:)
        character(len = *), intent(in) :: name
        type(resolve_t) :: res
        integer(kind = ik) :: status
        real(kind = rk) :: vol, zcm
        call resolve_fixture_s(params, res, status)
        call reference_volume_s(params, res, vol, zcm)
        call assert_abs_close(res%zm_raw, zcm, 1.0e-7_rk, 'V5 zm_raw vs Simpson ' // name)
        call assert_abs_close(res%z_shift * res%c, res%zm_raw, 1.0e-15_rk, 'V5 z_shift = zm_raw/c ' // name)
    end subroutine test_centre_of_mass

    !> V6: reflection-symmetric vectors have z_shift = 0 and equal poles.
    subroutine test_symmetric(params, name)
        real(kind = rk), intent(in) :: params(:)
        character(len = *), intent(in) :: name
        type(resolve_t) :: res
        integer(kind = ik) :: status
        call resolve_fixture_s(params, res, status)
        call assert_true(is_zero_f(res%z_shift) .or. abs(res%z_shift) <= 1.0e-16_rk, 'V6 z_shift = 0 ' // name)
        call assert_bits_eq(res%r_north, res%r_south, 'V6 r_north = r_south ' // name)
    end subroutine test_symmetric

    !> V7: pinched and separated flags, the gap.
    subroutine test_flags()
        type(resolve_t) :: res
        integer(kind = ik) :: status
        call resolve_fixture_s(PINCHED, res, status)
        call assert_true(res%is_pinched .and. .not. res%is_separated, 'V7 PINCHED is_pinched')
        call assert_true(is_zero_f(res%gap), 'V7 PINCHED gap = 0')
        call resolve_fixture_s(PINCHED_ASYM, res, status)
        call assert_true(res%is_pinched .and. .not. res%is_separated, 'V7 PINCHED_ASYM is_pinched')
        call assert_true(is_zero_f(res%gap), 'V7 PINCHED_ASYM gap = 0')
        call resolve_fixture_s(SEPARATED, res, status)
        call assert_true(res%is_separated .and. .not. res%is_pinched, 'V7 SEPARATED is_separated')
        call assert_abs_close(res%gap, 2.0_rk * sqrt(res%eps - 1.0_rk) / res%c, 1.0e-14_rk, 'V7 SEPARATED gap')
        call assert_abs_close(res%gap, 0.687685_rk, 1.0e-5_rk, 'V7 SEPARATED gap value')
        call resolve_fixture_s(NEAR_PINCH, res, status)
        call assert_true(.not. res%is_pinched .and. .not. res%is_separated, 'V7 NEAR_PINCH connected')
        call assert_true(is_zero_f(res%gap), 'V7 NEAR_PINCH gap = 0')
    end subroutine test_flags

    !> V8: 100 before 101; failures zero-fill the record.
    subroutine test_status()
        type(resolve_t) :: res
        integer(kind = ik) :: status
        call resolve_fixture_s(R_NEG, res, status)
        call assert_int_eq(status, CASSINI_ERROR_R_NONPOSITIVE, 'V8 R_NEG -> 100 (its pole is also toroidal)')
        call assert_true(is_zero_f(res%c) .and. is_zero_f(res%eps) .and. is_zero_f(res%r_north) &
                .and. is_zero_f(res%r_south) .and. is_zero_f(res%zm_raw) .and. res%n_active == 0_ik &
                .and. .not. res%is_pinched .and. .not. res%is_separated, 'V8 R_NEG zero-filled')
        call resolve_fixture_s(TOROIDAL, res, status)
        call assert_int_eq(status, CASSINI_ERROR_TOROIDAL, 'V8 TOROIDAL -> 101')
        call assert_true(is_zero_f(res%c), 'V8 TOROIDAL zero-filled')
        call resolve_fixture_s(TOROIDAL_ONE, res, status)
        call assert_int_eq(status, CASSINI_ERROR_TOROIDAL, 'V8 TOROIDAL_ONE -> 101')
    end subroutine test_status

    !> V9: 64 nodes in t hold 1e-13 through the pinch (against 512 in t).
    subroutine test_pinch_quadrature()
        real(kind = rk) :: walk(1)
        integer(kind = ik) :: k
        real(kind = rk), parameter :: ALPHAS(3) = [1.0_rk - 1.0e-9_rk, 1.0_rk, 1.0_rk + 1.0e-9_rk]
        type(resolve_t) :: res
        integer(kind = ik) :: status
        real(kind = rk) :: vbar, zm
        do k = 1_ik, 3_ik
            walk(1) = ALPHAS(k)
            call resolve_fixture_s(walk, res, status)
            call reference_gl_volume_s(walk, res%n_active, res%eps, 512_ik, .true., vbar, zm)
            call assert_close(FOUR_THIRDS_PI * res%c**3 / vbar, 1.0_rk, 1.0e-13_rk, 'V9 pinch walk volume')
        end do
    end subroutine test_pinch_quadrature

    !> V10: q2 values and the final-frame Simpson rule.
    subroutine test_quadrupole()
        type(resolve_t) :: res
        integer(kind = ik) :: status
        call resolve_fixture_s(SPHERE, res, status)
        call assert_abs_close(quadrupole_f(SPHERE, res, quad_tab, quad_w), 0.0_rk, 1.0e-15_rk, 'V10 sphere q2')
        call resolve_fixture_s(BARE_05, res, status)
        call assert_abs_close(quadrupole_f(BARE_05, res, quad_tab, quad_w), 0.288240_rk, 1.0e-5_rk, 'V10 BARE_05 q2')
        call resolve_fixture_s(SADDLE2, res, status)
        call assert_abs_close(quadrupole_f(SADDLE2, res, quad_tab, quad_w), 0.394688_rk, 1.0e-5_rk, 'V10 SADDLE2 q2')
        call assert_abs_close(quadrupole_f(SADDLE2, res, quad_tab, quad_w), &
                reference_quadrupole_f(SADDLE2, res), 1.0e-7_rk, 'V10 SADDLE2 q2 vs Simpson')
        call resolve_fixture_s(ASYM, res, status)
        call assert_abs_close(quadrupole_f(ASYM, res, quad_tab, quad_w), 0.687640_rk, 1.0e-5_rk, 'V10 ASYM q2')
        call resolve_fixture_s(SEPARATED, res, status)
        call assert_abs_close(quadrupole_f(SEPARATED, res, quad_tab, quad_w), 1.824705_rk, 1.0e-5_rk, &
                'V10 SEPARATED q2')
        call resolve_fixture_s(FOLD_BOX, res, status)
        call assert_abs_close(quadrupole_f(FOLD_BOX, res, quad_tab, quad_w), 0.083071_rk, 1.0e-5_rk, &
                'V10 FOLD_BOX q2')
        call assert_abs_close(quadrupole_f(FOLD_BOX, res, quad_tab, quad_w), &
                reference_quadrupole_f(FOLD_BOX, res), 1.0e-7_rk, 'V10 FOLD_BOX q2 vs Simpson')
    end subroutine test_quadrupole

end program cassini_param_resolve_test
