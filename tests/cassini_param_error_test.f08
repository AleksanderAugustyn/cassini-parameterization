!> The status matrix on every output of both tiers: valid rows, 100, 101,
!! 103, 102 (pinch and margin), folded shapes, precedence, and the bitwise
!! relations between outputs (at-thetas, unchecked, trimming).
program cassini_param_error_test
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use precision_utilities_mod, only: ik, rk
    use mathematical_and_physical_constants_mod, only: PI_C
    use cassini_parameterization_mod, only: cache_t, cache_init_s, &
            cache_radius_and_derivative_at_thetas_s, cache_radius_grid_s, cache_rho_profile_s, &
            compute_rho_profile_standalone_s, compute_radius_grid_standalone_s, compute_epsilon_f, &
            CASSINI_ERROR_R_NONPOSITIVE, CASSINI_ERROR_TOROIDAL, &
            CASSINI_ERROR_NOT_STAR_CONVEX, CASSINI_ERROR_SEPARATED, SHAPE_ERROR_TOO_MANY_PARAMS, &
            SHAPE_ERROR_CACHE_NOT_INITIALIZED, SHAPE_ERROR_INVALID_GRID, SHAPE_ERROR_WRONG_PARAM_COUNT
    use cassini_test_fixtures_mod, only: SPHERE, GS, SADDLE2, ASYM, DIMPLED, FOLD_BOX, STAR_FAIL, &
            PINCHED, SEPARATED, R_NEG, TOROIDAL, TOROIDAL_ONE
    use cassini_test_tiers_mod, only: result_t, run_cached, run_one_shot, build_cache, &
            output_zero_filled, results_bits_eq, N_OUT
    use test_utils_mod, only: assert_true, assert_int_eq, assert_bits_eq, &
            arrays_bits_eq_f, is_zero_f, all_zero_f, test_summary
    implicit none

    integer(kind = ik), parameter :: N_X = 129_ik, N_Z = 101_ik
    integer(kind = ik), parameter :: R_PATH = CASSINI_ERROR_NOT_STAR_CONVEX
    real(kind = rk) :: thetas(8)
    type(cache_t) :: cache
    integer(kind = ik) :: k, status

    do k = 1_ik, 8_ik
        thetas(k) = PI_C * (real(k, rk) - 0.5_rk) / 8.0_rk
    end do
    call build_cache(5_ik, N_X, thetas, cache, status)
    call assert_int_eq(status, 0_ik, 'cache init')

    call test_valid_row()
    call test_row(R_NEG, [(CASSINI_ERROR_R_NONPOSITIVE, k = 1, 9)], 'E2 R_NEG')
    call test_row(TOROIDAL, [(CASSINI_ERROR_TOROIDAL, k = 1, 9)], 'E2 TOROIDAL')
    call test_row(TOROIDAL_ONE, [(CASSINI_ERROR_TOROIDAL, k = 1, 9)], 'E2 TOROIDAL_ONE')
    call test_separated()
    call test_pinched()
    call test_star_fail()
    call test_folded(DIMPLED, 'DIMPLED')
    call test_folded(FOLD_BOX, 'FOLD_BOX')
    call test_precedence()
    call test_at_thetas_bits(SADDLE2, 'SADDLE2')
    call test_at_thetas_bits(DIMPLED, 'DIMPLED')
    call test_unchecked_bits(ASYM, 'ASYM')
    call test_unchecked_bits(FOLD_BOX, 'FOLD_BOX')
    call test_trimming()
    call test_summary()

contains

    !> Both tiers of one parameter vector; every status as expected, every
    !! failing output zero-filled, and the two tiers bitwise equal.
    subroutine test_row(params, expected, name)
        real(kind = rk), intent(in) :: params(:)
        integer(kind = ik), intent(in) :: expected(N_OUT)
        character(len = *), intent(in) :: name
        type(result_t) :: cached, one_shot
        integer(kind = ik) :: slot
        logical :: zero_ok
        call run_cached(cache, params, thetas, N_Z, cached)
        call run_one_shot(params, N_X, thetas, N_Z, one_shot)
        do slot = 1_ik, N_OUT
            call assert_int_eq(cached%status(slot), expected(slot), name // ' cached status')
            call assert_int_eq(one_shot%status(slot), expected(slot), name // ' one-shot status')
        end do
        zero_ok = .true.
        do slot = 1_ik, N_OUT
            if (expected(slot) /= 0_ik) then
                if (.not. (output_zero_filled(cached, slot) .and. output_zero_filled(one_shot, slot))) &
                        zero_ok = .false.
            end if
        end do
        call assert_true(zero_ok, name // ' failing outputs zero-filled')
        call assert_true(results_bits_eq(cached, one_shot), name // ' one-shot == cached bits')
    end subroutine test_row

    !> E1: GS is valid everywhere; the poles agree across outputs.
    subroutine test_valid_row()
        type(result_t) :: cached
        real(kind = rk) :: poles(2), radii(2), dr(2)
        poles = [0.0_rk, PI_C]
        call test_row(GS, [(0_ik, k = 1, 9)], 'E1 GS')
        call run_cached(cache, GS, thetas, N_Z, cached)
        call cache_radius_and_derivative_at_thetas_s(cache, GS, poles, radii, dr, status)
        call assert_bits_eq(cached%r_north, cached%geom%r_north, 'E1 shape r_north = geometry')
        call assert_bits_eq(cached%r_south, cached%geom%r_south, 'E1 shape r_south = geometry')
        call assert_bits_eq(radii(1), cached%r_north, 'E1 R(0) = r_north')
        call assert_bits_eq(radii(2), cached%r_south, 'E1 R(pi) = r_south')
    end subroutine test_valid_row

    !> E3: separated shapes are 103 on the R(theta) paths, shape and the
    !! checked profile; everything else draws them.
    subroutine test_separated()
        type(result_t) :: cached
        integer(kind = ik), parameter :: SEP = CASSINI_ERROR_SEPARATED
        call test_row(SEPARATED, [SEP, SEP, SEP, SEP, SEP, 0_ik, 0_ik, 0_ik, 0_ik], 'E3 SEPARATED')
        call run_cached(cache, SEPARATED, thetas, N_Z, cached)
        call assert_true(is_zero_f(cached%rhou(51)) .and. is_zero_f(cached%rhou(52)), &
                'E3 unchecked profile: tips at rho = 0')
        call assert_true(cached%rhou(50) > 0.0_rk .and. cached%rhou(53) > 0.0_rk, &
                'E3 unchecked profile: rho > 0 beside the tips')
        call assert_true(cached%found .and. is_zero_f(cached%rho_neck), 'E3 neck found, rho_neck = 0')
        call assert_true(cached%geom%is_separated .and. cached%geom%r_north > 0.0_rk &
                .and. cached%geom%r_south > 0.0_rk .and. cached%geom%gap > 0.0_rk &
                .and. is_zero_f(cached%geom%d_min), 'E3 geometry of a separated shape')
    end subroutine test_separated

    !> E4: the pinched shape is 102 on the R(theta) paths by rule.
    subroutine test_pinched()
        type(result_t) :: cached
        call test_row(PINCHED, [R_PATH, R_PATH, R_PATH, R_PATH, 0_ik, 0_ik, 0_ik, 0_ik, 0_ik], 'E4 PINCHED')
        call run_cached(cache, PINCHED, thetas, N_Z, cached)
        call assert_true(arrays_bits_eq_f(cached%z, cached%zu) .and. arrays_bits_eq_f(cached%rho, cached%rhou) &
                .and. arrays_bits_eq_f(cached%drho, cached%drhou), 'E4 checked == unchecked profile bits')
        call assert_true(is_zero_f(cached%rho(51)) .and. is_zero_f(cached%rho(52)), 'E4 pinch node rho = 0')
        call assert_true(cached%found .and. is_zero_f(cached%rho_neck), 'E4 neck found, rho_neck = 0')
        call assert_true(cached%geom%is_pinched .and. .not. cached%geom%is_separated, 'E4 geometry is_pinched')
    end subroutine test_pinched

    !> E5: STAR_FAIL is 102 on the R(theta) paths by the margin.
    subroutine test_star_fail()
        type(result_t) :: cached
        call test_row(STAR_FAIL, [R_PATH, R_PATH, R_PATH, R_PATH, 0_ik, 0_ik, 0_ik, 0_ik, 0_ik], 'E5 STAR_FAIL')
        call run_cached(cache, STAR_FAIL, thetas, N_Z, cached)
        call assert_true(cached%found .and. cached%rho_neck > 0.0_rk, 'E5 neck found')
        call assert_true(cached%geom%d_min < 0.0_rk .and. .not. cached%geom%is_folded, 'E5 d_min < 0, not folded')
    end subroutine test_star_fail

    !> E6: folded shapes are valid everywhere.
    subroutine test_folded(params, name)
        real(kind = rk), intent(in) :: params(:)
        character(len = *), intent(in) :: name
        type(result_t) :: cached
        call test_row(params, [(0_ik, k = 1, 9)], 'E6 ' // name)
        call run_cached(cache, params, thetas, N_Z, cached)
        call assert_true(cached%geom%is_folded, 'E6 is_folded ' // name)
        call assert_true(all(ieee_is_finite(cached%radii_d)) .and. all(ieee_is_finite(cached%dr_d)), &
                'E6 radius_and_derivative finite ' // name)
    end subroutine test_folded

    !> E7: usage codes come before value codes, in their own order.
    subroutine test_precedence()
        type(cache_t) :: cache21, blank
        real(kind = rk) :: long(22), empty(0), radii(8), bad(2), z(N_Z), rho(N_Z), drho(N_Z - 1_ik)
        real(kind = rk) :: radii2(2), dr2(2)
        long = 0.0_rk
        long(1) = SEPARATED(1)
        call cache_init_s(cache21, 21_ik, N_X, thetas, status)
        call compute_radius_grid_standalone_s(long, N_X, thetas, radii, status)
        call assert_int_eq(status, SHAPE_ERROR_TOO_MANY_PARAMS, 'E7 one-shot 22 params -> 1')
        call cache_radius_grid_s(cache21, long, radii, status)
        call assert_int_eq(status, SHAPE_ERROR_WRONG_PARAM_COUNT, 'E7 cached 22 params -> 4')
        call compute_radius_grid_standalone_s(empty, N_X, thetas, radii, status)
        call assert_int_eq(status, SHAPE_ERROR_WRONG_PARAM_COUNT, 'E7 one-shot empty -> 4')
        call cache_radius_grid_s(cache21, empty, radii, status)
        call assert_int_eq(status, SHAPE_ERROR_WRONG_PARAM_COUNT, 'E7 cached empty -> 4')
        bad = [0.5_rk, 1.5_rk * PI_C]
        call cache_radius_and_derivative_at_thetas_s(cache, PINCHED, bad, radii2, dr2, status)
        call assert_int_eq(status, SHAPE_ERROR_INVALID_GRID, 'E7 PINCHED at_thetas bad theta -> 3')
        call assert_true(all_zero_f(radii2) .and. all_zero_f(dr2), 'E7 bad theta zero-filled')
        call cache_rho_profile_s(cache, SEPARATED, z, rho, drho, status)
        call assert_int_eq(status, SHAPE_ERROR_INVALID_GRID, 'E7 SEPARATED profile extents -> 3 (cached)')
        call compute_rho_profile_standalone_s(SEPARATED, z, rho, drho, status)
        call assert_int_eq(status, SHAPE_ERROR_INVALID_GRID, 'E7 SEPARATED profile extents -> 3 (one-shot)')
        call cache_radius_grid_s(blank, long, radii, status)
        call assert_int_eq(status, SHAPE_ERROR_CACHE_NOT_INITIALIZED, 'E7 uninitialized -> 2 before 4')
        call cache_radius_and_derivative_at_thetas_s(cache, R_NEG, bad, radii2, dr2, status)
        call assert_int_eq(status, SHAPE_ERROR_INVALID_GRID, 'E7 R_NEG with a bad theta -> 3')
    end subroutine test_precedence

    !> E8: at_thetas at the cache's own thetas is radius_and_derivative.
    subroutine test_at_thetas_bits(params, name)
        real(kind = rk), intent(in) :: params(:)
        character(len = *), intent(in) :: name
        type(result_t) :: cached
        call run_cached(cache, params, thetas, N_Z, cached)
        call assert_true(arrays_bits_eq_f(cached%radii_a, cached%radii_d) &
                .and. arrays_bits_eq_f(cached%dr_a, cached%dr_d) &
                .and. arrays_bits_eq_f(cached%radii_g, cached%radii_d), 'E8 at_thetas == grid bits ' // name)
    end subroutine test_at_thetas_bits

    !> E9: unchecked == checked on connected shapes; the sphere is 1.
    subroutine test_unchecked_bits(params, name)
        real(kind = rk), intent(in) :: params(:)
        character(len = *), intent(in) :: name
        type(result_t) :: cached
        call run_cached(cache, params, thetas, N_Z, cached)
        call assert_true(arrays_bits_eq_f(cached%z, cached%zu) .and. arrays_bits_eq_f(cached%rho, cached%rhou) &
                .and. arrays_bits_eq_f(cached%drho, cached%drhou), 'E9 unchecked == checked ' // name)
        call run_cached(cache, SPHERE, thetas, N_Z, cached)
        call assert_true(maxval(abs(cached%radii_g - 1.0_rk)) <= 1.0e-14_rk, 'E9 sphere radius_grid = 1')
    end subroutine test_unchecked_bits

    !> E10: trailing zeros change nothing; compute_epsilon_f.
    subroutine test_trimming()
        type(result_t) :: a, b
        type(cache_t) :: cache7
        real(kind = rk), parameter :: GS_PADDED(7) = [0.25_rk, 0.0_rk, 0.0_rk, 0.0_rk, 0.075_rk, 0.0_rk, 0.0_rk]
        real(kind = rk), parameter :: HALF(1) = [0.5_rk], HALF_PADDED(3) = [0.5_rk, 0.0_rk, 0.0_rk]
        real(kind = rk) :: empty(0)
        call cache_init_s(cache7, 7_ik, N_X, thetas, status)
        call run_cached(cache7, GS, thetas, N_Z, a)
        call run_cached(cache7, GS_PADDED, thetas, N_Z, b)
        call assert_true(results_bits_eq(a, b), 'E10 GS padded == GS (cached)')
        call run_one_shot(GS, N_X, thetas, N_Z, a)
        call run_one_shot(GS_PADDED, N_X, thetas, N_Z, b)
        call assert_true(results_bits_eq(a, b), 'E10 GS padded == GS (one-shot)')
        call run_one_shot(HALF, N_X, thetas, N_Z, a)
        call run_one_shot(HALF_PADDED, N_X, thetas, N_Z, b)
        call assert_true(results_bits_eq(a, b), 'E10 [0.5] padded == [0.5]')
        call assert_bits_eq(compute_epsilon_f(empty), 0.0_rk, 'E10 epsilon([]) = 0')
        call assert_bits_eq(compute_epsilon_f(HALF), 0.5_rk, 'E10 epsilon([0.5]) = 0.5')
    end subroutine test_trimming

end program cassini_param_error_test
