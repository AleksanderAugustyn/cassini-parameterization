!> Contract family 1, equivalence: one-shot == cached bitwise on every output
!! and status branch, caches of different max_params agree, a short vector
!! equals its zero-padded form, and rejections agree across tiers.
program cassini_param_equivalence_test
    use precision_utilities_mod, only: ik, rk
    use mathematical_and_physical_constants_mod, only: PI_C
    use cassini_parameterization_mod, only: cache_t, cache_radius_grid_s, &
            cache_radius_and_derivative_at_thetas_s, compute_radius_grid_standalone_s, &
            SHAPE_ERROR_TOO_MANY_PARAMS, SHAPE_ERROR_WRONG_PARAM_COUNT, SHAPE_ERROR_INVALID_GRID, &
            CASSINI_ERROR_R_NONPOSITIVE, CASSINI_ERROR_TOROIDAL, CASSINI_ERROR_NOT_STAR_CONVEX, &
            CASSINI_ERROR_SEPARATED
    use cassini_test_fixtures_mod, only: SPHERE, BARE_05, GS, SADDLE2, ASYM, DIMPLED, FOLD_BOX, &
            TWO_NECKS, SEPARATED, PINCHED, STAR_FAIL, OBLATE, R_NEG, TOROIDAL
    use cassini_test_tiers_mod, only: result_t, run_cached, run_one_shot, build_cache, &
            results_bits_eq, output_zero_filled, N_OUT, OUT_RADIUS_GRID
    use test_utils_mod, only: assert_true, assert_int_eq, all_zero_f, test_summary
    implicit none

    integer(kind = ik), parameter :: N_X = 129_ik, N_Z = 257_ik
    real(kind = rk) :: thetas(32)
    type(cache_t) :: cache5, cache21
    integer(kind = ik) :: k, status

    do k = 1_ik, 32_ik
        thetas(k) = PI_C * (real(k, rk) - 0.5_rk) / 32.0_rk
    end do
    call build_cache(5_ik, N_X, thetas, cache5, status)
    call build_cache(21_ik, N_X, thetas, cache21, status)

    call test_tiers(SPHERE, 'SPHERE')
    call test_tiers(BARE_05, 'BARE_05')
    call test_tiers(GS, 'GS')
    call test_tiers(SADDLE2, 'SADDLE2')
    call test_tiers(ASYM, 'ASYM')
    call test_tiers(DIMPLED, 'DIMPLED')
    call test_tiers(FOLD_BOX, 'FOLD_BOX')
    call test_tiers(TWO_NECKS, 'TWO_NECKS')
    call test_tiers(SEPARATED, 'SEPARATED')
    call test_tiers(PINCHED, 'PINCHED')
    call test_tiers(STAR_FAIL, 'STAR_FAIL')
    call test_tiers(OBLATE, 'OBLATE')
    call test_cache_sizes(GS, 'GS')
    call test_cache_sizes(SADDLE2, 'SADDLE2')
    call test_long_vectors()
    call test_short_padded()
    call test_rejections(R_NEG, CASSINI_ERROR_R_NONPOSITIVE, OUT_RADIUS_GRID, 'R_NEG')
    call test_rejections(TOROIDAL, CASSINI_ERROR_TOROIDAL, OUT_RADIUS_GRID, 'TOROIDAL')
    call test_rejections(STAR_FAIL, CASSINI_ERROR_NOT_STAR_CONVEX, OUT_RADIUS_GRID, 'STAR_FAIL')
    call test_rejections(SEPARATED, CASSINI_ERROR_SEPARATED, OUT_RADIUS_GRID, 'SEPARATED')
    call test_bad_theta()
    call test_summary()

contains

    !> Q1: every output and status, one-shot against cached.
    subroutine test_tiers(params, name)
        real(kind = rk), intent(in) :: params(:)
        character(len = *), intent(in) :: name
        type(result_t) :: cached, one_shot
        call run_cached(cache5, params, thetas, N_Z, cached)
        call run_one_shot(params, N_X, thetas, N_Z, one_shot)
        call assert_true(results_bits_eq(cached, one_shot), 'Q1 one-shot == cached ' // name)
    end subroutine test_tiers

    !> Q2: a max_params = 21 cache and a max_params = 5 cache agree.
    subroutine test_cache_sizes(params, name)
        real(kind = rk), intent(in) :: params(:)
        character(len = *), intent(in) :: name
        type(result_t) :: a, b
        call run_cached(cache5, params, thetas, N_Z, a)
        call run_cached(cache21, params, thetas, N_Z, b)
        call assert_true(results_bits_eq(a, b), 'Q2 21-cache == 5-cache ' // name)
    end subroutine test_cache_sizes

    !> Q2: a 15-vector, its 21-padding and GS agree; 22 is never trimmed first.
    subroutine test_long_vectors()
        type(result_t) :: a, b, c
        real(kind = rk) :: v15(15), v21(21), v22(22), radii(32)
        v15 = 0.0_rk
        v15(1:5) = GS
        v21 = 0.0_rk
        v21(1:5) = GS
        v22 = 0.0_rk
        v22(1:5) = GS
        call run_cached(cache21, v15, thetas, N_Z, a)
        call run_cached(cache21, v21, thetas, N_Z, b)
        call run_cached(cache5, GS, thetas, N_Z, c)
        call assert_true(results_bits_eq(a, b), 'Q2 15-vector == 21-padded on the 21-cache')
        call assert_true(results_bits_eq(a, c), 'Q2 15-vector == GS on the 5-cache')
        call cache_radius_grid_s(cache21, v22, radii, status)
        call assert_int_eq(status, SHAPE_ERROR_WRONG_PARAM_COUNT, 'Q2 22-vector on the 21-cache -> 4')
        call compute_radius_grid_standalone_s(v22, N_X, thetas, radii, status)
        call assert_int_eq(status, SHAPE_ERROR_TOO_MANY_PARAMS, 'Q2 22-vector one-shot -> 1')
    end subroutine test_long_vectors

    !> Q3: a short vector and its zero-padded form, both tiers.
    subroutine test_short_padded()
        type(result_t) :: a, b
        real(kind = rk), parameter :: HALF(1) = [0.5_rk]
        real(kind = rk), parameter :: HALF_PADDED(5) = [0.5_rk, 0.0_rk, 0.0_rk, 0.0_rk, 0.0_rk]
        real(kind = rk) :: saddle_padded(8)
        type(cache_t) :: cache8
        saddle_padded = 0.0_rk
        saddle_padded(1:5) = SADDLE2
        call run_cached(cache5, HALF, thetas, N_Z, a)
        call run_cached(cache5, HALF_PADDED, thetas, N_Z, b)
        call assert_true(results_bits_eq(a, b), 'Q3 [0.5] == padded (cached)')
        call run_one_shot(HALF, N_X, thetas, N_Z, a)
        call run_one_shot(HALF_PADDED, N_X, thetas, N_Z, b)
        call assert_true(results_bits_eq(a, b), 'Q3 [0.5] == padded (one-shot)')
        call build_cache(8_ik, N_X, thetas, cache8, status)
        call run_cached(cache8, SADDLE2, thetas, N_Z, a)
        call run_cached(cache8, saddle_padded, thetas, N_Z, b)
        call assert_true(results_bits_eq(a, b), 'Q3 SADDLE2 == SADDLE2 // [0, 0, 0] (cached)')
        call run_one_shot(SADDLE2, N_X, thetas, N_Z, a)
        call run_one_shot(saddle_padded, N_X, thetas, N_Z, b)
        call assert_true(results_bits_eq(a, b), 'Q3 SADDLE2 == SADDLE2 // [0, 0, 0] (one-shot)')
    end subroutine test_short_padded

    !> Q4: a rejection has the same code and zero-filled outputs in both tiers.
    subroutine test_rejections(params, code, slot, name)
        real(kind = rk), intent(in) :: params(:)
        integer(kind = ik), intent(in) :: code, slot
        character(len = *), intent(in) :: name
        type(result_t) :: cached, one_shot
        integer(kind = ik) :: j
        logical :: zeroed
        call run_cached(cache5, params, thetas, N_Z, cached)
        call run_one_shot(params, N_X, thetas, N_Z, one_shot)
        call assert_int_eq(cached%status(slot), code, 'Q4 cached code ' // name)
        call assert_int_eq(one_shot%status(slot), code, 'Q4 one-shot code ' // name)
        zeroed = .true.
        do j = 1_ik, N_OUT
            if (cached%status(j) /= 0_ik) zeroed = zeroed .and. output_zero_filled(cached, j) &
                    .and. output_zero_filled(one_shot, j)
        end do
        call assert_true(zeroed .and. results_bits_eq(cached, one_shot), 'Q4 zero-filled both sides ' // name)
    end subroutine test_rejections

    !> Q4: a theta of 3 pi/2 is 3 in both tiers.
    subroutine test_bad_theta()
        real(kind = rk) :: bad(2), radii(2), dr(2)
        bad = [0.5_rk, 1.5_rk * PI_C]
        call cache_radius_and_derivative_at_thetas_s(cache5, GS, bad, radii, dr, status)
        call assert_int_eq(status, SHAPE_ERROR_INVALID_GRID, 'Q4 bad theta cached -> 3')
        call assert_true(all_zero_f(radii) .and. all_zero_f(dr), 'Q4 bad theta cached zero-filled')
        call compute_radius_grid_standalone_s(GS, N_X, bad, radii, status)
        call assert_int_eq(status, SHAPE_ERROR_INVALID_GRID, 'Q4 bad theta one-shot -> 3')
        call assert_true(all_zero_f(radii), 'Q4 bad theta one-shot zero-filled')
    end subroutine test_bad_theta

end program cassini_param_equivalence_test
