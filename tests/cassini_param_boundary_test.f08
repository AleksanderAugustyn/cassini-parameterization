!> Contract boundaries: init arguments, parameter counts in both tiers, the
!! last parameter is live, the profile floor, buffer extents, theta counts.
program cassini_param_boundary_test
    use precision_utilities_mod, only: ik, rk
    use mathematical_and_physical_constants_mod, only: PI_C
    use cassini_parameterization_mod, only: cache_t, cache_init_s, cache_radius_grid_s, &
            cache_radius_and_derivative_s, cache_radius_and_derivative_at_thetas_s, cache_shape_s, &
            cache_rho_profile_s, cache_rho_profile_unchecked_s, cache_geometry_s, cassini_geometry_t, &
            compute_radius_grid_standalone_s, compute_rho_profile_standalone_s, &
            compute_rho_profile_unchecked_standalone_s, cache_max_params_f, cache_n_x_grid_f, &
            cache_n_thetas_f, cache_is_initialized_f, cache_free_s, SHAPE_ERROR_TOO_MANY_PARAMS, &
            SHAPE_ERROR_INVALID_GRID, SHAPE_ERROR_WRONG_PARAM_COUNT, SHAPE_ERROR_INVALID_INIT
    use cassini_test_fixtures_mod, only: GS, SEPARATED
    use test_utils_mod, only: assert_true, assert_int_eq, arrays_bits_eq_f, is_zero_f, test_summary
    implicit none

    real(kind = rk) :: thetas(8)
    integer(kind = ik) :: k

    do k = 1_ik, 8_ik
        thetas(k) = PI_C * (real(k, rk) - 0.5_rk) / 8.0_rk
    end do

    call test_init()
    call test_cached_counts()
    call test_one_shot_counts()
    call test_last_parameter_live()
    call test_profile_floor()
    call test_extents()
    call test_theta_counts()
    call test_small_cache()
    call test_summary()

contains

    !> B1: init arguments.
    subroutine test_init()
        type(cache_t) :: cache
        integer(kind = ik) :: status
        real(kind = rk) :: empty(0), one(1)
        call cache_init_s(cache, 21_ik, 65_ik, thetas, status)
        call assert_int_eq(status, 0_ik, 'B1 max_params 21 -> 0')
        call assert_true(cache_is_initialized_f(cache) .and. cache_max_params_f(cache) == 21_ik &
                .and. cache_n_x_grid_f(cache) == 65_ik .and. cache_n_thetas_f(cache) == 8_ik, 'B1 accessors')
        call cache_free_s(cache)
        call assert_true(.not. cache_is_initialized_f(cache) .and. cache_n_thetas_f(cache) == 0_ik, 'B1 free')
        call cache_init_s(cache, 22_ik, 65_ik, thetas, status)
        call assert_int_eq(status, SHAPE_ERROR_TOO_MANY_PARAMS, 'B1 max_params 22 -> 1')
        call assert_true(.not. cache_is_initialized_f(cache), 'B1 rejected cache stays uninitialized')
        call cache_init_s(cache, 1_ik, 65_ik, thetas, status)
        call assert_int_eq(status, 0_ik, 'B1 max_params 1 -> 0')
        call cache_init_s(cache, 0_ik, 65_ik, thetas, status)
        call assert_int_eq(status, SHAPE_ERROR_INVALID_INIT, 'B1 max_params 0 -> 5')
        call cache_init_s(cache, 5_ik, 64_ik, thetas, status)
        call assert_int_eq(status, SHAPE_ERROR_INVALID_GRID, 'B1 n_x_grid 64 -> 3')
        call cache_init_s(cache, 5_ik, 2049_ik, thetas, status)
        call assert_int_eq(status, 0_ik, 'B1 n_x_grid 2049 -> 0')
        call cache_init_s(cache, 5_ik, 65_ik, empty, status)
        call assert_int_eq(status, SHAPE_ERROR_INVALID_GRID, 'B1 empty thetas -> 3')
        one(1) = 1.5_rk * PI_C
        call cache_init_s(cache, 5_ik, 65_ik, one, status)
        call assert_int_eq(status, SHAPE_ERROR_INVALID_GRID, 'B1 theta 3 pi/2 -> 3')
        one(1) = PI_C
        call cache_init_s(cache, 5_ik, 65_ik, one, status)
        call assert_int_eq(status, 0_ik, 'B1 theta pi -> 0')
        one(1) = 0.0_rk
        call cache_init_s(cache, 5_ik, 65_ik, one, status)
        call assert_int_eq(status, 0_ik, 'B1 theta 0 -> 0')
    end subroutine test_init

    !> B2: cached size(params) on a max_params = 5 cache.
    subroutine test_cached_counts()
        type(cache_t) :: cache
        type(cassini_geometry_t) :: geometry
        integer(kind = ik) :: status, n, s_grid, s_shape, s_profile, s_geom, want
        real(kind = rk) :: params(6), radii(8), z_shift, r_north, r_south, z(65), rho(65), drho(65)
        call cache_init_s(cache, 5_ik, 65_ik, thetas, status)
        params = 0.0_rk
        params(1) = 0.25_rk
        do n = 0_ik, 6_ik
            call cache_radius_grid_s(cache, params(1:n), radii, s_grid)
            call cache_shape_s(cache, params(1:n), z_shift, r_north, r_south, s_shape)
            call cache_rho_profile_s(cache, params(1:n), z, rho, drho, s_profile)
            call cache_geometry_s(cache, params(1:n), geometry, s_geom)
            want = 0_ik
            if (n == 0_ik .or. n == 6_ik) want = SHAPE_ERROR_WRONG_PARAM_COUNT
            call assert_true(s_grid == want .and. s_shape == want .and. s_profile == want &
                    .and. s_geom == want, 'B2 cached count')
        end do
    end subroutine test_cached_counts

    !> B3: one-shot size(params).
    subroutine test_one_shot_counts()
        real(kind = rk) :: params(22), radii(8)
        integer(kind = ik) :: status
        params = 0.0_rk
        params(1) = 0.25_rk
        call compute_radius_grid_standalone_s(params(1:1), 65_ik, thetas, radii, status)
        call assert_int_eq(status, 0_ik, 'B3 one-shot 1 -> 0')
        call compute_radius_grid_standalone_s(params(1:21), 65_ik, thetas, radii, status)
        call assert_int_eq(status, 0_ik, 'B3 one-shot 21 -> 0')
        call compute_radius_grid_standalone_s(params, 65_ik, thetas, radii, status)
        call assert_int_eq(status, SHAPE_ERROR_TOO_MANY_PARAMS, 'B3 one-shot 22 -> 1')
        call compute_radius_grid_standalone_s(params(1:0), 65_ik, thetas, radii, status)
        call assert_int_eq(status, SHAPE_ERROR_WRONG_PARAM_COUNT, 'B3 one-shot empty -> 4')
    end subroutine test_one_shot_counts

    !> B4: every position up to the cap is live.
    subroutine test_last_parameter_live()
        type(cache_t) :: cache
        real(kind = rk) :: base(21), moved(21), long(22), r0(8), r1(8)
        integer(kind = ik) :: status, pos
        call cache_init_s(cache, 21_ik, 65_ik, thetas, status)
        base = 0.0_rk
        base(1:5) = GS
        base(21) = 0.01_rk
        call cache_radius_grid_s(cache, base, r0, status)
        call assert_int_eq(status, 0_ik, 'B4 21-vector valid')
        moved = base
        moved(21) = moved(21) + 0.01_rk
        call cache_radius_grid_s(cache, moved, r1, status)
        call assert_true(.not. arrays_bits_eq_f(r0, r1), 'B4 alpha_20 is live')
        call cache_radius_grid_s(cache, GS, r0, status)
        do pos = 1_ik, 5_ik
            moved(1:5) = GS
            moved(pos) = moved(pos) + 0.01_rk
            call cache_radius_grid_s(cache, moved(1:5), r1, status)
            call assert_true(.not. arrays_bits_eq_f(r0, r1), 'B4 GS position is live')
        end do
        long = 0.0_rk
        long(1:5) = GS
        call cache_radius_grid_s(cache, long, r1, status)
        call assert_int_eq(status, SHAPE_ERROR_WRONG_PARAM_COUNT, 'B4 22-vector on the 21-cache -> 4')
    end subroutine test_last_parameter_live

    !> B5: the profile floor, both tiers.
    subroutine test_profile_floor()
        type(cache_t) :: cache
        real(kind = rk) :: z4(4), rho4(4), d4(4), z3(3), rho3(3), d3(3)
        integer(kind = ik) :: status
        call cache_init_s(cache, 5_ik, 65_ik, thetas, status)
        call cache_rho_profile_s(cache, GS, z4, rho4, d4, status)
        call assert_int_eq(status, 0_ik, 'B5 size 4 -> 0')
        call assert_true(is_zero_f(rho4(1)) .and. is_zero_f(rho4(4)) .and. rho4(2) > 0.0_rk &
                .and. rho4(3) > 0.0_rk, 'B5 size 4: poles and two interior nodes')
        call cache_rho_profile_s(cache, GS, z3, rho3, d3, status)
        call assert_int_eq(status, SHAPE_ERROR_INVALID_GRID, 'B5 size 3 -> 3 (cached)')
        call compute_rho_profile_standalone_s(GS, z3, rho3, d3, status)
        call assert_int_eq(status, SHAPE_ERROR_INVALID_GRID, 'B5 size 3 -> 3 (one-shot)')
        call cache_rho_profile_unchecked_s(cache, SEPARATED, z4, rho4, d4, status)
        call assert_int_eq(status, 0_ik, 'B5 SEPARATED size 4 unchecked -> 0')
        call assert_true(is_zero_f(rho4(2)) .and. is_zero_f(rho4(3)) .and. z4(2) > 0.0_rk &
                .and. z4(3) < 0.0_rk, 'B5 SEPARATED size 4: nodes 2 and 3 are the tips')
        call compute_rho_profile_unchecked_standalone_s(SEPARATED, z4, rho4, d4, status)
        call assert_int_eq(status, 0_ik, 'B5 SEPARATED size 4 unchecked one-shot -> 0')
    end subroutine test_profile_floor

    !> B6: output extents.
    subroutine test_extents()
        type(cache_t) :: cache
        real(kind = rk) :: r9(9), d9(9), r8(8), d8(8), d7(7), z(65), rho(65), drho(64)
        integer(kind = ik) :: status
        call cache_init_s(cache, 5_ik, 65_ik, thetas, status)
        call cache_radius_grid_s(cache, GS, r9, status)
        call assert_int_eq(status, SHAPE_ERROR_INVALID_GRID, 'B6 radius_grid n + 1 -> 3')
        call cache_radius_and_derivative_s(cache, GS, r9, d9, status)
        call assert_int_eq(status, SHAPE_ERROR_INVALID_GRID, 'B6 radius_and_derivative n + 1 -> 3')
        call cache_radius_and_derivative_s(cache, GS, r8, d7, status)
        call assert_int_eq(status, SHAPE_ERROR_INVALID_GRID, 'B6 radius_and_derivative dr short -> 3')
        call cache_radius_and_derivative_at_thetas_s(cache, GS, thetas, r9, d9, status)
        call assert_int_eq(status, SHAPE_ERROR_INVALID_GRID, 'B6 at_thetas size(radii) /= size(thetas) -> 3')
        call cache_rho_profile_s(cache, GS, z, rho, drho, status)
        call assert_int_eq(status, SHAPE_ERROR_INVALID_GRID, 'B6 size(drho_dz) /= size(z) -> 3')
        call compute_radius_grid_standalone_s(GS, 65_ik, thetas, r9, status)
        call assert_int_eq(status, SHAPE_ERROR_INVALID_GRID, 'B6 one-shot size(radii) /= size(thetas) -> 3')
        call cache_radius_and_derivative_s(cache, GS, r8, d8, status)
        call assert_int_eq(status, 0_ik, 'B6 matching extents -> 0')
    end subroutine test_extents

    !> B7: at_thetas theta counts.
    subroutine test_theta_counts()
        type(cache_t) :: cache
        real(kind = rk) :: one(1), r1(1), d1(1), empty(0), r0(0), d0(0)
        integer(kind = ik) :: status
        call cache_init_s(cache, 5_ik, 65_ik, thetas, status)
        one(1) = 1.0_rk
        call cache_radius_and_derivative_at_thetas_s(cache, GS, one, r1, d1, status)
        call assert_int_eq(status, 0_ik, 'B7 one theta -> 0')
        call cache_radius_and_derivative_at_thetas_s(cache, GS, empty, r0, d0, status)
        call assert_int_eq(status, SHAPE_ERROR_INVALID_GRID, 'B7 no theta -> 3')
    end subroutine test_theta_counts

    !> B8: a max_params = 3 cache.
    subroutine test_small_cache()
        type(cache_t) :: cache
        real(kind = rk), parameter :: SHORT(3) = [0.25_rk, 0.0_rk, 0.0_rk]
        real(kind = rk) :: radii(8)
        integer(kind = ik) :: status
        call cache_init_s(cache, 3_ik, 65_ik, thetas, status)
        call cache_radius_grid_s(cache, GS, radii, status)
        call assert_int_eq(status, SHAPE_ERROR_WRONG_PARAM_COUNT, 'B8 GS on a 3-cache -> 4')
        call cache_radius_grid_s(cache, SHORT, radii, status)
        call assert_int_eq(status, 0_ik, 'B8 [0.25, 0, 0] on a 3-cache -> 0')
    end subroutine test_small_cache

end program cassini_param_boundary_test
