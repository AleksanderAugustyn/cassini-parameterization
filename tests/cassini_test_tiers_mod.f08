!> Runs every output of one tier into a result record, for the error,
!! boundary, equivalence and statelessness suites. Output slots, in order:
!! radius_grid, radius_and_derivative, at_thetas, shape, rho_profile,
!! rho_profile_unchecked, neck, geometry, quadrupole.
module cassini_test_tiers_mod

    use precision_utilities_mod, only: ik, rk
    use cassini_parameterization_mod, only: cache_t, cassini_geometry_t, cache_init_s, &
            cache_radius_grid_s, cache_radius_and_derivative_s, &
            cache_radius_and_derivative_at_thetas_s, cache_shape_s, cache_rho_profile_s, &
            cache_rho_profile_unchecked_s, cache_neck_s, cache_geometry_s, cache_quadrupole_s, &
            compute_radius_grid_standalone_s, compute_radius_and_derivative_standalone_s, &
            compute_shape_standalone_s, compute_rho_profile_standalone_s, &
            compute_rho_profile_unchecked_standalone_s, compute_neck_standalone_s, &
            compute_geometry_standalone_s, compute_quadrupole_standalone_s
    use test_utils_mod, only: bits_eq_f, arrays_bits_eq_f, is_zero_f, all_zero_f

    implicit none

    private

    integer(kind = ik), parameter, public :: N_OUT = 9_ik
    integer(kind = ik), parameter, public :: OUT_RADIUS_GRID = 1_ik, OUT_RADIUS_DERIV = 2_ik, &
            OUT_AT_THETAS = 3_ik, OUT_SHAPE = 4_ik, OUT_PROFILE = 5_ik, OUT_PROFILE_UNCHECKED = 6_ik, &
            OUT_NECK = 7_ik, OUT_GEOMETRY = 8_ik, OUT_QUADRUPOLE = 9_ik

    type, public :: result_t
        integer(kind = ik) :: status(N_OUT) = -1_ik
        real(kind = rk), allocatable :: radii_g(:), radii_d(:), dr_d(:), radii_a(:), dr_a(:)
        real(kind = rk) :: z_shift = 0.0_rk, r_north = 0.0_rk, r_south = 0.0_rk
        real(kind = rk), allocatable :: z(:), rho(:), drho(:), zu(:), rhou(:), drhou(:)
        real(kind = rk) :: z_neck = 0.0_rk, rho_neck = 0.0_rk
        logical :: found = .false.
        type(cassini_geometry_t) :: geom
        real(kind = rk) :: q2 = 0.0_rk
    end type result_t

    public :: run_cached, run_one_shot, results_bits_eq, geometries_bits_eq, &
            output_zero_filled, build_cache

contains

    !> A cache on the given thetas; status returned.
    subroutine build_cache(max_params, n_x_grid, thetas, cache, status)
        integer(kind = ik), intent(in) :: max_params, n_x_grid
        real(kind = rk), intent(in) :: thetas(:)
        type(cache_t), intent(out) :: cache
        integer(kind = ik), intent(out) :: status
        call cache_init_s(cache, max_params, n_x_grid, thetas, status)
    end subroutine build_cache

    subroutine allocate_result(res, n_thetas, n_z)
        type(result_t), intent(out) :: res
        integer(kind = ik), intent(in) :: n_thetas, n_z
        allocate(res%radii_g(n_thetas), res%radii_d(n_thetas), res%dr_d(n_thetas), &
                res%radii_a(n_thetas), res%dr_a(n_thetas))
        allocate(res%z(n_z), res%rho(n_z), res%drho(n_z), res%zu(n_z), res%rhou(n_z), res%drhou(n_z))
    end subroutine allocate_result

    !> Every output of the cached tier on one cache; at_thetas at `thetas`
    !! (which must have the cache's theta count for the grid slots).
    subroutine run_cached(cache, params, thetas, n_z, res)
        type(cache_t), intent(in) :: cache
        real(kind = rk), intent(in) :: params(:), thetas(:)
        integer(kind = ik), intent(in) :: n_z
        type(result_t), intent(out) :: res
        call allocate_result(res, size(thetas, kind = ik), n_z)
        call cache_radius_grid_s(cache, params, res%radii_g, res%status(OUT_RADIUS_GRID))
        call cache_radius_and_derivative_s(cache, params, res%radii_d, res%dr_d, &
                res%status(OUT_RADIUS_DERIV))
        call cache_radius_and_derivative_at_thetas_s(cache, params, thetas, res%radii_a, res%dr_a, &
                res%status(OUT_AT_THETAS))
        call cache_shape_s(cache, params, res%z_shift, res%r_north, res%r_south, res%status(OUT_SHAPE))
        call cache_rho_profile_s(cache, params, res%z, res%rho, res%drho, res%status(OUT_PROFILE))
        call cache_rho_profile_unchecked_s(cache, params, res%zu, res%rhou, res%drhou, &
                res%status(OUT_PROFILE_UNCHECKED))
        call cache_neck_s(cache, params, res%z_neck, res%rho_neck, res%found, res%status(OUT_NECK))
        call cache_geometry_s(cache, params, res%geom, res%status(OUT_GEOMETRY))
        call cache_quadrupole_s(cache, params, res%q2, res%status(OUT_QUADRUPOLE))
    end subroutine run_cached

    !> Every output of the one-shot tier; the at_thetas slot holds the one-shot
    !! radius_and_derivative at `thetas` (there is no one-shot at-thetas form).
    subroutine run_one_shot(params, n_x_grid, thetas, n_z, res)
        real(kind = rk), intent(in) :: params(:), thetas(:)
        integer(kind = ik), intent(in) :: n_x_grid, n_z
        type(result_t), intent(out) :: res
        call allocate_result(res, size(thetas, kind = ik), n_z)
        call compute_radius_grid_standalone_s(params, n_x_grid, thetas, res%radii_g, &
                res%status(OUT_RADIUS_GRID))
        call compute_radius_and_derivative_standalone_s(params, n_x_grid, thetas, res%radii_d, &
                res%dr_d, res%status(OUT_RADIUS_DERIV))
        call compute_radius_and_derivative_standalone_s(params, n_x_grid, thetas, res%radii_a, &
                res%dr_a, res%status(OUT_AT_THETAS))
        call compute_shape_standalone_s(params, n_x_grid, res%z_shift, res%r_north, res%r_south, &
                res%status(OUT_SHAPE))
        call compute_rho_profile_standalone_s(params, res%z, res%rho, res%drho, res%status(OUT_PROFILE))
        call compute_rho_profile_unchecked_standalone_s(params, res%zu, res%rhou, res%drhou, &
                res%status(OUT_PROFILE_UNCHECKED))
        call compute_neck_standalone_s(params, n_x_grid, res%z_neck, res%rho_neck, res%found, &
                res%status(OUT_NECK))
        call compute_geometry_standalone_s(params, n_x_grid, res%geom, res%status(OUT_GEOMETRY))
        call compute_quadrupole_standalone_s(params, res%q2, res%status(OUT_QUADRUPOLE))
    end subroutine run_one_shot

    logical function geometries_bits_eq(a, b) result(same)
        type(cassini_geometry_t), intent(in) :: a, b
        same = bits_eq_f(a%eps, b%eps) .and. bits_eq_f(a%c, b%c) .and. bits_eq_f(a%z_shift, b%z_shift) &
                .and. bits_eq_f(a%r_north, b%r_north) .and. bits_eq_f(a%r_south, b%r_south) &
                .and. bits_eq_f(a%gap, b%gap) .and. bits_eq_f(a%d_min, b%d_min) &
                .and. (a%is_separated .eqv. b%is_separated) .and. (a%is_pinched .eqv. b%is_pinched) &
                .and. (a%is_folded .eqv. b%is_folded)
    end function geometries_bits_eq

    !> Every status and every output bit equal.
    logical function results_bits_eq(a, b) result(same)
        type(result_t), intent(in) :: a, b
        same = all(a%status == b%status) &
                .and. arrays_bits_eq_f(a%radii_g, b%radii_g) &
                .and. arrays_bits_eq_f(a%radii_d, b%radii_d) .and. arrays_bits_eq_f(a%dr_d, b%dr_d) &
                .and. arrays_bits_eq_f(a%radii_a, b%radii_a) .and. arrays_bits_eq_f(a%dr_a, b%dr_a) &
                .and. bits_eq_f(a%z_shift, b%z_shift) &
                .and. bits_eq_f(a%r_north, b%r_north) .and. bits_eq_f(a%r_south, b%r_south) &
                .and. arrays_bits_eq_f(a%z, b%z) .and. arrays_bits_eq_f(a%rho, b%rho) &
                .and. arrays_bits_eq_f(a%drho, b%drho) .and. arrays_bits_eq_f(a%zu, b%zu) &
                .and. arrays_bits_eq_f(a%rhou, b%rhou) .and. arrays_bits_eq_f(a%drhou, b%drhou) &
                .and. bits_eq_f(a%z_neck, b%z_neck) .and. bits_eq_f(a%rho_neck, b%rho_neck) &
                .and. (a%found .eqv. b%found) .and. geometries_bits_eq(a%geom, b%geom) &
                .and. bits_eq_f(a%q2, b%q2)
    end function results_bits_eq

    !> Whether output slot `which` of a result is zero-filled over its full extent.
    logical function output_zero_filled(res, which) result(zeroed)
        type(result_t), intent(in) :: res
        integer(kind = ik), intent(in) :: which
        type(cassini_geometry_t) :: g
        select case (which)
        case (OUT_RADIUS_GRID)
            zeroed = all_zero_f(res%radii_g)
        case (OUT_RADIUS_DERIV)
            zeroed = all_zero_f(res%radii_d) .and. all_zero_f(res%dr_d)
        case (OUT_AT_THETAS)
            zeroed = all_zero_f(res%radii_a) .and. all_zero_f(res%dr_a)
        case (OUT_SHAPE)
            zeroed = is_zero_f(res%z_shift) .and. is_zero_f(res%r_north) .and. is_zero_f(res%r_south)
        case (OUT_PROFILE)
            zeroed = all_zero_f(res%z) .and. all_zero_f(res%rho) .and. all_zero_f(res%drho)
        case (OUT_PROFILE_UNCHECKED)
            zeroed = all_zero_f(res%zu) .and. all_zero_f(res%rhou) .and. all_zero_f(res%drhou)
        case (OUT_NECK)
            zeroed = is_zero_f(res%z_neck) .and. is_zero_f(res%rho_neck) .and. .not. res%found
        case (OUT_GEOMETRY)
            g = res%geom
            zeroed = is_zero_f(g%eps) .and. is_zero_f(g%c) .and. is_zero_f(g%z_shift) &
                    .and. is_zero_f(g%r_north) .and. is_zero_f(g%r_south) .and. is_zero_f(g%gap) &
                    .and. is_zero_f(g%d_min) .and. .not. g%is_separated .and. .not. g%is_pinched &
                    .and. .not. g%is_folded
        case default
            zeroed = is_zero_f(res%q2)
        end select
    end function output_zero_filled

end module cassini_test_tiers_mod
