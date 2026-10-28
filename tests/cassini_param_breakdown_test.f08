!> Where the formulas break down: the walk through alpha = 1, R -> 0+, the
!! toroidal boundary, an exact fold node, the Codex bracket point, the
!! 21-parameter cap, the margin edges, branch continuity through the whole
!! pipeline, and the scan floor.
program cassini_param_breakdown_test
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use precision_utilities_mod, only: ik, ikl, rk
    use mathematical_and_physical_constants_mod, only: PI_C
    use cassini_parameterization_mod, only: cache_t, cassini_geometry_t, cache_init_s, &
            cache_geometry_s, cache_radius_grid_s, cache_rho_profile_unchecked_s, &
            cache_rho_profile_s, STAR_CONVEXITY_MARGIN, CASSINI_ERROR_R_NONPOSITIVE, &
            CASSINI_ERROR_TOROIDAL, CASSINI_ERROR_NOT_STAR_CONVEX, CASSINI_ERROR_SEPARATED
    use cassini_parameterization_workers_mod, only: resolve_t, profile_nodes_s, series_at_s, &
            lemniscate_point_s
    use cassini_test_fixtures_mod, only: PINCHED_ASYM, GS, ASYM
    use cassini_test_reference_mod, only: test_resolve_s, reference_graded_volume_s, gl_thetas_s, &
            reference_theta_volume_s
    use test_utils_mod, only: assert_true, assert_int_eq, assert_close, assert_abs_close, &
            is_zero_f, test_summary
    implicit none

    real(kind = rk), parameter :: FOUR_THIRDS_PI = 4.0_rk * PI_C / 3.0_rk
    real(kind = rk) :: thetas(16)
    integer(kind = ik) :: k
    integer(kind = ikl) :: lcg_state

    lcg_state = 1971_ikl
    do k = 1_ik, 16_ik
        thetas(k) = PI_C * (real(k, rk) - 0.5_rk) / 16.0_rk
    end do

    call test_bare_walk()
    call test_asym_walk()
    call test_r_to_zero()
    call test_toroidal_boundary()
    call test_exact_fold_node()
    call test_codex_point()
    call test_cap_21()
    call test_margin_edges()
    call test_pipeline_continuity()
    call test_scan_floor(GS, 'GS')
    call test_scan_floor(ASYM, 'ASYM')
    call test_gate_at_the_floor()
    call test_summary()

contains

    function uniform_f(lo, hi) result(u)
        real(kind = rk), intent(in) :: lo, hi
        real(kind = rk) :: u
        lcg_state = mod(lcg_state * 48271_ikl, 2147483647_ikl)
        u = lo + (hi - lo) * real(lcg_state, rk) / 2147483647.0_rk
    end function uniform_f

    !> A geometry on a fresh cache.
    subroutine geometry_s(params, n_x_grid, geometry, status)
        real(kind = rk), intent(in) :: params(:)
        integer(kind = ik), intent(in) :: n_x_grid
        type(cassini_geometry_t), intent(out) :: geometry
        integer(kind = ik), intent(out) :: status
        type(cache_t) :: cache
        call cache_init_s(cache, size(params, kind = ik), n_x_grid, thetas, status)
        call cache_geometry_s(cache, params, geometry, status)
    end subroutine geometry_s

    !> The resolve's volume against graded panels in x.
    subroutine check_volume(params, tol, label)
        real(kind = rk), intent(in) :: params(:), tol
        character(len = *), intent(in) :: label
        type(resolve_t) :: res
        integer(kind = ik) :: status
        real(kind = rk) :: vbar, zm
        call test_resolve_s(params, res, status)
        call reference_graded_volume_s(params, res%n_active, res%eps, 64_ik, vbar, zm)
        call assert_close(FOUR_THIRDS_PI * res%c**3 / vbar, 1.0_rk, tol, label)
    end subroutine check_volume

    !> K1: the bare oval through alpha = 1.
    subroutine test_bare_walk()
        real(kind = rk), parameter :: WALK(5) = [1.0_rk - 1.0e-12_rk, 1.0_rk - 1.0e-15_rk, 1.0_rk, &
                1.0_rk + 1.0e-15_rk, 1.0_rk + 1.0e-12_rk]
        logical, parameter :: PINCH(5) = [.false., .false., .true., .false., .false.]
        logical, parameter :: SEP(5) = [.false., .false., .false., .true., .true.]
        type(cassini_geometry_t) :: geometry
        type(cache_t) :: cache
        real(kind = rk) :: bare(1), z(65), rho(65), drho(65)
        integer(kind = ik) :: status, j
        do j = 1_ik, 5_ik
            bare(1) = WALK(j)
            call geometry_s(bare, 65_ik, geometry, status)
            call assert_int_eq(status, 0_ik, 'K1 geometry resolves')
            call assert_true((geometry%is_pinched .eqv. PINCH(j)) .and. (geometry%is_separated .eqv. SEP(j)), &
                    'K1 flags along the walk')
            call check_volume(bare, 1.0e-13_rk, 'K1 volume along the walk')
        end do
        bare(1) = WALK(1)
        call geometry_s(bare, 65_ik, geometry, status)
        call cache_init_s(cache, 1_ik, 65_ik, thetas, status)
        call cache_rho_profile_s(cache, bare, z, rho, drho, status)
        call assert_close(rho(33) * geometry%c / sqrt(1.0e-12_rk), 1.0_rk, 1.0e-4_rk, &
                'K1 [1 - 1e-12] equator radius = sqrt(1e-12)/c (the hypot case)')
        bare(1) = WALK(5)
        call geometry_s(bare, 65_ik, geometry, status)
        call assert_close(geometry%gap * geometry%c / (2.0_rk * sqrt(1.0e-12_rk)), 1.0_rk, 1.0e-4_rk, &
                'K1 [1 + 1e-12] gap = 2 sqrt(1e-12)/c')
    end subroutine test_bare_walk

    !> K2: the same walk on PINCHED_ASYM's alpha_n.
    subroutine test_asym_walk()
        real(kind = rk), parameter :: DELTA(5) = [-1.0e-12_rk, -1.0e-15_rk, 0.0_rk, 1.0e-15_rk, 1.0e-12_rk]
        type(cassini_geometry_t) :: geometry
        type(resolve_t) :: res
        real(kind = rk) :: params(5)
        integer(kind = ik) :: status, j
        do j = 1_ik, 5_ik
            params = PINCHED_ASYM
            params(1) = 1.0_rk + DELTA(j)
            call geometry_s(params, 65_ik, geometry, status)
            call test_resolve_s(params, res, status)
            call assert_true((geometry%is_pinched .eqv. j == 3_ik) .and. &
                    (geometry%is_separated .eqv. j > 3_ik), 'K2 flags along the asymmetric walk')
            call assert_true((res%r_zero**2 - res%eps < 0.0_rk .eqv. j > 3_ik) .and. &
                    (is_zero_f(res%r_zero**2 - res%eps) .eqv. j == 3_ik), 'K2 eps - R0^2 changes sign at 1')
            call check_volume(params, 1.0e-13_rk, 'K2 volume along the asymmetric walk')
        end do
    end subroutine test_asym_walk

    !> K3: R(1) -> 0+.
    subroutine test_r_to_zero()
        real(kind = rk), parameter :: ALPHA1(3) = [-0.999_rk, -0.9999_rk, -1.0_rk]
        integer(kind = ik), parameter :: WANT(3) = [0_ik, 0_ik, CASSINI_ERROR_R_NONPOSITIVE]
        type(cassini_geometry_t) :: geometry
        real(kind = rk) :: params(2)
        integer(kind = ik) :: status, j
        do j = 1_ik, 3_ik
            params = [0.5_rk, 0.0_rk]
            params(2) = ALPHA1(j)
            call geometry_s(params, 129_ik, geometry, status)
            call assert_int_eq(status, WANT(j), 'K3 R -> 0+ status')
            if (status == 0_ik) then
                call assert_true(geometry%r_north > 0.0_rk .and. geometry%r_south > 0.0_rk &
                        .and. ieee_is_finite(geometry%r_north), 'K3 poles finite and positive')
                write(*, '(A,F8.4,A,2ES14.6,A,ES12.4)') 'K3 alpha1 = ', ALPHA1(j), ' poles ', &
                        geometry%r_north, geometry%r_south, ' d_min ', geometry%d_min
            end if
        end do
    end subroutine test_r_to_zero

    !> K4: s -> -R(+-1)^2, the toroidal boundary of the bare oval at alpha = -1.
    subroutine test_toroidal_boundary()
        real(kind = rk), parameter :: WALK(4) = [-1.0_rk + 1.0e-6_rk, -1.0_rk + 1.0e-9_rk, -1.0_rk, &
                -1.0_rk - 1.0e-9_rk]
        integer(kind = ik), parameter :: WANT(4) = [0_ik, 0_ik, CASSINI_ERROR_TOROIDAL, CASSINI_ERROR_TOROIDAL]
        type(cassini_geometry_t) :: geometry
        real(kind = rk) :: bare(1)
        integer(kind = ik) :: status, j
        do j = 1_ik, 4_ik
            bare(1) = WALK(j)
            call geometry_s(bare, 129_ik, geometry, status)
            call assert_int_eq(status, WANT(j), 'K4 toroidal boundary status')
            if (j <= 2_ik) then
                call assert_close(geometry%r_north * geometry%c / sqrt(1.0_rk + WALK(j)), 1.0_rk, 1.0e-6_rk, &
                        'K4 pole radius sqrt(1 + alpha)/c')
                write(*, '(A,ES10.2,A,ES12.4)') 'K4 1 + alpha = ', 1.0_rk + WALK(j), ' d_min ', geometry%d_min
            end if
        end do
    end subroutine test_toroidal_boundary

    !> K5: bare -1/sqrt(2) has its fold boundary on a node of the 65-scan.
    subroutine test_exact_fold_node()
        type(cassini_geometry_t) :: geometry
        type(cache_t) :: cache
        real(kind = rk) :: bare(1), radii(16)
        integer(kind = ik) :: status
        bare(1) = -1.0_rk / sqrt(2.0_rk)
        call geometry_s(bare, 65_ik, geometry, status)
        call assert_int_eq(status, 0_ik, 'K5 geometry resolves')
        call assert_true(geometry%is_folded .and. ieee_is_finite(geometry%d_min), 'K5 folded, d_min finite')
        call cache_init_s(cache, 1_ik, 65_ik, thetas, status)
        call cache_radius_grid_s(cache, bare, radii, status)
        call assert_true((status == 0_ik .eqv. geometry%d_min >= STAR_CONVEXITY_MARGIN), &
                'K5 the gate decides by d_min, not by the fold')
        write(*, '(A,I0,A,ES12.4)') 'K5 R(theta) status ', status, ', d_min ', geometry%d_min
    end subroutine test_exact_fold_node

    !> K6: the bare oval at s = -3/4 (the Codex bracket point at u = 2/3).
    subroutine test_codex_point()
        type(cassini_geometry_t) :: geometry
        type(cache_t) :: cache
        real(kind = rk) :: bare(1), radii(16)
        integer(kind = ik) :: status
        bare(1) = -0.75_rk
        call geometry_s(bare, 2049_ik, geometry, status)
        call assert_int_eq(status, 0_ik, 'K6 [-0.75] geometry resolves')
        call cache_init_s(cache, 1_ik, 2049_ik, thetas, status)
        call cache_radius_grid_s(cache, bare, radii, status)
        call assert_true(status == 0_ik .or. status == CASSINI_ERROR_NOT_STAR_CONVEX, 'K6 R(theta) 0 or 102')
        call assert_true((status == 0_ik .eqv. geometry%d_min >= STAR_CONVEXITY_MARGIN) &
                .and. all(ieee_is_finite(radii)), 'K6 verdict by d_min, finite radii')
    end subroutine test_codex_point

    !> K7: 50 random 21-vectors on a max_params = 21 cache.
    subroutine test_cap_21()
        type(cache_t) :: cache
        type(cassini_geometry_t) :: geometry
        real(kind = rk) :: params(21), gl_t(200), gl_w(200), radii(200), vol, zcm, worst
        integer(kind = ik) :: status, trial, n, counts(0:3)
        call gl_thetas_s(200_ik, gl_t, gl_w)
        call cache_init_s(cache, 21_ik, 257_ik, gl_t, status)
        counts = 0_ik
        worst = 0.0_rk
        do trial = 1_ik, 50_ik
            params(1) = uniform_f(0.0_rk, 0.9_rk)
            do n = 2_ik, 21_ik
                params(n) = uniform_f(-0.02_rk, 0.02_rk) / real(n - 1_ik, rk)
            end do
            call cache_geometry_s(cache, params, geometry, status)
            call assert_int_eq(status, 0_ik, 'K7 geometry resolves (never 100/101)')
            call cache_radius_grid_s(cache, params, radii, status)
            if (status == 0_ik) then
                counts(0) = counts(0) + 1_ik
                call reference_theta_volume_s(gl_t, gl_w, radii, vol, zcm)
                worst = max(worst, abs(vol / FOUR_THIRDS_PI - 1.0_rk))
            else if (status == CASSINI_ERROR_NOT_STAR_CONVEX) then
                counts(2) = counts(2) + 1_ik
            else if (status == CASSINI_ERROR_SEPARATED) then
                counts(3) = counts(3) + 1_ik
            else
                counts(1) = counts(1) + 1_ik
            end if
        end do
        write(*, '(A,4I4,A,ES10.2)') 'K7 R(theta) counts 0 / other / 102 / 103:', counts, &
                '; worst volume ', worst
        call assert_int_eq(counts(1), 0_ik, 'K7 R(theta) only 0, 102 or 103')
        call assert_true(worst <= 1.0e-11_rk, 'K7 theta-frame volume 4 pi/3')
    end subroutine test_cap_21

    !> K8: bisect alpha to the margin crossing on STAR_FAIL's alpha_n.
    subroutine test_margin_edges()
        type(cassini_geometry_t) :: geometry
        type(cache_t) :: cache
        real(kind = rk) :: params(5), lo, hi, mid, radii(16)
        integer(kind = ik) :: status, step
        params = [0.0_rk, -0.2_rk, 0.0_rk, -0.05_rk, -0.05_rk]
        lo = 0.80_rk
        hi = 0.95_rk
        do step = 1_ik, 30_ik
            mid = 0.5_rk * (lo + hi)
            params(1) = mid
            call geometry_s(params, 129_ik, geometry, status)
            if (geometry%d_min >= STAR_CONVEXITY_MARGIN) then
                lo = mid
            else
                hi = mid
            end if
        end do
        write(*, '(A,F12.8)') 'K8 margin crossing at alpha = ', lo
        call cache_init_s(cache, 5_ik, 129_ik, thetas, status)
        params(1) = lo - 1.0e-3_rk
        call cache_radius_grid_s(cache, params, radii, status)
        call assert_int_eq(status, 0_ik, 'K8 just inside the margin -> 0')
        params(1) = lo + 1.0e-3_rk
        call cache_radius_grid_s(cache, params, radii, status)
        call assert_int_eq(status, CASSINI_ERROR_NOT_STAR_CONVEX, 'K8 just outside the margin -> 102')
    end subroutine test_margin_edges

    !> K9: no jump anywhere along a dense ASYM profile (|x| <= 0.9: the
    !! profile's own resolution, not a branch, sets the step near the poles
    !! and the dimple fold).
    subroutine test_pipeline_continuity()
        type(cache_t) :: cache
        type(resolve_t) :: res
        real(kind = rk) :: z(2001), rho(2001), drho_dz(2001), x(2001), lobe(2001)
        real(kind = rk) :: dx_rho(2001), dx_z(2001), r, dr, rb, zb, worst, worst_d
        integer(kind = ik) :: status, j
        call cache_init_s(cache, 5_ik, 129_ik, thetas, status)
        call cache_rho_profile_unchecked_s(cache, ASYM, z, rho, drho_dz, status)
        call test_resolve_s(ASYM, res, status)
        call profile_nodes_s(res, 2001_ik, x, lobe)
        do j = 1_ik, 2001_ik
            call series_at_s(ASYM, res%n_active, x(j), r, dr)
            call lemniscate_point_s(x(j), r, dr, res%eps, rb, zb, dx_rho(j), dx_z(j))
        end do
        dx_rho = dx_rho / res%c
        dx_z = dx_z / res%c
        worst = 0.0_rk
        worst_d = 0.0_rk
        do j = 2_ik, 2001_ik
            if (abs(x(j)) > 0.9_rk .or. abs(x(j - 1_ik)) > 0.9_rk) cycle
            worst = max(worst, abs(rho(j) - rho(j - 1_ik)), abs(z(j) - z(j - 1_ik)), &
                    abs(dx_rho(j) - dx_rho(j - 1_ik)), abs(dx_z(j) - dx_z(j - 1_ik)))
            ! relative: drho_dz is a steep smooth function toward ASYM's south
            ! fold (2.5 -> 2.7 by 0.02 per node near x = -0.84)
            if (dx_z(j) > 1.0e-3_rk .and. dx_z(j - 1_ik) > 1.0e-3_rk) &
                    worst_d = max(worst_d, abs(drho_dz(j) - drho_dz(j - 1_ik)) / (1.0_rk + abs(drho_dz(j))))
        end do
        call assert_true(worst <= 2.0e-2_rk, 'K9 rho, z, drho_dx, dz_dx continuous along the profile')
        call assert_true(worst_d <= 2.0e-2_rk, 'K9 drho_dz continuous where dz_dx > 1e-3')
    end subroutine test_pipeline_continuity

    !> K10: d_min and R(theta) barely move with the scan resolution.
    subroutine test_scan_floor(params, name)
        real(kind = rk), intent(in) :: params(:)
        character(len = *), intent(in) :: name
        integer(kind = ik), parameter :: GRIDS(3) = [65_ik, 129_ik, 2049_ik]
        type(cassini_geometry_t) :: geometry
        type(cache_t) :: cache
        real(kind = rk) :: d(3), radii(16, 3)
        integer(kind = ik) :: status, j
        do j = 1_ik, 3_ik
            call geometry_s(params, GRIDS(j), geometry, status)
            d(j) = geometry%d_min
            call cache_init_s(cache, 5_ik, GRIDS(j), thetas, status)
            call cache_radius_grid_s(cache, params, radii(:, j), status)
        end do
        call assert_true(maxval(d) - minval(d) <= 1.0e-3_rk, 'K10 d_min vs scan size ' // name)
        call assert_true(maxval(abs(radii(:, 1) - radii(:, 3))) <= 1.0e-12_rk &
                .and. maxval(abs(radii(:, 2) - radii(:, 3))) <= 1.0e-12_rk, 'K10 radius_grid vs scan size ' // name)
    end subroutine test_scan_floor

    !> K11: between scan nodes d can dip below its sampled minimum. The box
    !! corner [0.95, 0.25, 0.1, 0.25, 0.25] samples d_min = +4.3e-3 on the
    !! 65-node floor scan but is not star-convex (dense d_min = -1.25e-2):
    !! the gate must say 102 at every scan size, and d_min must not depend on it.
    subroutine test_gate_at_the_floor()
        real(kind = rk), parameter :: CORNER(5) = [0.95_rk, 0.25_rk, 0.1_rk, 0.25_rk, 0.25_rk]
        integer(kind = ik), parameter :: GRIDS(3) = [65_ik, 129_ik, 2049_ik]
        type(cassini_geometry_t) :: geometry
        type(cache_t) :: cache
        real(kind = rk) :: radii(16), d(3)
        integer(kind = ik) :: status, j
        do j = 1_ik, 3_ik
            call cache_init_s(cache, 5_ik, GRIDS(j), thetas, status)
            call cache_radius_grid_s(cache, CORNER, radii, status)
            call assert_int_eq(status, CASSINI_ERROR_NOT_STAR_CONVEX, 'K11 non-star-convex corner -> 102')
            call cache_geometry_s(cache, CORNER, geometry, status)
            d(j) = geometry%d_min
        end do
        call assert_true(d(1) < 0.0_rk, 'K11 d_min < 0 on the floor scan')
        call assert_true(maxval(d) - minval(d) <= 1.0e-6_rk, 'K11 d_min independent of the scan size')
    end subroutine test_gate_at_the_floor

end program cassini_param_breakdown_test
