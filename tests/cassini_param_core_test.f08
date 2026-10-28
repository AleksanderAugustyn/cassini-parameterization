!> The final-frame curve: the profile node sets and conventions, the
!! product identity and the polar line on the profile, the star margin and
!! the fold flag on the scan.
program cassini_param_core_test
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use precision_utilities_mod, only: ik, rk
    use cassini_parameterization_workers_mod, only: resolve_t, series_at_s, lemniscate_point_s, &
            classification_numbers_s, star_margin_s, profile_nodes_s, profile_s
    use cassini_test_fixtures_mod, only: SPHERE, BARE_05, GS, SADDLE2, ASYM, DIMPLED, FOLD_BOX, &
            TWO_NECKS, STAR_FAIL, PINCHED, SEPARATED
    use cassini_test_reference_mod, only: test_resolve_s, test_scan_s
    use test_utils_mod, only: assert_true, assert_int_eq, assert_abs_close, is_zero_f, all_zero_f, &
            test_summary
    implicit none

    call test_sphere_profile()
    call test_product_identity(SADDLE2, 'SADDLE2')
    call test_product_identity(ASYM, 'ASYM')
    call test_product_identity(DIMPLED, 'DIMPLED')
    call test_product_identity(FOLD_BOX, 'FOLD_BOX')
    call test_polar_line()
    call test_margin(SPHERE, 1.0_rk, 1.0e-14_rk, .false., 'SPHERE')
    call test_margin(BARE_05, 0.79698_rk, 1.0e-4_rk, .false., 'BARE_05')
    call test_margin(GS, 0.95103_rk, 1.0e-4_rk, .false., 'GS')
    call test_margin(SADDLE2, 0.74856_rk, 1.0e-4_rk, .false., 'SADDLE2')
    call test_margin(ASYM, 0.42733_rk, 1.0e-4_rk, .true., 'ASYM')
    ! d_min is the refined minimum (the plan's 0.62747, 0.24313 and 0.55891
    ! were 129-node samples; TWO_NECKS's minimum sits between scan nodes on
    ! its neck bump, 1.6e-3 below its sample)
    call test_margin(DIMPLED, 0.62734_rk, 1.0e-4_rk, .true., 'DIMPLED')
    call test_margin(TWO_NECKS, 0.24156_rk, 1.0e-4_rk, .false., 'TWO_NECKS')
    call test_margin(STAR_FAIL, -0.03960_rk, 1.0e-4_rk, .false., 'STAR_FAIL')
    call test_fold_box_folded()
    call test_margin_split()
    call test_profile_derivative(ASYM, 'ASYM')
    call test_profile_derivative(DIMPLED, 'DIMPLED')
    call test_two_lobes()
    call test_fold_boundary_node()
    call test_odd_grid_equator()
    call test_summary()

contains

    !> C1: the sphere's profile is the unit circle, poles exact, z descending.
    subroutine test_sphere_profile()
        type(resolve_t) :: res
        integer(kind = ik) :: status, j
        real(kind = rk) :: z(101), rho(101), drho_dz(101), worst
        logical :: descending
        call test_resolve_s(SPHERE, res, status)
        call profile_s(SPHERE, res, 101_ik, z, rho, drho_dz)
        worst = 0.0_rk
        descending = .true.
        do j = 2_ik, 100_ik
            worst = max(worst, abs(rho(j)**2 - (1.0_rk - z(j)**2)))
        end do
        do j = 1_ik, 100_ik
            if (.not. (z(j + 1_ik) < z(j))) descending = .false.
        end do
        call assert_true(worst <= 1.0e-14_rk, 'C1 sphere rho^2 = 1 - z^2')
        call assert_true(is_zero_f(rho(1)) .and. is_zero_f(rho(101)), 'C1 sphere pole rho = 0')
        call assert_true(is_zero_f(drho_dz(1)) .and. is_zero_f(drho_dz(101)), 'C1 sphere pole drho_dz = 0')
        ! 2e-15: c carries the 512-node quadrature's summation (5 ulps in Debug)
        call assert_abs_close(z(1), 1.0_rk, 2.0e-15_rk, 'C1 sphere z(1) = 1')
        call assert_abs_close(z(101), -1.0_rk, 2.0e-15_rk, 'C1 sphere z(101) = -1')
        call assert_true(descending, 'C1 sphere z descends')
    end subroutine test_sphere_profile

    !> C2: c^2 rho |z + z_shift| = R(x)^2 |x| sqrt(1 - x^2) at every profile node.
    subroutine test_product_identity(params, name)
        real(kind = rk), intent(in) :: params(:)
        character(len = *), intent(in) :: name
        type(resolve_t) :: res
        integer(kind = ik) :: status, j
        real(kind = rk) :: x(257), lobe(257), z(257), rho(257), drho_dz(257), r, dr, worst
        call test_resolve_s(params, res, status)
        call profile_nodes_s(res, 257_ik, x, lobe)
        call profile_s(params, res, 257_ik, z, rho, drho_dz)
        worst = 0.0_rk
        do j = 1_ik, 257_ik
            call series_at_s(params, res%n_active, x(j), r, dr)
            worst = max(worst, abs(res%c**2 * rho(j) * abs(z(j) + res%z_shift) &
                    - r**2 * abs(x(j)) * sqrt((1.0_rk - x(j)) * (1.0_rk + x(j)))))
        end do
        call assert_true(worst <= 1.0e-13_rk, 'C2 product identity on the profile ' // name)
    end subroutine test_product_identity

    !> C3: at alpha = alpha* the map is polar coordinates.
    subroutine test_polar_line()
        type(resolve_t) :: res
        integer(kind = ik) :: status, j
        real(kind = rk) :: params(5), r_plus, r_minus, r_zero, pole_plus, pole_minus, equator, sp
        real(kind = rk) :: x(257), lobe(257), z(257), rho(257), drho_dz(257), r, dr, worst
        params = SADDLE2
        call classification_numbers_s(params, 5_ik, 0.0_rk, r_plus, r_minus, r_zero, &
                pole_plus, pole_minus, equator)
        sp = r_plus**2 + r_minus**2
        params(1) = (sp - 2.0_rk * r_zero**2) / (sp + 2.0_rk * r_zero**2)
        call test_resolve_s(params, res, status)
        call profile_nodes_s(res, 257_ik, x, lobe)
        call profile_s(params, res, 257_ik, z, rho, drho_dz)
        worst = 0.0_rk
        do j = 1_ik, 257_ik
            call series_at_s(params, res%n_active, x(j), r, dr)
            worst = max(worst, abs(res%c * sqrt(rho(j)**2 + (z(j) + res%z_shift)**2) - r))
        end do
        call assert_true(worst <= 1.0e-13_rk, 'C3 polar line: c |r| = R(x)')
    end subroutine test_polar_line

    !> C4: d_min (refined) on the 129-node scan and the fold flag; C5: against 20001 nodes.
    subroutine test_margin(params, d_want, tol, folded, name)
        real(kind = rk), intent(in) :: params(:), d_want, tol
        logical, intent(in) :: folded
        character(len = *), intent(in) :: name
        type(resolve_t) :: res
        integer(kind = ik) :: status
        real(kind = rk) :: d_min, d_dense
        logical :: is_folded
        call test_resolve_s(params, res, status)
        call margin_s(params, res, 129_ik, d_min, is_folded)
        call assert_abs_close(d_min, d_want, tol, 'C4 d_min ' // name)
        call assert_true(is_folded .eqv. folded, 'C4 is_folded ' // name)
        call margin_s(params, res, 20001_ik, d_dense, is_folded)
        ! refined at both resolutions: the 129-node margin is the 20001-node one
        call assert_true(abs(d_min - d_dense) <= 1.0e-6_rk, 'C5 refined margin vs 20001 ' // name)
    end subroutine test_margin

    subroutine margin_s(params, res, n, d_min, is_folded)
        real(kind = rk), intent(in) :: params(:)
        type(resolve_t), intent(in) :: res
        integer(kind = ik), intent(in) :: n
        real(kind = rk), intent(out) :: d_min
        logical, intent(out) :: is_folded
        real(kind = rk) :: x(n), rho(n), z(n), drho_dx(n), dz_dx(n)
        call test_scan_s(params, res, n, x, rho, z, drho_dx, dz_dx)
        call star_margin_s(params, res, x, rho, z, drho_dx, dz_dx, d_min, is_folded)
    end subroutine margin_s

    !> C4: FOLD_BOX folds; separated d_min = 0; the pinch node is skipped.
    subroutine test_fold_box_folded()
        type(resolve_t) :: res
        integer(kind = ik) :: status
        real(kind = rk) :: d_min
        logical :: is_folded
        call test_resolve_s(FOLD_BOX, res, status)
        call margin_s(FOLD_BOX, res, 129_ik, d_min, is_folded)
        call assert_true(is_folded, 'C4 is_folded FOLD_BOX')
        call assert_abs_close(d_min, 0.55852_rk, 1.0e-4_rk, 'C4 d_min FOLD_BOX')
        call test_resolve_s(SEPARATED, res, status)
        call margin_s(SEPARATED, res, 129_ik, d_min, is_folded)
        call assert_true(is_zero_f(d_min), 'C4 SEPARATED d_min = 0')
        call test_resolve_s(PINCHED, res, status)
        call margin_s(PINCHED, res, 65_ik, d_min, is_folded)
        call assert_true(ieee_is_finite(d_min), 'C4 PINCHED d_min finite on the 65-node scan')
        call assert_true(.not. is_folded, 'C4 PINCHED not folded')
    end subroutine test_fold_box_folded

    !> C5 on a split-sensitive shape: STAR_FAIL is the margin failure.
    subroutine test_margin_split()
        type(resolve_t) :: res
        integer(kind = ik) :: status
        real(kind = rk) :: d_min
        logical :: is_folded
        call test_resolve_s(STAR_FAIL, res, status)
        call margin_s(STAR_FAIL, res, 129_ik, d_min, is_folded)
        call assert_true(d_min < 0.0_rk, 'C4 STAR_FAIL d_min < 0')
    end subroutine test_margin_split

    !> C6: drho_dz against the 3-point secant of the profile where dz_dx >
    !! 5e-2 (within ~2e-2 of a fold the secant itself is off by up to 6e-3 on
    !! ASYM), and against a fine secant (h = 1e-7 in x) wherever dz_dx > 1e-3.
    subroutine test_profile_derivative(params, name)
        real(kind = rk), intent(in) :: params(:)
        character(len = *), intent(in) :: name
        real(kind = rk), parameter :: H = 1.0e-7_rk
        type(resolve_t) :: res
        integer(kind = ik) :: status, j
        real(kind = rk) :: x(1001), lobe(1001), z(1001), rho(1001), drho_dz(1001)
        real(kind = rk) :: r, dr, rb, zb, drb, dzb, secant, rho_p, z_p, rho_m, z_m
        logical :: ok, ok_fine
        call test_resolve_s(params, res, status)
        call profile_nodes_s(res, 1001_ik, x, lobe)
        call profile_s(params, res, 1001_ik, z, rho, drho_dz)
        ok = .true.
        ok_fine = .true.
        do j = 5_ik, 997_ik
            call series_at_s(params, res%n_active, x(j), r, dr)
            call lemniscate_point_s(x(j), r, dr, res%eps, rb, zb, drb, dzb)
            if (.not. (dzb / res%c > 1.0e-3_rk)) cycle
            call series_at_s(params, res%n_active, x(j) + H, r, dr)
            call lemniscate_point_s(x(j) + H, r, dr, res%eps, rho_p, z_p, drb, dzb)
            call series_at_s(params, res%n_active, x(j) - H, r, dr)
            call lemniscate_point_s(x(j) - H, r, dr, res%eps, rho_m, z_m, drb, dzb)
            secant = (rho_p - rho_m) / (z_p - z_m)
            if (abs(drho_dz(j) - secant) > 1.0e-5_rk * (1.0_rk + abs(drho_dz(j)))) ok_fine = .false.
            call series_at_s(params, res%n_active, x(j), r, dr)
            call lemniscate_point_s(x(j), r, dr, res%eps, rb, zb, drb, dzb)
            if (.not. (dzb / res%c > 5.0e-2_rk)) cycle
            secant = (rho(j - 1_ik) - rho(j + 1_ik)) / (z(j - 1_ik) - z(j + 1_ik))
            if (abs(drho_dz(j) - secant) > 1.0e-3_rk * (1.0_rk + abs(drho_dz(j)))) ok = .false.
        end do
        call assert_true(ok, 'C6 drho_dz vs 3-point secant ' // name)
        call assert_true(ok_fine, 'C6 drho_dz vs fine secant near folds ' // name)
    end subroutine test_profile_derivative

    !> C7: the two-lobe node set puts both tips on exact nodes.
    subroutine test_two_lobes()
        type(resolve_t) :: res
        integer(kind = ik) :: status
        real(kind = rk) :: z(100), rho(100), drho_dz(100), z4(4), rho4(4), d4(4)
        real(kind = rk) :: zp(101), rhop(101), dp(101)
        call test_resolve_s(SEPARATED, res, status)
        call profile_s(SEPARATED, res, 100_ik, z, rho, drho_dz)
        call assert_true(is_zero_f(rho(50)) .and. is_zero_f(rho(51)), 'C7 SEPARATED tips rho = 0')
        call assert_abs_close(z(50), sqrt(0.05_rk) / res%c, 1.0e-15_rk, 'C7 north tip z')
        call assert_abs_close(z(51), -z(50), 1.0e-15_rk, 'C7 south tip z = -north tip')
        call assert_true(rho(49) > 0.0_rk .and. rho(52) > 0.0_rk, 'C7 rho > 0 beside the tips')
        call assert_true(is_zero_f(drho_dz(50)) .and. is_zero_f(drho_dz(51)), 'C7 tips drho_dz = 0')
        call assert_abs_close(z(1), res%r_north, 1.0e-15_rk, 'C7 z(1) = r_north')
        call assert_abs_close(z(100), -res%r_south, 1.0e-15_rk, 'C7 z(100) = -r_south')
        call test_resolve_s(PINCHED, res, status)
        call profile_s(PINCHED, res, 101_ik, zp, rhop, dp)
        call assert_true(is_zero_f(zp(51)) .and. is_zero_f(zp(52)), 'C7 PINCHED tips at z = 0')
        call assert_true(is_zero_f(rhop(51)) .and. is_zero_f(rhop(52)), 'C7 PINCHED tips rho = 0')
        call test_resolve_s(SEPARATED, res, status)
        call profile_s(SEPARATED, res, 4_ik, z4, rho4, d4)
        call assert_true(all_zero_f(rho4), 'C7 n = 4: pole, tip, tip, pole')
        call assert_true(z4(1) > z4(2) .and. z4(2) > 0.0_rk .and. z4(3) < 0.0_rk .and. z4(3) > z4(4), &
                'C7 n = 4 order')
    end subroutine test_two_lobes

    !> C8: bare alpha = -1/sqrt(2) folds for |x| > cos(pi/8); the 129-node
    !! profile has that x as a node, where dz_dx is 0 to rounding (s + p
    !! is 1.1e-16 in double, not 0), and the nodes north of it are folded.
    subroutine test_fold_boundary_node()
        real(kind = rk) :: bare(1)
        type(resolve_t) :: res
        integer(kind = ik) :: status, j
        real(kind = rk) :: x(129), lobe(129), z(129), rho(129), drho_dz(129)
        real(kind = rk) :: r, dr, rb, zb, drb, dzb
        logical :: north_folded
        bare(1) = -1.0_rk / sqrt(2.0_rk)
        call test_resolve_s(bare, res, status)
        call assert_int_eq(status, 0_ik, 'C8 bare -1/sqrt(2) resolves')
        call profile_nodes_s(res, 129_ik, x, lobe)
        call profile_s(bare, res, 129_ik, z, rho, drho_dz)
        call assert_abs_close(x(17), cos(acos(-1.0_rk) / 8.0_rk), 1.0e-15_rk, 'C8 node 17 is cos(pi/8)')
        call series_at_s(bare, 1_ik, x(17), r, dr)
        call lemniscate_point_s(x(17), r, dr, res%eps, rb, zb, drb, dzb)
        call assert_true(abs(dzb) <= 1.0e-14_rk, 'C8 dz_dx = 0 to rounding at cos(pi/8)')
        call assert_true(rho(17) > 0.0_rk .and. ieee_is_finite(drho_dz(17)), 'C8 node regular, finite')
        north_folded = .true.
        do j = 2_ik, 16_ik
            call series_at_s(bare, 1_ik, x(j), r, dr)
            call lemniscate_point_s(x(j), r, dr, res%eps, rb, zb, drb, dzb)
            if (.not. (dzb < 0.0_rk)) north_folded = .false.
        end do
        call assert_true(north_folded, 'C8 dz_dx < 0 north of the fold boundary')
        call assert_true(abs(drho_dz(16)) > 10.0_rk .and. abs(drho_dz(18)) > 10.0_rk, &
                'C8 |drho_dz| > 10 beside the boundary')
    end subroutine test_fold_boundary_node

    !> C9: the odd profile's middle node is the exact equator.
    subroutine test_odd_grid_equator()
        type(resolve_t) :: res
        integer(kind = ik) :: status
        real(kind = rk) :: z(65), rho(65), drho_dz(65)
        call test_resolve_s(GS, res, status)
        call profile_s(GS, res, 65_ik, z, rho, drho_dz)
        call assert_true(is_zero_f(z(33) + res%z_shift), 'C9 z(33) = -z_shift exactly')
        call assert_abs_close(rho(33), sqrt(res%r_zero**2 - res%eps) / res%c, 1.0e-15_rk, &
                'C9 rho(33) = sqrt(R0^2 - s)/c')
        call test_resolve_s(ASYM, res, status)
        call profile_s(ASYM, res, 65_ik, z, rho, drho_dz)
        call assert_true(is_zero_f(z(33) + res%z_shift), 'C9 ASYM z(33) = -z_shift exactly')
    end subroutine test_odd_grid_equator

end program cassini_param_core_test
