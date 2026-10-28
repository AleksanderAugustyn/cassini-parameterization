!> R(theta) and dR/dtheta by root-finding along the curve: the round trip
!! through the profile, anchors, derivatives against central differences,
!! the analytic poles, a dimple theta, and monotone theta on the scan.
program cassini_param_rtheta_test
    use precision_utilities_mod, only: ik, rk
    use mathematical_and_physical_constants_mod, only: PI_C
    use cassini_parameterization_workers_mod, only: resolve_t, profile_s, theta_on_scan_s, &
            radius_at_theta_s, radius_at_thetas_s, star_margin_s
    use cassini_test_fixtures_mod, only: SPHERE, BARE_05, BARE_095, OBLATE, GS, SADDLE1, SD_MIN, &
            SADDLE2, ASYM, DIMPLED, FOLD_BOX, TWO_NECKS, STAR_FAIL, THETA_ANCHORS
    use cassini_test_reference_mod, only: test_resolve_s, test_scan_s
    use test_utils_mod, only: assert_true, assert_abs_close, assert_bits_eq, is_zero_f, test_summary
    implicit none

    integer(kind = ik), parameter :: N_SCAN = 129_ik

    call test_round_trip(SPHERE, 'SPHERE')
    call test_round_trip(BARE_05, 'BARE_05')
    call test_round_trip(GS, 'GS')
    call test_round_trip(SADDLE2, 'SADDLE2')
    call test_round_trip(ASYM, 'ASYM')
    call test_round_trip(DIMPLED, 'DIMPLED')
    call test_round_trip(FOLD_BOX, 'FOLD_BOX')
    call test_round_trip(TWO_NECKS, 'TWO_NECKS')
    call test_round_trip(OBLATE, 'OBLATE')
    call test_anchors(GS, [1.129948_rk, 0.951027_rk, 1.129948_rk], 'GS')
    call test_anchors(SADDLE2, [1.316662_rk, 0.774451_rk, 1.363059_rk], 'SADDLE2')
    call test_anchors(ASYM, [1.473126_rk, 0.585242_rk, 1.562033_rk], 'ASYM')
    call test_anchors(DIMPLED, [1.285575_rk, 0.687014_rk, 1.285575_rk], 'DIMPLED')
    call test_anchors(FOLD_BOX, [1.082596_rk, 0.873401_rk, 1.006061_rk], 'FOLD_BOX')
    call test_anchors(TWO_NECKS, [0.918086_rk, 1.097093_rk, 0.918086_rk], 'TWO_NECKS')
    call test_bare_095_anchors()
    call test_derivative(SADDLE2, 'SADDLE2')
    call test_derivative(ASYM, 'ASYM')
    call test_derivative(DIMPLED, 'DIMPLED')
    call test_poles()
    call test_dimple_theta()
    call test_monotone(SPHERE, 'SPHERE')
    call test_monotone(BARE_05, 'BARE_05')
    call test_monotone(BARE_095, 'BARE_095')
    call test_monotone(OBLATE, 'OBLATE')
    call test_monotone(GS, 'GS')
    call test_monotone(SADDLE1, 'SADDLE1')
    call test_monotone(SD_MIN, 'SD_MIN')
    call test_monotone(SADDLE2, 'SADDLE2')
    call test_monotone(ASYM, 'ASYM')
    call test_monotone(DIMPLED, 'DIMPLED')
    call test_monotone(FOLD_BOX, 'FOLD_BOX')
    call test_monotone(TWO_NECKS, 'TWO_NECKS')
    call test_monotone(STAR_FAIL, 'STAR_FAIL')
    call test_summary()

contains

    !> The resolve and the scan's x and theta.
    subroutine setup_s(params, res, x, theta, d_min)
        real(kind = rk), intent(in) :: params(:)
        type(resolve_t), intent(out) :: res
        real(kind = rk), intent(out) :: x(N_SCAN), theta(N_SCAN), d_min
        real(kind = rk) :: rho(N_SCAN), z(N_SCAN), drho_dx(N_SCAN), dz_dx(N_SCAN), dtheta_dx(N_SCAN)
        integer(kind = ik) :: status
        logical :: is_folded
        call test_resolve_s(params, res, status)
        call test_scan_s(params, res, N_SCAN, x, rho, z, drho_dx, dz_dx)
        call theta_on_scan_s(rho, z, drho_dx, dz_dx, theta, dtheta_dx)
        call star_margin_s(params, res, x, rho, z, drho_dx, dz_dx, d_min, is_folded)
    end subroutine setup_s

    !> T1: R(atan2(rho, z)) = hypot(z, rho) at 401 profile nodes.
    subroutine test_round_trip(params, name)
        real(kind = rk), intent(in) :: params(:)
        character(len = *), intent(in) :: name
        type(resolve_t) :: res
        real(kind = rk) :: x(N_SCAN), theta(N_SCAN), d_min
        real(kind = rk) :: z(401), rho(401), drho_dz(401), r, dr, worst
        integer(kind = ik) :: j
        call setup_s(params, res, x, theta, d_min)
        call profile_s(params, res, 401_ik, z, rho, drho_dz)
        worst = 0.0_rk
        do j = 1_ik, 401_ik
            call radius_at_theta_s(params, res, x, theta, atan2(rho(j), z(j)), r, dr)
            worst = max(worst, abs(r - hypot(z(j), rho(j))))
        end do
        call assert_true(worst <= 1.0e-13_rk, 'T1 round trip ' // name)
        if (.not. (worst <= 1.0e-13_rk)) write(*, '(A,ES10.3)') '  worst ', worst
    end subroutine test_round_trip

    !> T2: R at pi/8, pi/2, 7 pi/8.
    subroutine test_anchors(params, want, name)
        real(kind = rk), intent(in) :: params(:), want(3)
        character(len = *), intent(in) :: name
        type(resolve_t) :: res
        real(kind = rk) :: x(N_SCAN), theta(N_SCAN), d_min, radii(3), dr(3)
        integer(kind = ik) :: k
        call setup_s(params, res, x, theta, d_min)
        call radius_at_thetas_s(params, res, x, theta, THETA_ANCHORS, radii, dr)
        do k = 1_ik, 3_ik
            call assert_abs_close(radii(k), want(k), 1.0e-5_rk, 'T2 anchor ' // name)
        end do
    end subroutine test_anchors

    subroutine test_bare_095_anchors()
        type(resolve_t) :: res
        real(kind = rk) :: x(N_SCAN), theta(N_SCAN), d_min, r, dr
        call setup_s(BARE_095, res, x, theta, d_min)
        call radius_at_theta_s(BARE_095, res, x, theta, THETA_ANCHORS(1), r, dr)
        call assert_abs_close(r, 1.68623_rk, 1.0e-5_rk, 'T2 BARE_095 R(pi/8)')
        call radius_at_theta_s(BARE_095, res, x, theta, THETA_ANCHORS(2), r, dr)
        call assert_abs_close(r, 0.317251_rk, 1.0e-5_rk, 'T2 BARE_095 R(pi/2)')
    end subroutine test_bare_095_anchors

    !> T3: dR/dtheta against central differences of R(theta), dimple rims included.
    subroutine test_derivative(params, name)
        real(kind = rk), intent(in) :: params(:)
        character(len = *), intent(in) :: name
        real(kind = rk), parameter :: H = 1.0e-6_rk
        type(resolve_t) :: res
        real(kind = rk) :: x(N_SCAN), theta(N_SCAN), d_min, th, r, dr, r_p, r_m, dummy
        integer(kind = ik) :: i
        logical :: ok
        call setup_s(params, res, x, theta, d_min)
        ok = .true.
        do i = 1_ik, 183_ik
            th = PI_C * real(i, rk) / 182.0_rk
            if (i == 182_ik) th = 0.05_rk
            if (i == 183_ik) th = PI_C - 0.05_rk
            call radius_at_theta_s(params, res, x, theta, th, r, dr)
            call radius_at_theta_s(params, res, x, theta, th + H, r_p, dummy)
            call radius_at_theta_s(params, res, x, theta, th - H, r_m, dummy)
            if (abs((r_p - r_m) / (2.0_rk * H) - dr) > 1.0e-6_rk * (1.0_rk + abs(dr))) then
                ok = .false.
                write(*, '(A,3ES16.8)') '  dR/dtheta mismatch', th, dr, (r_p - r_m) / (2.0_rk * H)
            end if
        end do
        call assert_true(ok, 'T3 dR/dtheta vs central difference ' // name)
    end subroutine test_derivative

    !> T4: the poles are analytic.
    subroutine test_poles()
        type(resolve_t) :: res
        real(kind = rk) :: x(N_SCAN), theta(N_SCAN), d_min, r, dr
        call setup_s(SADDLE2, res, x, theta, d_min)
        call radius_at_theta_s(SADDLE2, res, x, theta, 0.0_rk, r, dr)
        call assert_bits_eq(r, res%r_north, 'T4 theta = 0 -> r_north')
        call assert_true(is_zero_f(dr), 'T4 theta = 0 -> dR = 0')
        call radius_at_theta_s(SADDLE2, res, x, theta, PI_C, r, dr)
        call assert_bits_eq(r, res%r_south, 'T4 theta = pi -> r_south')
        call assert_true(is_zero_f(dr), 'T4 theta = pi -> dR = 0')
        call radius_at_theta_s(SADDLE2, res, x, theta, 1.0e-12_rk, r, dr)
        call assert_abs_close(r, res%r_north, 1.0e-10_rk, 'T4 theta = 1e-12 near r_north')
        call radius_at_theta_s(SADDLE2, res, x, theta, PI_C - 1.0e-12_rk, r, dr)
        call assert_abs_close(r, res%r_south, 1.0e-10_rk, 'T4 theta = pi - 1e-12 near r_south')
    end subroutine test_poles

    !> T5: a theta inside DIMPLED's polar dimple has one crossing.
    subroutine test_dimple_theta()
        type(resolve_t) :: res
        real(kind = rk) :: x(N_SCAN), theta(N_SCAN), d_min, r, dr, r_p, r_m, dummy
        call setup_s(DIMPLED, res, x, theta, d_min)
        call radius_at_theta_s(DIMPLED, res, x, theta, 0.05_rk, r, dr)
        call assert_true(r > 0.0_rk .and. r < res%r_north + 0.1_rk, 'T5 DIMPLED R(0.05) is the dimple floor')
        call radius_at_theta_s(DIMPLED, res, x, theta, 0.05_rk + 1.0e-6_rk, r_p, dummy)
        call radius_at_theta_s(DIMPLED, res, x, theta, 0.05_rk - 1.0e-6_rk, r_m, dummy)
        call assert_true(abs(r_p - r_m) < 1.0e-4_rk, 'T5 DIMPLED R continuous at 0.05')
    end subroutine test_dimple_theta

    !> T6: theta strictly increases along the scan whenever d_min > 0.
    subroutine test_monotone(params, name)
        real(kind = rk), intent(in) :: params(:)
        character(len = *), intent(in) :: name
        type(resolve_t) :: res
        real(kind = rk) :: x(N_SCAN), theta(N_SCAN), d_min
        call setup_s(params, res, x, theta, d_min)
        if (d_min > 0.0_rk) then
            call assert_true(all(theta(2:) > theta(:N_SCAN - 1_ik)), 'T6 theta monotone ' // name)
        else
            call assert_true(.not. all(theta(2:) > theta(:N_SCAN - 1_ik)), &
                    'T6 theta not monotone when d_min < 0 ' // name)
        end if
    end subroutine test_monotone

end program cassini_param_rtheta_test
