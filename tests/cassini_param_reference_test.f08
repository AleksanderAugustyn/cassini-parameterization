!> Integration paths independent of the library's x-quadrature: volume and
!! centre of mass from R(theta), q2 from R(theta), plain Gauss-Legendre in
!! x, and graded panels through the near-pinch band.
program cassini_param_reference_test
    use precision_utilities_mod, only: ik, rk
    use mathematical_and_physical_constants_mod, only: PI_C
    use cassini_parameterization_workers_mod, only: resolve_t, tables_t, theta_on_scan_s, &
            radius_at_thetas_s, quadrupole_f, quadrature_nodes_s, build_tables_s, &
            CASSINI_N_GL_HALF, CASSINI_MAX_ORDER
    use cassini_test_fixtures_mod, only: SPHERE, BARE_05, BARE_095, OBLATE, GS, SADDLE1, SD_MIN, &
            SADDLE2, ASYM, DIMPLED, FOLD_BOX, TWO_NECKS
    use cassini_test_reference_mod, only: test_resolve_s, test_scan_s, gl_thetas_s, &
            reference_theta_volume_s, reference_theta_q2_f, reference_gl_volume_s, &
            reference_graded_volume_s
    use test_utils_mod, only: assert_true, assert_abs_close, assert_close, test_summary
    implicit none

    real(kind = rk), parameter :: FOUR_THIRDS_PI = 4.0_rk * PI_C / 3.0_rk

    call test_theta_volume(SPHERE, 400_ik, 1.0e-11_rk, 'SPHERE')
    call test_theta_volume(BARE_05, 400_ik, 1.0e-11_rk, 'BARE_05')
    call test_theta_volume(BARE_095, 400_ik, 1.0e-11_rk, 'BARE_095')
    call test_theta_volume(OBLATE, 400_ik, 1.0e-11_rk, 'OBLATE')
    call test_theta_volume(GS, 400_ik, 1.0e-11_rk, 'GS')
    call test_theta_volume(SADDLE1, 400_ik, 1.0e-11_rk, 'SADDLE1')
    call test_theta_volume(SD_MIN, 400_ik, 1.0e-11_rk, 'SD_MIN')
    call test_theta_volume(SADDLE2, 400_ik, 1.0e-11_rk, 'SADDLE2')
    call test_theta_volume(ASYM, 400_ik, 1.0e-11_rk, 'ASYM')
    call test_theta_volume(DIMPLED, 400_ik, 1.0e-11_rk, 'DIMPLED')
    call test_theta_volume(FOLD_BOX, 400_ik, 1.0e-11_rk, 'FOLD_BOX')
    ! 400 thetas resolve TWO_NECKS's neck bump only to 2.4e-11; 800 and 1600 hold 1e-11
    call test_theta_volume(TWO_NECKS, 800_ik, 1.0e-11_rk, 'TWO_NECKS')
    call test_theta_volume([0.999_rk], 4000_ik, 1.0e-10_rk, '[0.999]')
    call test_near_pinch_graded()
    call test_q2(GS, 'GS')
    call test_q2(SADDLE2, 'SADDLE2')
    call test_q2(ASYM, 'ASYM')
    call test_plain_gl(BARE_05, 'BARE_05')
    call test_plain_gl(SADDLE2, 'SADDLE2')
    call test_band()
    call test_summary()

contains

    !> R(theta) at n Gauss-Legendre thetas on a 129-node scan.
    subroutine radii_s(params, res, n, thetas, weights, radii)
        real(kind = rk), intent(in) :: params(:)
        type(resolve_t), intent(out) :: res
        integer(kind = ik), intent(in) :: n
        real(kind = rk), intent(out) :: thetas(n), weights(n), radii(n)
        real(kind = rk) :: x(129), rho(129), z(129), drho_dx(129), dz_dx(129), theta(129), dth(129)
        real(kind = rk) :: dr(n)
        integer(kind = ik) :: status
        call test_resolve_s(params, res, status)
        call test_scan_s(params, res, 129_ik, x, rho, z, drho_dx, dz_dx)
        call theta_on_scan_s(rho, z, drho_dx, dz_dx, theta, dth)
        call gl_thetas_s(n, thetas, weights)
        call radius_at_thetas_s(params, res, x, theta, thetas, radii, dr)
    end subroutine radii_s

    !> F1, F2: V = 4 pi/3 and zcm = 0 from R(theta) alone.
    subroutine test_theta_volume(params, n, tol, name)
        real(kind = rk), intent(in) :: params(:), tol
        integer(kind = ik), intent(in) :: n
        character(len = *), intent(in) :: name
        type(resolve_t) :: res
        real(kind = rk) :: thetas(n), weights(n), radii(n), vol, zcm
        call radii_s(params, res, n, thetas, weights, radii)
        call reference_theta_volume_s(thetas, weights, radii, vol, zcm)
        call assert_close(vol / FOUR_THIRDS_PI, 1.0_rk, tol, 'F1 theta-frame volume ' // name)
        call assert_abs_close(zcm, 0.0_rk, tol, 'F1 theta-frame zcm ' // name)
    end subroutine test_theta_volume

    !> F2: at 1 - 1e-9 the resolve's Vbar against graded plain GL in x.
    subroutine test_near_pinch_graded()
        type(resolve_t) :: res
        integer(kind = ik) :: status
        real(kind = rk) :: walk(1), vbar, zm
        walk(1) = 1.0_rk - 1.0e-9_rk
        call test_resolve_s(walk, res, status)
        call reference_graded_volume_s(walk, res%n_active, res%eps, 64_ik, vbar, zm)
        call assert_close(FOUR_THIRDS_PI * res%c**3 / vbar, 1.0_rk, 1.0e-12_rk, 'F2 [1 - 1e-9] Vbar vs graded GL in x')
    end subroutine test_near_pinch_graded

    !> F3: the library's q2 against 3/10 int R^5 P2 sin from R(theta).
    subroutine test_q2(params, name)
        real(kind = rk), intent(in) :: params(:)
        character(len = *), intent(in) :: name
        type(resolve_t) :: res
        type(tables_t) :: quad_tab
        real(kind = rk) :: thetas(400), weights(400), radii(400)
        real(kind = rk) :: quad_x(2 * CASSINI_N_GL_HALF), quad_w(2 * CASSINI_N_GL_HALF)
        call radii_s(params, res, 400_ik, thetas, weights, radii)
        call quadrature_nodes_s(quad_x, quad_w)
        call build_tables_s(quad_x, CASSINI_MAX_ORDER, quad_tab)
        call assert_abs_close(quadrupole_f(params, res, quad_tab, quad_w), &
                reference_theta_q2_f(thetas, weights, radii), 1.0e-10_rk, 'F3 q2 vs theta frame ' // name)
    end subroutine test_q2

    !> F4: the t-substituted quadrature against plain GL-512 in x away from the pinch.
    subroutine test_plain_gl(params, name)
        real(kind = rk), intent(in) :: params(:)
        character(len = *), intent(in) :: name
        type(resolve_t) :: res
        integer(kind = ik) :: status
        real(kind = rk) :: vbar, zm
        call test_resolve_s(params, res, status)
        call reference_gl_volume_s(params, res%n_active, res%eps, 512_ik, .false., vbar, zm)
        call assert_close(FOUR_THIRDS_PI * res%c**3 / vbar, 1.0_rk, 1.0e-13_rk, 'F4 Vbar vs plain GL-512 ' // name)
        call assert_abs_close(res%zm_raw, zm, 1.0e-13_rk, 'F4 zm_raw vs plain GL-512 ' // name)
    end subroutine test_plain_gl

    !> F5: the library's quadrature order across the near-pinch band, both
    !! sides of alpha = 1, bare and asymmetric, against graded panels.
    subroutine test_band()
        type(resolve_t) :: res
        integer(kind = ik) :: status, k
        real(kind = rk), parameter :: PINCH_ASYM_TAIL(5) = [0.0_rk, 0.2_rk, 0.0_rk, 0.025_rk, 0.05_rk]
        real(kind = rk) :: params(5), vbar, zm, worst_v, worst_z, delta
        worst_v = 0.0_rk
        worst_z = 0.0_rk
        do k = 0_ik, 120_ik
            delta = 10.0_rk**(-1.0_rk - 11.0_rk * real(k, rk) / 120.0_rk)
            params = PINCH_ASYM_TAIL
            params(1) = 1.0_rk - delta
            if (mod(k, 2_ik) == 1_ik) params(1) = 1.0_rk + delta
            if (mod(k, 4_ik) >= 2_ik) params(2:) = 0.0_rk
            call test_resolve_s(params, res, status)
            call reference_graded_volume_s(params, res%n_active, res%eps, 64_ik, vbar, zm)
            worst_v = max(worst_v, abs(FOUR_THIRDS_PI * res%c**3 / vbar - 1.0_rk))
            worst_z = max(worst_z, abs(res%zm_raw - zm))
        end do
        write(*, '(A,2ES10.2)') 'F5 near-pinch band worst |dV/V|, |d zm_raw|: ', worst_v, worst_z
        call assert_true(worst_v <= 1.0e-12_rk, 'F5 near-pinch band volume')
        call assert_true(worst_z <= 1.0e-12_rk, 'F5 near-pinch band centre of mass')
    end subroutine test_band

end program cassini_param_reference_test
