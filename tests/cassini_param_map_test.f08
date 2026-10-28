!> The lemniscate map: values, the product identity, exact poles and
!! equator, the hypot case, total derivatives against central differences,
!! branch continuity, and the Codex bracket example.
program cassini_param_map_test
    use precision_utilities_mod, only: ik, rk
    use cassini_parameterization_workers_mod, only: lemniscate_point_s, lemniscate_nodes_s, &
            series_at_s, epsilon_f, n_active_f
    use cassini_test_fixtures_mod, only: SPHERE, BARE_05, SADDLE2, ASYM, DIMPLED, OBLATE, &
            TWO_NECKS, GS, PINCHED, SEPARATED, NEAR_PINCH
    use test_utils_mod, only: assert_true, assert_abs_close, assert_close, is_zero_f, test_summary
    implicit none

    real(kind = rk), parameter :: X_LIST(13) = [-0.999_rk, -0.9_rk, -0.5_rk, -0.3_rk, -0.01_rk, &
            -1.0e-6_rk, 0.0_rk, 1.0e-6_rk, 0.01_rk, 0.3_rk, 0.5_rk, 0.9_rk, 0.999_rk]

    call test_sphere()
    call test_product_identity(BARE_05, 'BARE_05')
    call test_product_identity(SADDLE2, 'SADDLE2')
    call test_product_identity(ASYM, 'ASYM')
    call test_product_identity(DIMPLED, 'DIMPLED')
    call test_product_identity(OBLATE, 'OBLATE')
    call test_poles(BARE_05, 'BARE_05')
    call test_poles(SADDLE2, 'SADDLE2')
    call test_poles(ASYM, 'ASYM')
    call test_poles(DIMPLED, 'DIMPLED')
    call test_poles(OBLATE, 'OBLATE')
    call test_equator()
    call test_near_pinch()
    call test_derivatives(SADDLE2, 'SADDLE2')
    call test_derivatives(ASYM, 'ASYM')
    call test_derivatives(DIMPLED, 'DIMPLED')
    call test_derivatives(OBLATE, 'OBLATE')
    call test_derivatives(TWO_NECKS, 'TWO_NECKS')
    call test_branch_continuity()
    call test_codex_bracket()
    call test_equator_bracket()
    call test_nodes_loop()
    call test_summary()

contains

    !> The map at one x for a parameter vector: R, R' from the series, s = eps.
    subroutine point_s(params, x, rho, z, drho_dx, dz_dx)
        real(kind = rk), intent(in) :: params(:), x
        real(kind = rk), intent(out) :: rho, z, drho_dx, dz_dx
        real(kind = rk) :: r, dr, s
        integer(kind = ik) :: n_active
        n_active = n_active_f(params)
        s = epsilon_f(params, n_active)
        call series_at_s(params, n_active, x, r, dr)
        call lemniscate_point_s(x, r, dr, s, rho, z, drho_dx, dz_dx)
    end subroutine point_s

    !> L1: s = 0 is polar coordinates.
    subroutine test_sphere()
        real(kind = rk) :: rho, z, drho_dx, dz_dx, x
        integer(kind = ik) :: k
        do k = 1_ik, size(X_LIST, kind = ik)
            x = X_LIST(k)
            call point_s(SPHERE, x, rho, z, drho_dx, dz_dx)
            ! (1 - x)(1 + x): 1 - x**2 alone loses 1.6e-13 of drho_dx at x = 0.999
            call assert_abs_close(rho, sqrt((1.0_rk - x) * (1.0_rk + x)), 1.0e-15_rk, 'L1 sphere rho')
            call assert_abs_close(z, x, 1.0e-15_rk, 'L1 sphere z')
            call assert_abs_close(drho_dx, -x / sqrt((1.0_rk - x) * (1.0_rk + x)), 1.0e-13_rk, &
                    'L1 sphere drho_dx')
            call assert_abs_close(dz_dx, 1.0_rk, 1.0e-13_rk, 'L1 sphere dz_dx')
        end do
    end subroutine test_sphere

    !> L2: rho |z| = R^2 |x| sqrt(1 - x^2).
    subroutine test_product_identity(params, name)
        real(kind = rk), intent(in) :: params(:)
        character(len = *), intent(in) :: name
        real(kind = rk) :: rho, z, drho_dx, dz_dx, x, r, dr, lhs, rhs, worst
        integer(kind = ik) :: k
        worst = 0.0_rk
        do k = 1_ik, size(X_LIST, kind = ik)
            x = X_LIST(k)
            call point_s(params, x, rho, z, drho_dx, dz_dx)
            call series_at_s(params, n_active_f(params), x, r, dr)
            lhs = rho * abs(z)
            rhs = r**2 * abs(x) * sqrt(1.0_rk - x**2)
            if (abs(rhs) > 0.0_rk) worst = max(worst, abs(lhs - rhs) / rhs)
            if (.not. (abs(rhs) > 0.0_rk)) worst = max(worst, abs(lhs))
        end do
        call assert_true(worst <= 1.0e-14_rk, 'L2 product identity ' // name)
    end subroutine test_product_identity

    !> L3: x = +-1 put the poles on the axis at z = +-sqrt(R(+-1)^2 + s).
    subroutine test_poles(params, name)
        real(kind = rk), intent(in) :: params(:)
        character(len = *), intent(in) :: name
        real(kind = rk) :: rho, z, drho_dx, dz_dx, r, dr, s
        s = epsilon_f(params, n_active_f(params))
        call point_s(params, 1.0_rk, rho, z, drho_dx, dz_dx)
        call series_at_s(params, n_active_f(params), 1.0_rk, r, dr)
        call assert_true(is_zero_f(rho), 'L3 north pole rho = 0 ' // name)
        call assert_abs_close(z, sqrt(r**2 + s), 1.0e-15_rk, 'L3 north pole z ' // name)
        call point_s(params, -1.0_rk, rho, z, drho_dx, dz_dx)
        call series_at_s(params, n_active_f(params), -1.0_rk, r, dr)
        call assert_true(is_zero_f(rho), 'L3 south pole rho = 0 ' // name)
        call assert_abs_close(z, -sqrt(r**2 + s), 1.0e-15_rk, 'L3 south pole z ' // name)
    end subroutine test_poles

    !> L4: the equator of connected, pinched and separated shapes.
    subroutine test_equator()
        real(kind = rk) :: rho, z, drho_dx, dz_dx, r, dr, s
        integer(kind = ik) :: k
        real(kind = rk), parameter :: CONNECTED(4, 5) = reshape([ &
                0.5_rk, 0.0_rk, 0.0_rk, 0.0_rk, 0.0_rk, &
                0.65_rk, 0.2_rk, 0.0_rk, 0.025_rk, 0.05_rk, &
                0.8_rk, 0.15_rk, -0.05_rk, 0.1_rk, -0.05_rk, &
                -0.25_rk, 0.0_rk, 0.0_rk, 0.0_rk, 0.0_rk], [4, 5], order = [2, 1])
        do k = 1_ik, 4_ik
            s = epsilon_f(CONNECTED(k, :), 5_ik)
            call series_at_s(CONNECTED(k, :), 5_ik, 0.0_rk, r, dr)
            call lemniscate_point_s(0.0_rk, r, dr, s, rho, z, drho_dx, dz_dx)
            call assert_true(is_zero_f(z), 'L4 connected equator z = 0')
            call assert_abs_close(rho, sqrt(r**2 - s), 1.0e-15_rk, 'L4 connected equator rho')
        end do
        call point_s(PINCHED, 0.0_rk, rho, z, drho_dx, dz_dx)
        call assert_true(is_zero_f(rho) .and. is_zero_f(z), 'L4 pinch point rho = z = 0')
        call assert_true(is_zero_f(drho_dx) .and. is_zero_f(dz_dx), 'L4 pinch point derivatives 0')
        s = epsilon_f(SEPARATED, 1_ik)
        call point_s(SEPARATED, 1.0e-300_rk, rho, z, drho_dx, dz_dx)
        call assert_abs_close(z, sqrt(s - 1.0_rk), 1.0e-15_rk, 'L4 separated north tip z')
        call assert_abs_close(rho, 0.0_rk, 1.0e-150_rk, 'L4 separated north tip rho')
        call point_s(SEPARATED, -1.0e-300_rk, rho, z, drho_dx, dz_dx)
        call assert_abs_close(z, -sqrt(s - 1.0_rk), 1.0e-15_rk, 'L4 separated south tip z')
        call assert_abs_close(rho, 0.0_rk, 1.0e-150_rk, 'L4 separated south tip rho')
        call point_s(SEPARATED, 0.0_rk, rho, z, drho_dx, dz_dx)
        call assert_true(is_zero_f(rho), 'L4 separated x = 0 rho = 0')
        call assert_abs_close(z, sqrt(s - 1.0_rk), 1.0e-15_rk, 'L4 separated x = 0 is the north tip')
    end subroutine test_equator

    !> L5: the hypot case; the printed A.5 rounds p to 0 here.
    subroutine test_near_pinch()
        real(kind = rk) :: rho, z, drho_dx, dz_dx
        call point_s(NEAR_PINCH, 0.0_rk, rho, z, drho_dx, dz_dx)
        call assert_close(rho / 1.0e-6_rk, 1.0_rk, 1.0e-4_rk, 'L5 near-pinch equator rho = 1e-6')
    end subroutine test_near_pinch

    !> L6: total derivatives against central differences of the map itself.
    subroutine test_derivatives(params, name)
        real(kind = rk), intent(in) :: params(:)
        character(len = *), intent(in) :: name
        real(kind = rk), parameter :: H = 1.0e-6_rk
        real(kind = rk) :: rho, z, drho_dx, dz_dx, rho_p, z_p, rho_m, z_m, d1, d2, x
        real(kind = rk) :: fd_rho, fd_z
        logical :: ok_rho, ok_z
        integer(kind = ik) :: k
        ok_rho = .true.
        ok_z = .true.
        do k = 1_ik, size(X_LIST, kind = ik)
            x = X_LIST(k)
            call point_s(params, x, rho, z, drho_dx, dz_dx)
            call point_s(params, x + H, rho_p, z_p, d1, d2)
            call point_s(params, x - H, rho_m, z_m, d1, d2)
            fd_rho = (rho_p - rho_m) / (2.0_rk * H)
            fd_z = (z_p - z_m) / (2.0_rk * H)
            if (abs(fd_z - dz_dx) > 1.0e-6_rk * (1.0_rk + abs(dz_dx))) then
                ok_z = .false.
                write(*, '(A,ES12.4,2ES24.16)') '  dz_dx mismatch at x =', x, dz_dx, fd_z
            end if
            if (abs(x) < 0.99_rk) then
                if (abs(fd_rho - drho_dx) > 1.0e-6_rk * (1.0_rk + abs(drho_dx))) then
                    ok_rho = .false.
                    write(*, '(A,ES12.4,2ES24.16)') '  drho_dx mismatch at x =', x, drho_dx, fd_rho
                end if
            end if
        end do
        call assert_true(ok_z, 'L6 dz_dx = central difference ' // name)
        call assert_true(ok_rho, 'L6 drho_dx = central difference ' // name)
    end subroutine test_derivatives

    !> L7: no jump at any branch switch. On ASYM the radicand switches
    !! (R^2 u + s = 0) sit at x = -0.307 and 0.328 and the bracket switches
    !! (u = 0) at +-1/sqrt(2), so the walk spans [-0.8, 0.8] at step 1e-4.
    subroutine test_branch_continuity()
        real(kind = rk) :: prev(4), cur(4), x, worst
        integer(kind = ik) :: k
        worst = 0.0_rk
        call point_s(ASYM, -0.8_rk, prev(1), prev(2), prev(3), prev(4))
        do k = 1_ik, 16000_ik
            x = -0.8_rk + 1.6_rk * real(k, rk) / 16000.0_rk
            call point_s(ASYM, x, cur(1), cur(2), cur(3), cur(4))
            worst = max(worst, maxval(abs(cur - prev)))
            prev = cur
        end do
        call assert_true(worst <= 1.0e-3_rk, 'L7 ASYM continuity across both branch switches')
    end subroutine test_branch_continuity

    !> L8: at R = 1, s = -3/4, u = 2/3 both R^2 + 2us and a - up vanish while
    !! a + up = 1; dz_dx is finite and equals the printed A.1.
    subroutine test_codex_bracket()
        real(kind = rk) :: x, s, u, p, rho, z, drho_dx, dz_dx, a1
        s = -0.75_rk
        x = sqrt((2.0_rk / 3.0_rk + 1.0_rk) / 2.0_rk)
        u = 2.0_rk * x**2 - 1.0_rk
        call lemniscate_point_s(x, 1.0_rk, 0.0_rk, s, rho, z, drho_dx, dz_dx)
        p = sqrt(1.0_rk + 2.0_rk * s * u + s**2)
        a1 = sqrt(2.0_rk) * abs(x) * (s + p) / (p * sqrt(p + u + s))
        call assert_abs_close(dz_dx, a1, 1.0e-13_rk, 'L8 Codex example dz_dx = printed A.1')
        call assert_true(rho > 0.0_rk .and. z > 0.0_rk, 'L8 Codex example point is regular')
    end subroutine test_codex_bracket

    !> L9: the regular equator form: dz_dx at x = 1e-9 equals its value at 0.
    subroutine test_equator_bracket()
        real(kind = rk) :: rho0, z0, drho0, dz0, rho, z, drho_dx, dz_dx
        call point_s(GS, 0.0_rk, rho0, z0, drho0, dz0)
        call point_s(GS, 1.0e-9_rk, rho, z, drho_dx, dz_dx)
        call assert_abs_close(dz_dx, dz0, 1.0e-8_rk, 'L9 GS dz_dx regular through x = 0')
        call assert_true(rho > 0.0_rk, 'L9 GS rho > 0 at x = 1e-9')
    end subroutine test_equator_bracket

    !> The node-set form equals the point form.
    subroutine test_nodes_loop()
        real(kind = rk) :: r(13), dr(13), rho(13), z(13), drho_dx(13), dz_dx(13)
        real(kind = rk) :: rho1, z1, drho1, dz1, s
        integer(kind = ik) :: k
        logical :: ok
        s = epsilon_f(SADDLE2, 5_ik)
        do k = 1_ik, 13_ik
            call series_at_s(SADDLE2, 5_ik, X_LIST(k), r(k), dr(k))
        end do
        call lemniscate_nodes_s(X_LIST, r, dr, s, rho, z, drho_dx, dz_dx)
        ok = .true.
        do k = 1_ik, 13_ik
            call lemniscate_point_s(X_LIST(k), r(k), dr(k), s, rho1, z1, drho1, dz1)
            if (abs(rho1 - rho(k)) > 0.0_rk .or. abs(z1 - z(k)) > 0.0_rk) ok = .false.
            if (abs(drho1 - drho_dx(k)) > 0.0_rk .or. abs(dz1 - dz_dx(k)) > 0.0_rk) ok = .false.
        end do
        call assert_true(ok, 'nodes loop == point form')
    end subroutine test_nodes_loop

end program cassini_param_map_test
