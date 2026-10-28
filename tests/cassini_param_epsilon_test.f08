!> Legendre tables, trimming, the collapsed eps against the printed double
!! sum, and the three classification numbers.
program cassini_param_epsilon_test
    use precision_utilities_mod, only: ik, ikl, rk
    use cassini_parameterization_workers_mod, only: tables_t, build_tables_s, n_active_f, &
            series_on_tables_s, series_at_s, epsilon_f, classification_numbers_s
    use cassini_test_fixtures_mod, only: GS, SADDLE2, PINCHED, PINCHED_ASYM, SEPARATED, &
            TOROIDAL, TOROIDAL_ONE, R_NEG
    use cassini_test_reference_mod, only: reference_epsilon_f
    use test_utils_mod, only: assert_true, assert_int_eq, assert_abs_close, assert_bits_eq, &
            is_zero_f, test_summary
    implicit none

    integer(kind = ikl) :: lcg_state

    lcg_state = 20261027_ikl

    call test_bare_oval()
    call test_collapsed_equals_printed()
    call test_alpha_one_is_pinch_bits()
    call test_polar_line()
    call test_trimming()
    call test_tables()
    call test_classification()
    call test_summary()

contains

    !> Park-Miller minimal standard generator: uniform in [lo, hi].
    function uniform_f(lo, hi) result(u)
        real(kind = rk), intent(in) :: lo, hi
        real(kind = rk) :: u
        lcg_state = mod(lcg_state * 48271_ikl, 2147483647_ikl)
        u = lo + (hi - lo) * real(lcg_state, rk) / 2147483647.0_rk
    end function uniform_f

    !> P1: the bare oval has eps = alpha exactly.
    subroutine test_bare_oval()
        real(kind = rk), parameter :: alphas(6) = [-0.25_rk, 0.0_rk, 0.5_rk, 0.95_rk, 1.0_rk, 1.05_rk]
        real(kind = rk) :: bare(1)
        integer(kind = ik) :: k
        do k = 1_ik, size(alphas, kind = ik)
            bare(1) = alphas(k)
            call assert_bits_eq(epsilon_f(bare, 1_ik), alphas(k), 'P1 bare oval eps = alpha')
        end do
    end subroutine test_bare_oval

    !> P2: the collapsed eq. 10 against the printed J eq. 6 double sum.
    subroutine test_collapsed_equals_printed()
        real(kind = rk) :: params(21), eps, ref
        integer(kind = ik) :: trial, n, k
        real(kind = rk) :: worst
        worst = 0.0_rk
        do trial = 1_ik, 200_ik
            n = 2_ik + mod(trial, 20_ik)
            params(1) = uniform_f(-0.5_rk, 1.2_rk)
            do k = 2_ik, n
                params(k) = uniform_f(-0.3_rk, 0.3_rk)
            end do
            eps = epsilon_f(params(1:n), n_active_f(params(1:n)))
            ref = reference_epsilon_f(params(1:n))
            worst = max(worst, abs(eps - ref) / (1.0_rk + abs(eps)))
        end do
        call assert_true(worst <= 1.0e-15_rk, 'P2 collapsed == printed on 200 random vectors')
        write(*, '(A,ES10.3)') 'P2 worst relative difference: ', worst
    end subroutine test_collapsed_equals_printed

    !> P3: at alpha = 1, eps = R(0)^2 bitwise for any alpha_n.
    subroutine test_alpha_one_is_pinch_bits()
        real(kind = rk) :: params(21), r_plus, r_minus, r_zero, pole_plus, pole_minus, equator, eps
        integer(kind = ik) :: trial, n, k
        logical :: all_bits
        n = size(PINCHED_ASYM, kind = ik)
        eps = epsilon_f(PINCHED_ASYM, n)
        call classification_numbers_s(PINCHED_ASYM, n, eps, r_plus, r_minus, r_zero, &
                pole_plus, pole_minus, equator)
        call assert_bits_eq(eps, r_zero**2, 'P3 PINCHED_ASYM eps = R(0)^2 bits')
        call assert_true(is_zero_f(equator), 'P3 PINCHED_ASYM equator number = 0')
        all_bits = .true.
        do trial = 1_ik, 20_ik
            n = 2_ik + mod(7_ik * trial, 20_ik)
            params(1) = 1.0_rk
            do k = 2_ik, n
                params(k) = uniform_f(-0.3_rk, 0.3_rk)
            end do
            eps = epsilon_f(params(1:n), n)
            call classification_numbers_s(params(1:n), n, eps, r_plus, r_minus, r_zero, &
                    pole_plus, pole_minus, equator)
            if (.not. is_zero_f(equator)) all_bits = .false.
        end do
        call assert_true(all_bits, 'P3 alpha = 1: R(0)^2 - eps = 0 exactly on 20 random vectors')
    end subroutine test_alpha_one_is_pinch_bits

    !> P4: eps = 0 on the polar line alpha*(alpha_n).
    subroutine test_polar_line()
        real(kind = rk) :: params(5), r_plus, r_minus, r_zero, pole_plus, pole_minus, equator
        real(kind = rk) :: sum_poles, alpha_star
        params = SADDLE2
        call classification_numbers_s(params, 5_ik, 0.0_rk, r_plus, r_minus, r_zero, &
                pole_plus, pole_minus, equator)
        sum_poles = r_plus**2 + r_minus**2
        alpha_star = (sum_poles - 2.0_rk * r_zero**2) / (sum_poles + 2.0_rk * r_zero**2)
        call assert_abs_close(alpha_star, 0.052613_rk, 1.0e-5_rk, 'P4 SADDLE2 alpha*')
        params(1) = alpha_star
        call assert_abs_close(epsilon_f(params, 5_ik), 0.0_rk, 1.0e-15_rk, 'P4 eps(alpha*) = 0')
    end subroutine test_polar_line

    !> P5: trailing zeros are trimmed; a trimmed vector gives the same bits.
    subroutine test_trimming()
        call assert_int_eq(n_active_f([0.5_rk, 0.0_rk, 0.0_rk]), 1_ik, 'P5 [0.5, 0, 0] -> 1')
        call assert_int_eq(n_active_f([0.5_rk, 0.0_rk, 0.1_rk, 0.0_rk]), 3_ik, 'P5 [0.5, 0, 0.1, 0] -> 3')
        call assert_int_eq(n_active_f([0.0_rk]), 1_ik, 'P5 [0] -> 1')
        call assert_bits_eq(epsilon_f([0.5_rk, 0.1_rk, 0.0_rk], n_active_f([0.5_rk, 0.1_rk, 0.0_rk])), &
                epsilon_f([0.5_rk, 0.1_rk], n_active_f([0.5_rk, 0.1_rk])), 'P5 padded eps bits')
    end subroutine test_trimming

    !> P6: tables on seven nodes with the exact poles and equator.
    subroutine test_tables()
        type(tables_t) :: tab
        real(kind = rk) :: x(7), r(7), dr(7)
        integer(kind = ik) :: j
        logical :: ok
        x = [1.0_rk, sqrt(3.0_rk) / 2.0_rk, 0.5_rk, 0.0_rk, -0.5_rk, -sqrt(3.0_rk) / 2.0_rk, -1.0_rk]
        call build_tables_s(x, 20_ik, tab)
        call assert_int_eq(tab%order, 20_ik, 'P6 order')
        call assert_int_eq(tab%n_nodes, 7_ik, 'P6 n_nodes')
        ok = .true.
        do j = 1_ik, 7_ik
            if (abs(tab%p(0, j) - 1.0_rk) > 0.0_rk) ok = .false.
            if (abs(tab%p(1, j) - x(j)) > 0.0_rk) ok = .false.
        end do
        call assert_true(ok, 'P6 p0 = 1, p1 = x')
        call assert_abs_close(tab%p(4, 4), 0.375_rk, 1.0e-15_rk, 'P6 P4(0) = 3/8')
        call assert_abs_close(tab%dp(4, 1), 10.0_rk, 1.0e-13_rk, "P6 P4'(1) = 10")
        call series_on_tables_s(GS, 5_ik, tab, r, dr)
        call assert_abs_close(r(4), 1.0_rk + 0.075_rk * 0.375_rk, 1.0e-15_rk, 'P6 R_GS(0)')
        call series_at_s(GS, 5_ik, 0.0_rk, r(1), dr(1))
        call assert_abs_close(r(1), 1.0_rk + 0.075_rk * 0.375_rk, 1.0e-15_rk, 'P6 series_at_s R_GS(0)')
    end subroutine test_tables

    !> P7: the three classification numbers on the structural fixtures.
    subroutine test_classification()
        real(kind = rk) :: r_plus, r_minus, r_zero, pole_plus, pole_minus, equator, eps
        eps = epsilon_f(PINCHED, 1_ik)
        call classification_numbers_s(PINCHED, 1_ik, eps, r_plus, r_minus, r_zero, &
                pole_plus, pole_minus, equator)
        call assert_true(is_zero_f(equator), 'P7 PINCHED equator = 0')
        eps = epsilon_f(SEPARATED, 1_ik)
        call classification_numbers_s(SEPARATED, 1_ik, eps, r_plus, r_minus, r_zero, &
                pole_plus, pole_minus, equator)
        call assert_true(equator < 0.0_rk, 'P7 SEPARATED equator < 0')
        eps = epsilon_f(TOROIDAL, 1_ik)
        call classification_numbers_s(TOROIDAL, 1_ik, eps, r_plus, r_minus, r_zero, &
                pole_plus, pole_minus, equator)
        call assert_true(pole_plus < 0.0_rk, 'P7 TOROIDAL pole_plus < 0')
        eps = epsilon_f(TOROIDAL_ONE, 2_ik)
        call classification_numbers_s(TOROIDAL_ONE, 2_ik, eps, r_plus, r_minus, r_zero, &
                pole_plus, pole_minus, equator)
        call assert_true(pole_minus < 0.0_rk .and. pole_plus > 0.0_rk, 'P7 TOROIDAL_ONE one-sided')
        eps = epsilon_f(R_NEG, 2_ik)
        call classification_numbers_s(R_NEG, 2_ik, eps, r_plus, r_minus, r_zero, &
                pole_plus, pole_minus, equator)
        call assert_abs_close(r_plus, -0.1_rk, 1.0e-15_rk, 'P7 R_NEG R(1) = -0.1')
    end subroutine test_classification

end program cassini_param_epsilon_test
