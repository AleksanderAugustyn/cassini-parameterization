!> Independent reference computations for the cassini test suites: the
!! printed formulas evaluated literally, and integrations along paths the
!! library does not take.
module cassini_test_reference_mod

    use precision_utilities_mod, only: ik, rk
    use mathematical_and_physical_constants_mod, only: PI_C
    use mathematical_utilities_mod, only: compute_gauss_legendre_quadrature_s
    use cassini_parameterization_workers_mod, only: resolve_t, tables_t, series_at_s, &
            lemniscate_point_s, quadrature_nodes_s, cosine_nodes_s, build_tables_s, resolve_s, &
            curve_on_tables_s, CASSINI_N_GL_HALF, CASSINI_N_R_CHECK, CASSINI_MAX_ORDER

    implicit none

    private

    public :: reference_epsilon_f
    public :: reference_volume_s
    public :: reference_quadrupole_f
    public :: reference_gl_volume_s
    public :: test_resolve_s
    public :: test_scan_s

    !> Simpson intervals per half of [-1, 1] (10,001 nodes per half, x = 0 shared).
    integer(kind = ik), parameter :: SIMPSON_HALF = 10000_ik

contains

    !> J eq. 6 as printed: the pole brackets from the plain sums of alpha_n and
    !! (-1)^n alpha_n, the equator bracket from (2m - 1)!! / (2^m m!) built by a
    !! loop. params = [alpha, alpha1, ..., alphaN], any length >= 1.
    pure function reference_epsilon_f(params) result(eps)
        real(kind = rk), intent(in) :: params(:)
        real(kind = rk) :: eps
        real(kind = rk) :: alpha, sum_plus, sum_minus, sum_zero, coeff, sign_n
        integer(kind = ik) :: n, m, k, n_max

        alpha = params(1)
        n_max = size(params, kind = ik) - 1_ik
        sum_plus = 1.0_rk
        sum_minus = 1.0_rk
        sign_n = 1.0_rk
        do n = 1_ik, n_max
            sign_n = -sign_n
            sum_plus = sum_plus + params(n + 1_ik)
            sum_minus = sum_minus + sign_n * params(n + 1_ik)
        end do

        sum_zero = 1.0_rk
        do m = 1_ik, n_max / 2_ik
            coeff = 1.0_rk
            do k = 1_ik, m
                ! (2k - 1) / (2 k) builds (2m - 1)!! / (2^m m!)
                coeff = coeff * real(2_ik * k - 1_ik, rk) / real(2_ik * k, rk)
            end do
            if (mod(m, 2_ik) == 1_ik) coeff = -coeff
            sum_zero = sum_zero + coeff * params(2_ik * m + 1_ik)
        end do

        eps = (alpha - 1.0_rk) / 4.0_rk * (sum_plus**2 + sum_minus**2) &
                + (alpha + 1.0_rk) / 2.0_rk * sum_zero**2
    end function reference_epsilon_f

    !> Raw-frame volume and centre of mass by composite Simpson on the curve in
    !! traversal order: each half x = +-t^2 with t uniform on [0, 1] (the
    !! integrand pi rho_bar^2 dz_bar/dx 2t is smooth in t even at the exact
    !! pinch, where it is sqrt|x| in x). Independent of the library's
    !! Gauss-Legendre nodes and weights.
    subroutine reference_volume_s(params, res, vol, zcm)
        real(kind = rk), intent(in) :: params(:)
        type(resolve_t), intent(in) :: res
        real(kind = rk), intent(out) :: vol, zcm
        real(kind = rk) :: f0, f1, f2
        call simpson_moments_s(params, res, .false., f0, f1, f2)
        vol = PI_C * f0
        zcm = f1 / f0
    end subroutine reference_volume_s

    !> Final-frame q2 = 3/4 int (z^2 rho^2 - rho^4/4) dz by the same Simpson
    !! rule, using the resolve's c and zm_raw.
    function reference_quadrupole_f(params, res) result(q2)
        real(kind = rk), intent(in) :: params(:)
        type(resolve_t), intent(in) :: res
        real(kind = rk) :: q2
        real(kind = rk) :: f0, f1, f2
        call simpson_moments_s(params, res, .true., f0, f1, f2)
        q2 = 0.75_rk * f2
    end function reference_quadrupole_f

    !> Simpson sums along the curve: f0 = int rho^2 dz, f1 = int z rho^2 dz
    !! (raw frame); with final = .true., f2 = int (z^2 rho^2 - rho^4/4) dz in
    !! the final frame.
    subroutine simpson_moments_s(params, res, final, f0, f1, f2)
        real(kind = rk), intent(in) :: params(:)
        type(resolve_t), intent(in) :: res
        logical, intent(in) :: final
        real(kind = rk), intent(out) :: f0, f1, f2
        real(kind = rk) :: t, x, r, dr, rho, z, drho_dx, dz_dx, weight, h, jac, zf, rf, side
        integer(kind = ik) :: k, half
        f0 = 0.0_rk
        f1 = 0.0_rk
        f2 = 0.0_rk
        h = 1.0_rk / real(SIMPSON_HALF, rk)
        do half = 1_ik, 2_ik
            side = merge(1.0_rk, -1.0_rk, half == 1_ik)
            do k = 0_ik, SIMPSON_HALF
                t = real(k, rk) * h
                x = side * t * t
                weight = merge(1.0_rk, merge(4.0_rk, 2.0_rk, mod(k, 2_ik) == 1_ik), &
                        k == 0_ik .or. k == SIMPSON_HALF) * h / 3.0_rk
                call series_at_s(params, res%n_active, x, r, dr)
                call lemniscate_point_s(x, r, dr, res%eps, rho, z, drho_dx, dz_dx)
                ! dx = side 2t dt; both halves integrate over increasing x
                jac = 2.0_rk * t * weight
                f0 = f0 + rho**2 * dz_dx * jac
                f1 = f1 + z * rho**2 * dz_dx * jac
                if (final) then
                    zf = (z - res%zm_raw) / res%c
                    rf = rho / res%c
                    f2 = f2 + (zf**2 * rf**2 - rf**4 / 4.0_rk) * dz_dx / res%c * jac
                end if
            end do
        end do
    end subroutine simpson_moments_s

    !> Raw-frame Vbar and zm_raw by Gauss-Legendre with n_half nodes per half:
    !! in t (x = +-t^2) when substituted, plainly in x otherwise.
    subroutine reference_gl_volume_s(params, n_active, s, n_half, substituted, vbar, zm_raw)
        real(kind = rk), intent(in) :: params(:)
        integer(kind = ik), intent(in) :: n_active, n_half
        real(kind = rk), intent(in) :: s
        logical, intent(in) :: substituted
        real(kind = rk), intent(out) :: vbar, zm_raw
        real(kind = rk) :: nodes(n_half), weights(n_half)
        real(kind = rk) :: t, x, w, r, dr, rho, z, drho_dx, dz_dx, f0, f1, side
        integer(kind = ik) :: j, half
        call compute_gauss_legendre_quadrature_s(n_half, nodes, weights)
        f0 = 0.0_rk
        f1 = 0.0_rk
        do half = 1_ik, 2_ik
            side = merge(1.0_rk, -1.0_rk, half == 1_ik)
            do j = 1_ik, n_half
                t = (nodes(j) + 1.0_rk) / 2.0_rk
                if (substituted) then
                    x = side * t * t
                    w = weights(j) / 2.0_rk * 2.0_rk * t
                else
                    x = side * t
                    w = weights(j) / 2.0_rk
                end if
                call series_at_s(params, n_active, x, r, dr)
                call lemniscate_point_s(x, r, dr, s, rho, z, drho_dx, dz_dx)
                f0 = f0 + w * rho**2 * dz_dx
                f1 = f1 + w * z * rho**2 * dz_dx
            end do
        end do
        vbar = PI_C * f0
        zm_raw = f1 / f0
    end subroutine reference_gl_volume_s

    !> Test plumbing: resolves params on freshly built order-20 tables.
    subroutine test_resolve_s(params, res, status)
        real(kind = rk), intent(in) :: params(:)
        type(resolve_t), intent(out) :: res
        integer(kind = ik), intent(out) :: status
        type(tables_t) :: quad_tab, check_tab
        real(kind = rk) :: quad_x(2 * CASSINI_N_GL_HALF), quad_w(2 * CASSINI_N_GL_HALF)
        real(kind = rk) :: check_x(CASSINI_N_R_CHECK)
        call quadrature_nodes_s(quad_x, quad_w)
        call build_tables_s(quad_x, CASSINI_MAX_ORDER, quad_tab)
        call cosine_nodes_s(CASSINI_N_R_CHECK, check_x)
        call build_tables_s(check_x, CASSINI_MAX_ORDER, check_tab)
        call resolve_s(params, quad_tab, quad_w, check_tab, res, status)
    end subroutine test_resolve_s

    !> Test plumbing: the final-frame curve on an n-node cosine scan.
    subroutine test_scan_s(params, res, n, x, rho, z, drho_dx, dz_dx)
        real(kind = rk), intent(in) :: params(:)
        type(resolve_t), intent(in) :: res
        integer(kind = ik), intent(in) :: n
        real(kind = rk), intent(out) :: x(n), rho(n), z(n), drho_dx(n), dz_dx(n)
        type(tables_t) :: tab
        call cosine_nodes_s(n, x)
        call build_tables_s(x, CASSINI_MAX_ORDER, tab)
        call curve_on_tables_s(params, res, tab, rho, z, drho_dx, dz_dx)
    end subroutine test_scan_s

end module cassini_test_reference_mod
