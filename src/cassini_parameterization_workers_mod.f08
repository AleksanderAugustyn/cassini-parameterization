!> Arithmetic of the Cassini-oval nuclear shape library.
!!
!! Owns the constants of the library and every piece of arithmetic: the
!! Legendre tables and the series R(x), eps from alpha, the lemniscate map
!! with its conditioning rules and derivatives, the node sets, the
!! Gauss-Legendre resolve in t (volume, centre of mass, quadrupole), the scan
!! with the star margin, the neck, the profile and the R(theta) root-find.
module cassini_parameterization_workers_mod

    use precision_utilities_mod, only: ik, rk
    use mathematical_utilities_mod, only: compute_legendre_polynomials_s

    implicit none

    private

    !---------------------------------------------------------------------------
    ! Expansion caps
    !---------------------------------------------------------------------------
    !> params = [alpha, alpha1, ..., alpha20]: alpha plus Legendre orders 1 .. 20.
    integer(kind = ik), parameter, public :: CASSINI_MAX_PARAMS = 21_ik
    integer(kind = ik), parameter, public :: CASSINI_MAX_ORDER = CASSINI_MAX_PARAMS - 1_ik

    !---------------------------------------------------------------------------
    ! Grids and floors
    !---------------------------------------------------------------------------
    !> Gauss-Legendre nodes in t per half of [-1, 1] (x = +-t^2); conserves the
    !! bare oval's volume to 1e-13 from eps = 0.95 through the pinch.
    integer(kind = ik), parameter, public :: CASSINI_N_GL_HALF = 64_ik
    !> Cosine-spaced nodes on which code 100 (R(x) <= 0) is judged.
    integer(kind = ik), parameter, public :: CASSINI_N_R_CHECK = 1025_ik
    !> Smallest scan grid (star margin, neck search, R(theta) brackets).
    integer(kind = ik), parameter, public :: CASSINI_N_X_GRID_FLOOR = 65_ik
    !> Smallest profile: the two-lobe set needs a pole and a tip per lobe.
    integer(kind = ik), parameter, public :: CASSINI_N_POINTS_FLOOR = 4_ik

    !---------------------------------------------------------------------------
    ! Tolerances
    !---------------------------------------------------------------------------
    !> Minimum distance from the centre of mass to a tangent line (R0 units)
    !! for the R(theta) outputs. Provisional and inclusive (owner, 2026-10-27).
    real(kind = rk), parameter, public :: STAR_CONVEXITY_MARGIN = 1.0e-3_rk
    !> Two neck candidates whose rho agree to this relative tolerance tie;
    !! the tie goes to the larger x (north).
    real(kind = rk), parameter, public :: NECK_TIE_TOL = 1.0e-12_rk
    !> Bracket width in x at which the root-finds stop.
    real(kind = rk), parameter, public :: X_SOLVE_TOL = 1.0e-14_rk

    !> Length of every string returned by status_message.
    integer(kind = ik), parameter, public :: STATUS_MESSAGE_LEN = 64_ik

    !---------------------------------------------------------------------------
    ! Legendre tables
    !---------------------------------------------------------------------------
    !> Legendre values and derivatives on one node set, orders 0 .. order:
    !! p(n, j) = P_n(x(j)), dp(n, j) = P_n'(x(j)). Parameter-independent; the
    !! cache holds one per grid.
    type, public :: tables_t
        integer(kind = ik) :: order = -1_ik
        integer(kind = ik) :: n_nodes = 0_ik
        real(kind = rk), allocatable :: x(:)
        real(kind = rk), allocatable :: p(:, :)
        real(kind = rk), allocatable :: dp(:, :)
    end type tables_t

    public :: build_tables_s
    public :: n_active_f
    public :: series_on_tables_s
    public :: series_at_s
    public :: epsilon_f
    public :: classification_numbers_s
    public :: lemniscate_point_s
    public :: lemniscate_nodes_s

contains

    !===========================================================================
    ! LEGENDRE SERIES
    !===========================================================================

    !> Tabulates P_n and P_n' (n = 0 .. order) at every node of x.
    !!
    !! The recurrence is the same for every order, so a table of a higher
    !! order carries the bits of a lower one in its leading rows. On an
    !! allocation failure tab is left empty (order = -1, n_nodes = 0).
    !!
    !! @param[in]  x      Nodes in [-1, 1]
    !! @param[in]  order  Highest Legendre order, >= 0
    !! @param[out] tab    The table
    pure subroutine build_tables_s(x, order, tab)
        real(kind = rk), intent(in) :: x(:)
        integer(kind = ik), intent(in) :: order
        type(tables_t), intent(out) :: tab
        integer(kind = ik) :: j, n
        integer :: alloc_stat

        n = size(x, kind = ik)
        allocate(tab%x(n), tab%p(0:order, n), tab%dp(0:order, n), stat = alloc_stat)
        if (alloc_stat /= 0) then
            if (allocated(tab%x)) deallocate(tab%x)
            if (allocated(tab%p)) deallocate(tab%p)
            if (allocated(tab%dp)) deallocate(tab%dp)
            return
        end if
        tab%x = x
        do j = 1_ik, n
            call compute_legendre_polynomials_s(order + 1_ik, x(j), tab%p(:, j), tab%dp(:, j))
        end do
        tab%order = order
        tab%n_nodes = n
    end subroutine build_tables_s

    !> Active length after trimming trailing zeros: max(1, last k with
    !! |params(k)| > 0), capped at size(params) (so 0 for an empty vector).
    pure function n_active_f(params) result(n_active)
        real(kind = rk), intent(in) :: params(:)
        integer(kind = ik) :: n_active
        integer(kind = ik) :: k
        n_active = 1_ik
        do k = size(params, kind = ik), 2_ik, -1_ik
            if (abs(params(k)) > 0.0_rk) then
                n_active = k
                exit
            end if
        end do
        n_active = min(n_active, size(params, kind = ik))
    end function n_active_f

    !> R(x) = 1 + sum_{n=1}^{n_active-1} params(n+1) P_n(x) and R'(x) at every
    !! node of a table of order >= n_active - 1.
    pure subroutine series_on_tables_s(params, n_active, tab, r, dr)
        real(kind = rk), intent(in) :: params(:)
        integer(kind = ik), intent(in) :: n_active
        type(tables_t), intent(in) :: tab
        real(kind = rk), intent(out) :: r(:), dr(:)
        integer(kind = ik) :: j, n
        do j = 1_ik, tab%n_nodes
            r(j) = 1.0_rk
            dr(j) = 0.0_rk
            do n = 1_ik, n_active - 1_ik
                r(j) = r(j) + params(n + 1_ik) * tab%p(n, j)
                dr(j) = dr(j) + params(n + 1_ik) * tab%dp(n, j)
            end do
        end do
    end subroutine series_on_tables_s

    !> R and R' at one x, Legendre values evaluated on the spot (the poles,
    !! the equator and the root-finds).
    pure subroutine series_at_s(params, n_active, x, r, dr)
        real(kind = rk), intent(in) :: params(:)
        integer(kind = ik), intent(in) :: n_active
        real(kind = rk), intent(in) :: x
        real(kind = rk), intent(out) :: r, dr
        real(kind = rk) :: p(CASSINI_MAX_PARAMS), dp(CASSINI_MAX_PARAMS)
        integer(kind = ik) :: n
        r = 1.0_rk
        dr = 0.0_rk
        if (n_active < 2_ik) return
        call compute_legendre_polynomials_s(n_active, x, p(1:n_active), dp(1:n_active))
        do n = 1_ik, n_active - 1_ik
            r = r + params(n + 1_ik) * p(n + 1_ik)
            dr = dr + params(n + 1_ik) * dp(n + 1_ik)
        end do
    end subroutine series_at_s

    !> eps from alpha, eq. 10 collapsed: (alpha - 1)/4 [R(1)^2 + R(-1)^2]
    !! + (alpha + 1)/2 R(0)^2, with R at exactly 1, -1 and 0, evaluated as
    !! R(0)^2 + (alpha - 1) ([R(1)^2 + R(-1)^2]/4 + R(0)^2/2). At alpha = 1 the
    !! correction is multiplied by an exact 0, so eps = R(0)^2 bitwise (the
    !! pinch test is exact); on the bare oval it is 1 + (alpha - 1), which is
    !! alpha exactly whenever alpha - 1 is (the printed grouping rounds
    !! (alpha + 1)/2 and misses alpha = 1.05 by an ulp).
    !! n_active = 0 (an empty vector) is the sphere: eps = 0.
    pure function epsilon_f(params, n_active) result(eps)
        real(kind = rk), intent(in) :: params(:)
        integer(kind = ik), intent(in) :: n_active
        real(kind = rk) :: eps
        real(kind = rk) :: r_plus, r_minus, r_zero, dr
        eps = 0.0_rk
        if (n_active < 1_ik) return
        call series_at_s(params, n_active, 1.0_rk, r_plus, dr)
        call series_at_s(params, n_active, -1.0_rk, r_minus, dr)
        call series_at_s(params, n_active, 0.0_rk, r_zero, dr)
        eps = r_zero**2 + (params(1) - 1.0_rk) &
                * ((r_plus**2 + r_minus**2) / 4.0_rk + r_zero**2 / 2.0_rk)
    end function epsilon_f

    !> R at the poles and the equator, and the three numbers that decide
    !! connectivity: R(1)^2 + s and R(-1)^2 + s (a pole on the axis iff > 0)
    !! and R(0)^2 - s (connected > 0, pinched = 0, separated < 0).
    pure subroutine classification_numbers_s(params, n_active, s, r_plus, r_minus, r_zero, &
            pole_plus, pole_minus, equator)
        real(kind = rk), intent(in) :: params(:)
        integer(kind = ik), intent(in) :: n_active
        real(kind = rk), intent(in) :: s
        real(kind = rk), intent(out) :: r_plus, r_minus, r_zero
        real(kind = rk), intent(out) :: pole_plus, pole_minus, equator
        real(kind = rk) :: dr
        call series_at_s(params, n_active, 1.0_rk, r_plus, dr)
        call series_at_s(params, n_active, -1.0_rk, r_minus, dr)
        call series_at_s(params, n_active, 0.0_rk, r_zero, dr)
        pole_plus = r_plus**2 + s
        pole_minus = r_minus**2 + s
        equator = r_zero**2 - s
    end subroutine classification_numbers_s

    !===========================================================================
    ! LEMNISCATE MAP
    !===========================================================================

    !> Raw coordinates (rho_bar, z_bar) and their TOTAL x-derivatives at one x
    !! for given R, R' = dR/dx and s (eq. 2, A.1-A.5).
    !!
    !! With u = 2x^2 - 1, a = R^2 + s u, b = 2 s |x| sqrt(1 - x^2):
    !!   p = hypot(a, b)                       (never the printed A.5 radicand)
    !!   rho_bar = sqrt(Q-/2), z_bar = sign(x) sqrt(Q+/2), Q+- = p +- (R^2 u + s),
    !! the radicands by `radicands_s` (Q_big direct, Q_small from the product
    !! identity Q+ Q- = 4 R^4 x^2 (1 - x^2)).
    !!
    !! Derivatives. The sign test on a and up picks the form: on the pole side
    !! (a, up of the same sign; a + up is a same-sign sum and Q+ is bounded
    !! away from 0) A.1 and A.3 are printed and
    !!   d rho_bar/dR = 2 |x| sqrt(1 - x^2) R (R^2 + 2us) |z_bar| / (p (a + up))
    !!   d rho_bar/dx = sign(x) (s - p) |z_bar| / (p sqrt(1 - x^2))
    !! on the equator side (opposite signs; a - up is a same-sign sum and Q- is
    !! bounded away from 0) A.2 and A.4 are printed and
    !!   d z_bar/dx = (s + p) rho_bar / (p sqrt(1 - x^2))
    !!   d z_bar/dR = 2 x sqrt(1 - x^2) R (R^2 + 2us) rho_bar / (p (a - up))
    !! Each regular form is its printed form with the vanishing radicand and
    !! bracket rationalized (a +- up = 4 x^2 (1 - x^2) R^2 (R^2 + 2us) / (a -+ up))
    !! and the common factor cancelled, so no form ever divides a vanishing
    !! bracket or radicand: the bracket that appears is always the same-sign
    !! sum, and the rationalized one is folded into the regular forms. The
    !! d rho_bar/dx form above is A.2 through the same product identity; it is
    !! finite at a separated tip (x = 0, where A.2 is 0/0 and the profile has
    !! a kink; sign(0) = +1 gives the north lobe's value).
    !!
    !! Conventions: at the pinch point (p = 0) rho_bar = z_bar = 0 and both
    !! derivatives are 0. At |x| = 1 the partial that diverges (d rho_bar/dx at
    !! a pole on the axis, a true vertical tangent; d z_bar/dx at a toroidal
    !! lift-off) is returned as 0; every caller treats the poles by its own
    !! rule.
    !!
    !! @param[in]  x        Lemniscate coordinate in [-1, 1]
    !! @param[in]  r, dr    R(x) and R'(x)
    !! @param[in]  s        Squared focal distance (eps in R0 = 1 units)
    !! @param[out] rho, z   rho_bar, z_bar (raw frame)
    !! @param[out] drho_dx, dz_dx  Total derivatives along the curve R = R(x)
    pure subroutine lemniscate_point_s(x, r, dr, s, rho, z, drho_dx, dz_dx)
        real(kind = rk), intent(in) :: x, r, dr, s
        real(kind = rk), intent(out) :: rho, z, drho_dx, dz_dx

        real(kind = rk), parameter :: ROOT2 = sqrt(2.0_rk)
        real(kind = rk) :: one_minus_x2, sq, u, r2, a, up, p, q_plus, q_minus
        real(kind = rk) :: bracket, dz_x, dz_r, drho_x, drho_r, sign_x

        rho = 0.0_rk
        z = 0.0_rk
        drho_dx = 0.0_rk
        dz_dx = 0.0_rk

        one_minus_x2 = (1.0_rk - x) * (1.0_rk + x)
        sq = sqrt(max(one_minus_x2, 0.0_rk))
        u = 2.0_rk * x * x - 1.0_rk
        r2 = r * r
        a = r2 + s * u
        p = hypot(a, 2.0_rk * s * abs(x) * sq)
        if (.not. (p > 0.0_rk)) return

        call radicands_s(r2, u, s, p, x, one_minus_x2, q_plus, q_minus)
        rho = sqrt(q_minus / 2.0_rk)
        z = sign(sqrt(q_plus / 2.0_rk), x)
        sign_x = sign(1.0_rk, x)

        dz_x = 0.0_rk
        dz_r = 0.0_rk
        drho_x = 0.0_rk
        drho_r = 0.0_rk
        up = u * p
        if ((a >= 0.0_rk) .eqv. (up >= 0.0_rk)) then
            bracket = a + up
            if (q_plus > 0.0_rk) then
                dz_x = ROOT2 * abs(x) * r2 * (s + p) / (p * sqrt(q_plus))
                dz_r = sign_x * r * bracket / (ROOT2 * p * sqrt(q_plus))
            end if
            if (abs(bracket) > 0.0_rk) then
                drho_r = 2.0_rk * abs(x) * sq * r * (r2 + 2.0_rk * u * s) * abs(z) / (p * bracket)
            end if
            if (sq > 0.0_rk) drho_x = sign_x * (s - p) * abs(z) / (p * sq)
        else
            bracket = a - up
            dz_r = 2.0_rk * x * sq * r * (r2 + 2.0_rk * u * s) * rho / (p * bracket)
            if (sq > 0.0_rk) dz_x = (s + p) * rho / (p * sq)
            drho_x = ROOT2 * x * r2 * (s - p) / (p * sqrt(q_minus))
            drho_r = r * bracket / (ROOT2 * p * sqrt(q_minus))
        end if

        dz_dx = dz_x + dz_r * dr
        drho_dx = drho_x + drho_r * dr

    end subroutine lemniscate_point_s

    !> `lemniscate_point_s` over a node set; every array has the length of x.
    pure subroutine lemniscate_nodes_s(x, r, dr, s, rho, z, drho_dx, dz_dx)
        real(kind = rk), intent(in) :: x(:), r(:), dr(:), s
        real(kind = rk), intent(out) :: rho(:), z(:), drho_dx(:), dz_dx(:)
        integer(kind = ik) :: j
        do j = 1_ik, size(x, kind = ik)
            call lemniscate_point_s(x(j), r(j), dr(j), s, rho(j), z(j), drho_dx(j), dz_dx(j))
        end do
    end subroutine lemniscate_nodes_s

    !> The two radicands without cancellation. Q_big = p + |R^2 u + s| is a
    !! same-sign sum; Q_small = 4 R^4 x^2 (1 - x^2) / Q_big from the product
    !! identity (0 when Q_big = 0, the pinch point only). Q_big goes to Q+ when
    !! R^2 u + s >= 0 and to Q- otherwise. This keeps rho_bar exact at the
    !! poles, the pinch and in the void, and z_bar exact at a connected equator.
    pure subroutine radicands_s(r2, u, s, p, x, one_minus_x2, q_plus, q_minus)
        real(kind = rk), intent(in) :: r2, u, s, p, x, one_minus_x2
        real(kind = rk), intent(out) :: q_plus, q_minus
        real(kind = rk) :: w, q_big, q_small
        w = r2 * u + s
        q_big = p + abs(w)
        q_small = 0.0_rk
        if (q_big > 0.0_rk) q_small = 4.0_rk * r2 * r2 * x * x * one_minus_x2 / q_big
        if (w >= 0.0_rk) then
            q_plus = q_big
            q_minus = q_small
        else
            q_plus = q_small
            q_minus = q_big
        end if
    end subroutine radicands_s

end module cassini_parameterization_workers_mod
