!> Arithmetic of the Cassini-oval nuclear shape library.
!!
!! Owns the constants of the library and every piece of arithmetic: the
!! Legendre tables and the series R(x), eps from alpha, the lemniscate map
!! with its conditioning rules and derivatives, the node sets, the
!! Gauss-Legendre resolve in t (volume, centre of mass, quadrupole), the scan
!! with the star margin, the neck, the profile and the R(theta) root-find.
module cassini_parameterization_workers_mod

    use precision_utilities_mod, only: ik, ikl, rk
    use mathematical_and_physical_constants_mod, only: PI_C
    use mathematical_utilities_mod, only: compute_legendre_polynomials_s, &
            compute_gauss_legendre_quadrature_s
    use shape_core_mod, only: SHAPE_VALID

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
    ! Library status codes (100+); the shared 0-5 live in shape_core_mod
    !---------------------------------------------------------------------------
    !! Owned here because the resolve raises them and cannot depend on the
    !! front module; the front module re-exports them.
    !> R(x) <= 0 at a node of the fixed check grid.
    integer(kind = ik), parameter, public :: CASSINI_ERROR_R_NONPOSITIVE = 100_ik
    !> R(1)^2 + s <= 0 or R(-1)^2 + s <= 0: a pole has left the axis.
    integer(kind = ik), parameter, public :: CASSINI_ERROR_TOROIDAL = 101_ik
    !> R(theta) outputs and shape only: pinched, or d_min below the margin.
    integer(kind = ik), parameter, public :: CASSINI_ERROR_NOT_STAR_CONVEX = 102_ik
    !> R(theta) outputs, shape and the checked profile: two fragments.
    integer(kind = ik), parameter, public :: CASSINI_ERROR_SEPARATED = 103_ik

    real(kind = rk), parameter :: FOUR_THIRDS_PI = 4.0_rk * PI_C / 3.0_rk

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

    !---------------------------------------------------------------------------
    ! Resolved shape
    !---------------------------------------------------------------------------
    !> One parameter vector, resolved: eps = s, the volume scale c, the raw
    !! centre of mass zm_raw and z_shift = zm_raw / c, the final-frame poles,
    !! the gap (separated only), R at x = 1, -1, 0, the trimmed length and
    !! the connectivity flags. Zero-filled on a failed resolve.
    type, public :: resolve_t
        real(kind = rk) :: eps = 0.0_rk
        real(kind = rk) :: c = 0.0_rk
        real(kind = rk) :: zm_raw = 0.0_rk
        real(kind = rk) :: z_shift = 0.0_rk
        real(kind = rk) :: r_north = 0.0_rk
        real(kind = rk) :: r_south = 0.0_rk
        real(kind = rk) :: gap = 0.0_rk
        real(kind = rk) :: r_plus = 0.0_rk
        real(kind = rk) :: r_minus = 0.0_rk
        real(kind = rk) :: r_zero = 0.0_rk
        integer(kind = ik) :: n_active = 0_ik
        logical :: is_pinched = .false.
        logical :: is_separated = .false.
    end type resolve_t

    public :: build_tables_s
    public :: n_active_f
    public :: series_on_tables_s
    public :: series_at_s
    public :: epsilon_f
    public :: classification_numbers_s
    public :: lemniscate_point_s
    public :: lemniscate_nodes_s
    public :: quadrature_nodes_s
    public :: cosine_nodes_s
    public :: resolve_s
    public :: quadrupole_f

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

    !===========================================================================
    ! NODE SETS
    !===========================================================================

    !> Quadrature nodes for x-integrals over [-1, 1]: CASSINI_N_GL_HALF
    !! Gauss-Legendre nodes in t on [0, 1] per half, mapped to x = +t^2 (north,
    !! first half of the arrays) and x = -t^2 (south, second half), weight 2 t w.
    !! The halves are exact mirrors, so a reflection-symmetric integrand's odd
    !! moment sums to exactly 0 when each half is summed on its own.
    !!
    !! The substitution: the integrands carry |x| and x and the pinch and gap
    !! live at x = 0, where at the exact pinch A.6's integrand goes like
    !! sqrt|x|; in t it is smooth, and near the pinch the branch points of p
    !! move away from the real axis.
    !!
    !! The weights are recomputed from the nodes as 2 / ((1 - xi^2) P_n'(xi)^2)
    !! with `compute_legendre_polynomials_s`: the weights of
    !! `compute_gauss_legendre_quadrature_s` (fortran-foundations 3.0.0) use
    !! P_n' at the iterate before the last Newton step and are off by up to
    !! 1.9e-14 at n = 64 (their sum by 4.8e-14); its nodes are good to an ulp.
    pure subroutine quadrature_nodes_s(x, w)
        real(kind = rk), intent(out) :: x(2 * CASSINI_N_GL_HALF), w(2 * CASSINI_N_GL_HALF)
        real(kind = rk) :: nodes(CASSINI_N_GL_HALF), weights(CASSINI_N_GL_HALF), t
        real(kind = rk) :: p(CASSINI_N_GL_HALF + 1_ik), dp(CASSINI_N_GL_HALF + 1_ik)
        integer(kind = ik) :: j
        call compute_gauss_legendre_quadrature_s(CASSINI_N_GL_HALF, nodes, weights)
        do j = 1_ik, CASSINI_N_GL_HALF
            call compute_legendre_polynomials_s(CASSINI_N_GL_HALF + 1_ik, nodes(j), p, dp)
            weights(j) = 2.0_rk / ((1.0_rk - nodes(j)) * (1.0_rk + nodes(j)) &
                    * dp(CASSINI_N_GL_HALF + 1_ik)**2)
            t = (nodes(j) + 1.0_rk) / 2.0_rk
            x(j) = t * t
            w(j) = t * weights(j)
            x(CASSINI_N_GL_HALF + j) = -x(j)
            w(CASSINI_N_GL_HALF + j) = w(j)
        end do
    end subroutine quadrature_nodes_s

    !> Symmetrized cosine grid of n >= 2 nodes: x_j = cos(pi j / (n - 1)),
    !! j = 0 .. n - 1, built for j < n/2 and mirrored (x_{n-1-j} = -x_j), the
    !! middle node exactly 0 when n is odd. x(1) = 1, x(n) = -1, descending.
    pure subroutine cosine_nodes_s(n, x)
        integer(kind = ik), intent(in) :: n
        real(kind = rk), intent(out) :: x(n)
        integer(kind = ik) :: j, half
        half = n / 2_ik
        do j = 0_ik, half - 1_ik
            x(j + 1_ik) = cos(PI_C * real(j, rk) / real(n - 1_ik, rk))
            x(n - j) = -x(j + 1_ik)
        end do
        if (mod(n, 2_ik) == 1_ik) x(half + 1_ik) = 0.0_rk
    end subroutine cosine_nodes_s

    !===========================================================================
    ! RESOLVE
    !===========================================================================

    !> Resolves one parameter vector: trim, eps, the code-100 check on the
    !! check grid, the code-101 pole test, then A.6/A.7 on the quadrature
    !! nodes (Vbar, zm_raw, c), the final-frame poles, the gap and the flags.
    !!
    !! Pinched is s = R(0)^2 bit for bit (exact at alpha = 1, see `epsilon_f`);
    !! separated is s > R(0)^2. Both are tests on one number.
    !!
    !! @param[in]  params     [alpha, alpha1, ...], at least one element
    !! @param[in]  quad_tab   Legendre table on `quadrature_nodes_s`, order >= n_active - 1
    !! @param[in]  quad_w     The quadrature weights
    !! @param[in]  check_tab  Legendre table on the check grid, same order rule
    !! @param[out] res        The resolved shape; zero-filled on failure
    !! @param[out] status     SHAPE_VALID, CASSINI_ERROR_R_NONPOSITIVE or CASSINI_ERROR_TOROIDAL
    pure subroutine resolve_s(params, quad_tab, quad_w, check_tab, res, status)
        real(kind = rk), intent(in) :: params(:)
        type(tables_t), intent(in) :: quad_tab, check_tab
        real(kind = rk), intent(in) :: quad_w(:)
        type(resolve_t), intent(out) :: res
        integer(kind = ik), intent(out) :: status

        real(kind = rk) :: r_check(check_tab%n_nodes), dr_check(check_tab%n_nodes)
        real(kind = rk) :: eps, r_plus, r_minus, r_zero, pole_plus, pole_minus, equator
        real(kind = rk) :: f0, f1, f2, zm_raw, c
        integer(kind = ik) :: n_active

        status = SHAPE_VALID
        n_active = n_active_f(params)

        call series_on_tables_s(params, n_active, check_tab, r_check, dr_check)
        if (any(.not. (r_check > 0.0_rk))) then
            status = CASSINI_ERROR_R_NONPOSITIVE
            return
        end if

        eps = epsilon_f(params, n_active)
        call classification_numbers_s(params, n_active, eps, r_plus, r_minus, r_zero, &
                pole_plus, pole_minus, equator)
        if (.not. (pole_plus > 0.0_rk) .or. .not. (pole_minus > 0.0_rk)) then
            status = CASSINI_ERROR_TOROIDAL
            return
        end if

        call raw_moments_s(params, n_active, eps, quad_tab, quad_w, .false., f0, f1, f2)
        zm_raw = f1 / f0
        c = (PI_C * f0 / FOUR_THIRDS_PI)**(1.0_rk / 3.0_rk)

        res%eps = eps
        res%c = c
        res%zm_raw = zm_raw
        res%z_shift = zm_raw / c
        res%r_north = (sqrt(pole_plus) - zm_raw) / c
        res%r_south = (sqrt(pole_minus) + zm_raw) / c
        res%r_plus = r_plus
        res%r_minus = r_minus
        res%r_zero = r_zero
        res%n_active = n_active
        res%is_pinched = transfer(r_zero**2, 0_ikl) == transfer(eps, 0_ikl)
        res%is_separated = equator < 0.0_rk
        if (res%is_separated) res%gap = 2.0_rk * sqrt(-equator) / c

    end subroutine resolve_s

    !> q2 = <r^2 P2(cos theta)> in R0 = 1 units = 3/16 c^-5 (Qbar - 4 zm_raw^2 Vbar/pi),
    !! Qbar = int rho_bar^2 (4 z_bar^2 - rho_bar^2) dz_bar (A.8, A.15) on the
    !! quadrature nodes. Pashkevich's Q is 16 q2 / 3.
    pure function quadrupole_f(params, res, quad_tab, quad_w) result(q2)
        real(kind = rk), intent(in) :: params(:)
        type(resolve_t), intent(in) :: res
        type(tables_t), intent(in) :: quad_tab
        real(kind = rk), intent(in) :: quad_w(:)
        real(kind = rk) :: q2
        real(kind = rk) :: f0, f1, qbar
        call raw_moments_s(params, res%n_active, res%eps, quad_tab, quad_w, .true., f0, f1, qbar)
        q2 = 3.0_rk / 16.0_rk * (qbar - 4.0_rk * res%zm_raw**2 * f0) / res%c**5
    end function quadrupole_f

    !> Raw-frame moments on the quadrature nodes, each half summed on its own
    !! (so an odd moment of a reflection-symmetric shape is exactly 0):
    !! f0 = int rho_bar^2 dz_bar = Vbar / pi, f1 = int z_bar rho_bar^2 dz_bar
    !! and, with with_q, f2 = int rho_bar^2 (4 z_bar^2 - rho_bar^2) dz_bar = Qbar.
    !! A.6 is rho_bar^2 dz_bar/dx with the A.1/A.3 derivatives (the printed
    !! bracket is the same integrand); A.7 and A.8 add their z_bar factors.
    pure subroutine raw_moments_s(params, n_active, s, quad_tab, quad_w, with_q, f0, f1, f2)
        real(kind = rk), intent(in) :: params(:)
        integer(kind = ik), intent(in) :: n_active
        real(kind = rk), intent(in) :: s
        type(tables_t), intent(in) :: quad_tab
        real(kind = rk), intent(in) :: quad_w(:)
        logical, intent(in) :: with_q
        real(kind = rk), intent(out) :: f0, f1, f2
        real(kind = rk) :: r(quad_tab%n_nodes), dr(quad_tab%n_nodes)
        real(kind = rk) :: rho(quad_tab%n_nodes), z(quad_tab%n_nodes)
        real(kind = rk) :: drho_dx(quad_tab%n_nodes), dz_dx(quad_tab%n_nodes)
        real(kind = rk) :: dv(quad_tab%n_nodes)
        integer(kind = ik) :: half

        call series_on_tables_s(params, n_active, quad_tab, r, dr)
        call lemniscate_nodes_s(quad_tab%x, r, dr, s, rho, z, drho_dx, dz_dx)
        dv = quad_w * rho**2 * dz_dx
        half = quad_tab%n_nodes / 2_ik
        f0 = sum(dv(1:half)) + sum(dv(half + 1_ik:))
        f1 = sum(z(1:half) * dv(1:half)) + sum(z(half + 1_ik:) * dv(half + 1_ik:))
        f2 = 0.0_rk
        if (with_q) then
            f2 = sum((4.0_rk * z(1:half)**2 - rho(1:half)**2) * dv(1:half)) &
                    + sum((4.0_rk * z(half + 1_ik:)**2 - rho(half + 1_ik:)**2) * dv(half + 1_ik:))
        end if
    end subroutine raw_moments_s

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
