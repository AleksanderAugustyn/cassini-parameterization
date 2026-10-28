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
    !> Gauss-Legendre nodes in t per half of [-1, 1] (x = +-sin^2(pi t/2), see
    !! `quadrature_nodes_s`). Worst relative volume error against panels in x
    !! graded toward both ends, over a distance delta in [1e-12, 1e-1] from
    !! each boundary: near the pinch and the gap 2.0e-14 at 512 (9.0e-13 at
    !! 256); near the toroidal boundary 5.8e-11 for an asymmetric box shape
    !! (4.7e-10 at 256) and 2.7e-9 for the bare oval at alpha = -1 + delta, where
    !! the pole radius and its conditioning both go to zero.
    integer(kind = ik), parameter, public :: CASSINI_N_GL_HALF = 512_ik
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
    !> Bracket width in x at which the golden-section refinement of the star
    !! margin stops: d is flat at its minimum, so x is not resolvable beyond
    !! ~sqrt(epsilon), while d there is already exact to rounding.
    real(kind = rk), parameter :: MARGIN_X_TOL = 1.0e-9_rk

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
    public :: curve_on_tables_s
    public :: star_margin_s
    public :: neck_s
    public :: profile_nodes_s
    public :: profile_s
    public :: theta_on_scan_s
    public :: radius_at_theta_s
    public :: radius_at_thetas_s

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
        integer(kind = ik) :: j
        do j = 1_ik, tab%n_nodes
            call series_at_node_s(params, n_active, tab, j, r(j), dr(j))
        end do
    end subroutine series_on_tables_s

    !> R and R' at node j of a table (the one summation every table user shares).
    pure subroutine series_at_node_s(params, n_active, tab, j, r, dr)
        real(kind = rk), intent(in) :: params(:)
        integer(kind = ik), intent(in) :: n_active, j
        type(tables_t), intent(in) :: tab
        real(kind = rk), intent(out) :: r, dr
        integer(kind = ik) :: n
        r = 1.0_rk
        dr = 0.0_rk
        do n = 1_ik, n_active - 1_ik
            r = r + params(n + 1_ik) * tab%p(n, j)
            dr = dr + params(n + 1_ik) * tab%dp(n, j)
        end do
    end subroutine series_at_node_s

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
    !! Gauss-Legendre nodes in t on [0, 1] per half, mapped to
    !! x = +sin^2(pi t/2) (north, first half of the arrays) and
    !! x = -sin^2(pi t/2) (south, second half), weight (pi/2) sin(pi t) w.
    !! The halves are exact mirrors, so a reflection-symmetric integrand's odd
    !! moment sums to exactly 0 when each half is summed on its own.
    !!
    !! The substitution: p vanishes at the pinch point (x = 0) and at the
    !! toroidal boundary (x = +-1), and near either the integrands have branch
    !! points close to the real axis (at the exact pinch A.6's integrand goes
    !! like sqrt|x|). The map is quadratic at both ends of each half --
    !! x ~ (pi t/2)^2 near 0 and 1 - |x| ~ (pi (1 - t)/2)^2 near 1 -- which
    !! makes the pinch smooth and moves both sets of branch points away from
    !! the real axis. The spec's x = +-t^2 treats only x = 0: near a pole close
    !! to lifting off it misses by 1.2e-6 (5.6e-5 on the bare oval) at 256.
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
            x(j) = sin(PI_C * t / 2.0_rk)**2
            w(j) = PI_C / 4.0_rk * sin(PI_C * t) * weights(j)
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
        integer(kind = ik) :: j
        do j = 1_ik, n
            x(j) = cosine_node_f(n, j)
        end do
    end subroutine cosine_nodes_s

    !> Node j (1-based) of `cosine_nodes_s(n)`.
    pure function cosine_node_f(n, j) result(x)
        integer(kind = ik), intent(in) :: n, j
        real(kind = rk) :: x
        integer(kind = ik) :: k, half
        k = j - 1_ik
        half = n / 2_ik
        if (k < half) then
            x = cos(PI_C * real(k, rk) / real(n - 1_ik, rk))
        else if (mod(n, 2_ik) == 1_ik .and. k == half) then
            x = 0.0_rk
        else
            x = -cos(PI_C * real(n - 1_ik - k, rk) / real(n - 1_ik, rk))
        end if
    end function cosine_node_f

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

    !===========================================================================
    ! SCAN, STAR MARGIN, NECK, PROFILE
    !===========================================================================

    !> The final-frame curve on a table's nodes: rho = rho_bar / c,
    !! z = z_bar / c - z_shift (= (z_bar - zm_raw) / c, written so that an exact
    !! z_bar = 0 lands on -z_shift bit for bit), derivatives / c.
    pure subroutine curve_on_tables_s(params, res, tab, rho, z, drho_dx, dz_dx)
        real(kind = rk), intent(in) :: params(:)
        type(resolve_t), intent(in) :: res
        type(tables_t), intent(in) :: tab
        real(kind = rk), intent(out) :: rho(:), z(:), drho_dx(:), dz_dx(:)
        real(kind = rk) :: r, dr
        integer(kind = ik) :: j
        ! node by node: no work array of the caller's size (Release builds put
        ! automatic arrays on the stack)
        do j = 1_ik, tab%n_nodes
            call series_at_node_s(params, res%n_active, tab, j, r, dr)
            call lemniscate_point_s(tab%x(j), r, dr, res%eps, rho(j), z(j), drho_dx(j), dz_dx(j))
        end do
        rho = rho / res%c
        z = z / res%c - res%z_shift
        drho_dx = drho_dx / res%c
        dz_dx = dz_dx / res%c
    end subroutine curve_on_tables_s

    !> The star-convexity margin and the fold flag on a final-frame scan.
    !!
    !! d = (rho dz_dx - z drho_dx) / hypot(dz_dx, drho_dx), the signed distance
    !! from the origin (the centre of mass) to the tangent line, at every
    !! interior node off the axis (rho > 0: the poles, a separated shape's
    !! x = 0 tip and the pinch point are on the axis and skipped; the pinch
    !! point also has hypot = 0). d > 0 everywhere is monotone theta(x), i.e.
    !! star-convexity.
    !!
    !! Between scan nodes d can dip below its sampled minimum (by 1.7e-2 on
    !! the 65-node floor scan, where the sampled gate passed shapes that are
    !! not star-convex). So every sampled local minimum of d is refined by a
    !! golden-section search for the minimum of d(x) between its two scan
    !! neighbours (`refined_margin_f`, analytic d, no derivatives of d).
    !! d_min is the smallest refined value: over the paper's lattice at
    !! alpha2 = 0, +-0.1 its gate verdict agrees with a 100,001-node scan on
    !! every point at 65 and 129 nodes. d_min is 0 for a separated shape,
    !! where it is meaningless. is_folded: any interior off-axis scan node
    !! with dz_dx <= 0 (rho(z) multi-valued; a polar dimple).
    !!
    !! @param[in] x  The scan's nodes (the table's x), descending
    pure subroutine star_margin_s(params, res, x, rho, z, drho_dx, dz_dx, d_min, is_folded)
        real(kind = rk), intent(in) :: params(:)
        type(resolve_t), intent(in) :: res
        real(kind = rk), intent(in) :: x(:), rho(:), z(:), drho_dx(:), dz_dx(:)
        real(kind = rk), intent(out) :: d_min
        logical, intent(out) :: is_folded
        real(kind = rk) :: d_left, d_here, d_right
        integer(kind = ik) :: j, n
        logical :: first, ok_left, ok_here, ok_right
        d_min = 0.0_rk
        is_folded = .false.
        first = .true.
        n = size(rho, kind = ik)
        do j = 2_ik, n - 1_ik
            if (rho(j) > 0.0_rk .and. .not. (dz_dx(j) > 0.0_rk)) is_folded = .true.
            call node_distance_s(rho(j), z(j), drho_dx(j), dz_dx(j), d_here, ok_here)
            if (.not. ok_here) cycle
            if (first .or. d_here < d_min) d_min = d_here
            first = .false.
        end do
        if (res%is_separated) then
            d_min = 0.0_rk
            return
        end if
        do j = 2_ik, n - 1_ik
            call node_distance_s(rho(j), z(j), drho_dx(j), dz_dx(j), d_here, ok_here)
            if (.not. ok_here) cycle
            ok_left = .false.
            ok_right = .false.
            if (j > 2_ik) call node_distance_s(rho(j - 1_ik), z(j - 1_ik), drho_dx(j - 1_ik), &
                    dz_dx(j - 1_ik), d_left, ok_left)
            if (j < n - 1_ik) call node_distance_s(rho(j + 1_ik), z(j + 1_ik), drho_dx(j + 1_ik), &
                    dz_dx(j + 1_ik), d_right, ok_right)
            if (ok_left) then
                if (d_left < d_here) cycle
            end if
            if (ok_right) then
                if (d_right < d_here) cycle
            end if
            d_min = min(d_min, refined_margin_f(params, res, x(j + 1_ik), x(j - 1_ik)))
        end do
    end subroutine star_margin_s

    !> d at one scan node, and whether it is evaluated (off the axis, hypot > 0).
    pure subroutine node_distance_s(rho, z, drho_dx, dz_dx, d, ok)
        real(kind = rk), intent(in) :: rho, z, drho_dx, dz_dx
        real(kind = rk), intent(out) :: d
        logical, intent(out) :: ok
        real(kind = rk) :: h
        d = 0.0_rk
        h = hypot(dz_dx, drho_dx)
        ok = rho > 0.0_rk .and. h > 0.0_rk
        if (ok) d = (rho * dz_dx - z * drho_dx) / h
    end subroutine node_distance_s

    !> The minimum of d(x) on [lo, hi] by golden-section search to
    !! MARGIN_X_TOL in x (interior points only; the ends are scan nodes the
    !! caller has already sampled).
    pure function refined_margin_f(params, res, lo, hi) result(d_min)
        real(kind = rk), intent(in) :: params(:)
        type(resolve_t), intent(in) :: res
        real(kind = rk), intent(in) :: lo, hi
        real(kind = rk) :: d_min
        real(kind = rk), parameter :: GOLDEN = 0.6180339887498949_rk
        integer(kind = ik), parameter :: MAX_ITER = 100_ik
        real(kind = rk) :: a, b, c, e, fc, fe
        integer(kind = ik) :: iter
        a = lo
        b = hi
        c = b - GOLDEN * (b - a)
        e = a + GOLDEN * (b - a)
        fc = distance_at_f(params, res, c)
        fe = distance_at_f(params, res, e)
        do iter = 1_ik, MAX_ITER
            if (b - a <= MARGIN_X_TOL) exit
            if (fc < fe) then
                b = e
                e = c
                fe = fc
                c = b - GOLDEN * (b - a)
                fc = distance_at_f(params, res, c)
            else
                a = c
                c = e
                fc = fe
                e = a + GOLDEN * (b - a)
                fe = distance_at_f(params, res, e)
            end if
        end do
        d_min = min(fc, fe)
    end function refined_margin_f

    !> d at one x of the final-frame curve; huge where it is undefined (the
    !! pinch point, hypot = 0).
    pure function distance_at_f(params, res, x) result(d)
        real(kind = rk), intent(in) :: params(:)
        type(resolve_t), intent(in) :: res
        real(kind = rk), intent(in) :: x
        real(kind = rk) :: d
        real(kind = rk) :: r, dr, rho, z, drho_dx, dz_dx, h
        call series_at_s(params, res%n_active, x, r, dr)
        call lemniscate_point_s(x, r, dr, res%eps, rho, z, drho_dx, dz_dx)
        rho = rho / res%c
        z = z / res%c - res%z_shift
        drho_dx = drho_dx / res%c
        dz_dx = dz_dx / res%c
        h = hypot(dz_dx, drho_dx)
        d = huge(1.0_rk)
        if (h > 0.0_rk) d = (rho * dz_dx - z * drho_dx) / h
    end function distance_at_f

    !> The neck (spec "Neck") from a final-frame scan rho on scan_tab's nodes.
    !!
    !! Candidates are the interior local minima of rho (rho(k) < rho(k-1) and
    !! rho(k) <= rho(k+1)); rho = 0 at both poles, so each has a local maximum
    !! on either side, and a polar dimple creates none. Each candidate is
    !! refined to the root of d rho_bar/dx between its scan neighbours
    !! (`neck_root_f`), and the smallest refined rho wins; candidates within
    !! NECK_TIE_TOL relative tie and the tie goes to the larger x (the scan
    !! runs north to south, so the first one is kept). Pinched and separated shapes:
    !! found, rho_neck = 0, z_neck = -z_shift (the pinch, or the gap's middle).
    !!
    !! @param[out] x_neck  Optional: the lemniscate coordinate of the neck
    pure subroutine neck_s(params, res, scan_tab, rho, found, z_neck, rho_neck, x_neck)
        real(kind = rk), intent(in) :: params(:)
        type(resolve_t), intent(in) :: res
        type(tables_t), intent(in) :: scan_tab
        real(kind = rk), intent(in) :: rho(:)
        logical, intent(out) :: found
        real(kind = rk), intent(out) :: z_neck, rho_neck
        real(kind = rk), intent(out), optional :: x_neck
        real(kind = rk) :: x, r, dr, rb, zb, drb, dzb
        integer(kind = ik) :: k

        found = .false.
        z_neck = 0.0_rk
        rho_neck = 0.0_rk
        if (present(x_neck)) x_neck = 0.0_rk

        if (res%is_pinched .or. res%is_separated) then
            found = .true.
            z_neck = -res%z_shift
            return
        end if

        ! every candidate is refined before ranking: two necks whose sampled
        ! rho differ by less than the sampling error rank wrongly otherwise
        do k = 2_ik, size(rho, kind = ik) - 1_ik
            if (.not. (rho(k) < rho(k - 1_ik) .and. rho(k) <= rho(k + 1_ik))) cycle
            x = neck_root_f(params, res, scan_tab%x(k + 1_ik), scan_tab%x(k), scan_tab%x(k - 1_ik))
            call series_at_s(params, res%n_active, x, r, dr)
            call lemniscate_point_s(x, r, dr, res%eps, rb, zb, drb, dzb)
            if (found) then
                if (.not. (rb / res%c < rho_neck - NECK_TIE_TOL * rho_neck)) cycle
            end if
            found = .true.
            rho_neck = rb / res%c
            z_neck = zb / res%c - res%z_shift
            if (present(x_neck)) x_neck = x
        end do
    end subroutine neck_s

    !> The minimum of rho_bar near a scan minimum at x_mid, scan neighbours
    !! x_lo < x_mid < x_hi, as a root of g = d rho_bar/dx.
    !!
    !! The half-bracket on the side g(x_mid) points to is closed when g
    !! changes sign across it; g(x_mid) = 0 exactly (x_mid = 0 on a
    !! reflection-symmetric shape) returns x_mid. Otherwise the minimum shares
    !! a scan interval with a maximum (TWO_NECKS: a 1e-5-deep minimum, g < 0
    !! at all three scan nodes): g is sampled on NECK_SUBDIVISIONS equal
    !! sub-intervals of [x_lo, x_hi], every - to + change is closed, and the
    !! smallest rho_bar wins. With no sign change anywhere x_mid stands.
    pure function neck_root_f(params, res, x_lo, x_mid, x_hi) result(x)
        real(kind = rk), intent(in) :: params(:)
        type(resolve_t), intent(in) :: res
        real(kind = rk), intent(in) :: x_lo, x_mid, x_hi
        real(kind = rk) :: x
        integer(kind = ik), parameter :: NECK_SUBDIVISIONS = 16_ik
        real(kind = rk) :: rho_mid, g_mid, rho_end, g_end, xa, ga, xb, gb, xr, rho_r, g_r, best_rho
        integer(kind = ik) :: i

        x = x_mid
        call neck_point_s(params, res, x_mid, rho_mid, g_mid)
        if (.not. (abs(g_mid) > 0.0_rk)) return
        if (g_mid > 0.0_rk) then
            call neck_point_s(params, res, x_lo, rho_end, g_end)
            if (g_end < 0.0_rk) then
                x = illinois_root_f(params, res, x_lo, g_end, x_mid, g_mid)
                return
            end if
        else
            call neck_point_s(params, res, x_hi, rho_end, g_end)
            if (g_end > 0.0_rk) then
                x = illinois_root_f(params, res, x_mid, g_mid, x_hi, g_end)
                return
            end if
        end if

        best_rho = rho_mid
        xa = x_lo
        call neck_point_s(params, res, xa, rho_end, ga)
        do i = 1_ik, NECK_SUBDIVISIONS
            xb = x_lo + (x_hi - x_lo) * real(i, rk) / real(NECK_SUBDIVISIONS, rk)
            call neck_point_s(params, res, xb, rho_end, gb)
            if (ga < 0.0_rk .and. gb > 0.0_rk) then
                xr = illinois_root_f(params, res, xa, ga, xb, gb)
                call neck_point_s(params, res, xr, rho_r, g_r)
                if (rho_r < best_rho) then
                    best_rho = rho_r
                    x = xr
                end if
            end if
            xa = xb
            ga = gb
        end do
    end function neck_root_f

    !> The root of g = d rho_bar/dx in [lo, hi] with g(lo) < 0 < g(hi), by
    !! regula falsi (Illinois: an end that stays twice has its g halved) to
    !! X_SOLVE_TOL in x, with bisection whenever the secant point leaves the
    !! bracket.
    pure function illinois_root_f(params, res, lo_in, g_lo_in, hi_in, g_hi_in) result(x)
        real(kind = rk), intent(in) :: params(:)
        type(resolve_t), intent(in) :: res
        real(kind = rk), intent(in) :: lo_in, g_lo_in, hi_in, g_hi_in
        real(kind = rk) :: x
        integer(kind = ik), parameter :: MAX_ITER = 200_ik
        real(kind = rk) :: lo, hi, g_lo, g_hi, g_new, rho_new
        integer(kind = ik) :: iter, last
        lo = lo_in
        hi = hi_in
        g_lo = g_lo_in
        g_hi = g_hi_in
        x = 0.5_rk * (lo + hi)
        last = 0_ik
        do iter = 1_ik, MAX_ITER
            if (hi - lo <= X_SOLVE_TOL) exit
            x = hi - g_hi * (hi - lo) / (g_hi - g_lo)
            if (.not. (x > lo .and. x < hi)) x = 0.5_rk * (lo + hi)
            call neck_point_s(params, res, x, rho_new, g_new)
            if (.not. (abs(g_new) > 0.0_rk)) exit
            if (g_new > 0.0_rk) then
                hi = x
                g_hi = g_new
                if (last == 1_ik) g_lo = 0.5_rk * g_lo
                last = 1_ik
            else
                lo = x
                g_lo = g_new
                if (last == -1_ik) g_hi = 0.5_rk * g_hi
                last = -1_ik
            end if
        end do
    end function illinois_root_f

    !> rho_bar and g = d rho_bar/dx (total, raw frame) at one x.
    pure subroutine neck_point_s(params, res, x, rho_bar, g)
        real(kind = rk), intent(in) :: params(:)
        type(resolve_t), intent(in) :: res
        real(kind = rk), intent(in) :: x
        real(kind = rk), intent(out) :: rho_bar, g
        real(kind = rk) :: r, dr, zb, dzb
        call series_at_s(params, res%n_active, x, r, dr)
        call lemniscate_point_s(x, r, dr, res%eps, rho_bar, zb, g, dzb)
    end subroutine neck_point_s

    !> The profile's node set, north pole first.
    !!
    !! Connected: `cosine_nodes_s(n)`. Pinched or separated: two lobes. The
    !! north lobe has m1 = ceil(n/2) nodes x = cos(phi), phi uniform on
    !! [0, pi/2], its last node exactly x = 0 (the north tip); the south lobe
    !! the other m2 nodes, mirrored the same way, its first node exactly
    !! x = 0 (the south tip). Both tips are nodes, so the polyline runs along
    !! the axis across the gap. lobe_sign is +1 on north nodes, -1 on south.
    !! n >= CASSINI_N_POINTS_FLOOR.
    pure subroutine profile_nodes_s(res, n, x, lobe_sign)
        type(resolve_t), intent(in) :: res
        integer(kind = ik), intent(in) :: n
        real(kind = rk), intent(out) :: x(n), lobe_sign(n)
        integer(kind = ik) :: j
        do j = 1_ik, n
            call profile_node_s(res, n, j, x(j), lobe_sign(j))
        end do
    end subroutine profile_nodes_s

    !> Node j (1-based) of `profile_nodes_s(res, n)`.
    pure subroutine profile_node_s(res, n, j, x, lobe_sign)
        type(resolve_t), intent(in) :: res
        integer(kind = ik), intent(in) :: n, j
        real(kind = rk), intent(out) :: x, lobe_sign
        integer(kind = ik) :: m1, m2, k
        if (res%is_pinched .or. res%is_separated) then
            m1 = (n + 1_ik) / 2_ik
            m2 = n - m1
            if (j < m1) then
                x = cos(PI_C / 2.0_rk * real(j - 1_ik, rk) / real(max(m1 - 1_ik, 1_ik), rk))
                lobe_sign = 1.0_rk
            else if (j == m1) then
                x = 0.0_rk
                lobe_sign = 1.0_rk
            else if (j == m1 + 1_ik) then
                x = 0.0_rk
                lobe_sign = -1.0_rk
            else
                k = j - m1 - 1_ik
                x = -cos(PI_C / 2.0_rk * real(m2 - 1_ik - k, rk) / real(max(m2 - 1_ik, 1_ik), rk))
                lobe_sign = -1.0_rk
            end if
        else
            x = cosine_node_f(n, j)
            lobe_sign = merge(1.0_rk, -1.0_rk, x >= 0.0_rk)
        end if
    end subroutine profile_node_s

    !> The parametric profile (z, rho, drho/dz) in the final frame on
    !! `profile_nodes_s`, north pole first, z descending except where the
    !! curve folds. drho/dz = drho_dx / dz_dx; at the poles, at the tips (the
    !! x = 0 nodes of a pinched or separated shape, whose z_bar takes its
    !! lobe's sign) and at a node with dz_dx = 0 exactly it is 0, with rho as
    !! computed (0 at poles and tips). The arrays are n long.
    pure subroutine profile_s(params, res, n, z, rho, drho_dz)
        real(kind = rk), intent(in) :: params(:)
        type(resolve_t), intent(in) :: res
        integer(kind = ik), intent(in) :: n
        real(kind = rk), intent(out) :: z(:), rho(:), drho_dz(:)
        real(kind = rk) :: x, lobe_sign, r, dr, rb, zb, drb, dzb
        integer(kind = ik) :: j
        logical :: tip
        ! node by node: no work array of the caller's size (Release builds put
        ! automatic arrays on the stack)
        do j = 1_ik, n
            call profile_node_s(res, n, j, x, lobe_sign)
            call series_at_s(params, res%n_active, x, r, dr)
            call lemniscate_point_s(x, r, dr, res%eps, rb, zb, drb, dzb)
            tip = (res%is_pinched .or. res%is_separated) .and. .not. (abs(x) > 0.0_rk)
            if (tip) zb = lobe_sign * abs(zb)
            rho(j) = rb / res%c
            z(j) = zb / res%c - res%z_shift
            if (j == 1_ik .or. j == n .or. tip .or. .not. (abs(dzb) > 0.0_rk)) then
                drho_dz(j) = 0.0_rk
            else
                drho_dz(j) = drb / dzb
            end if
        end do
    end subroutine profile_s

    !===========================================================================
    ! R(THETA)
    !===========================================================================

    !> theta(x) = atan2(rho, z) and dtheta/dx = (z drho_dx - rho dz_dx) / (rho^2 + z^2)
    !! on a final-frame scan (0 where rho = z = 0, the pinch point).
    pure subroutine theta_on_scan_s(rho, z, drho_dx, dz_dx, theta, dtheta_dx)
        real(kind = rk), intent(in) :: rho(:), z(:), drho_dx(:), dz_dx(:)
        real(kind = rk), intent(out) :: theta(:), dtheta_dx(:)
        real(kind = rk) :: r2
        integer(kind = ik) :: j
        do j = 1_ik, size(rho, kind = ik)
            theta(j) = atan2(rho(j), z(j))
            r2 = rho(j)**2 + z(j)**2
            dtheta_dx(j) = 0.0_rk
            if (r2 > 0.0_rk) dtheta_dx(j) = (z(j) * drho_dx(j) - rho(j) * dz_dx(j)) / r2
        end do
    end subroutine theta_on_scan_s

    !> R and dR/dtheta at one theta, by root-finding theta(x) = theta along
    !! the curve. For a shape that passed the star gate theta(x) decreases
    !! strictly from pi at x = -1 to 0 at x = 1 (folds included).
    !!
    !! theta <= 0 and theta >= pi return the analytic poles with dR/dtheta = 0.
    !! Otherwise the scan nodes whose theta straddle it bracket x (bisection on
    !! the scan index) and a safeguarded Newton iteration on theta(x) - theta,
    !! with the analytic dtheta/dx, runs to X_SOLVE_TOL in x (bisection when a
    !! step leaves the bracket or dtheta/dx is 0; at most 60 iterations, then
    !! the bracket midpoint). Then R = hypot(rho, z), dR/dx = (rho rho_x +
    !! z z_x) / R and dR/dtheta = (dR/dx) / (dtheta/dx); 0 if the root lands on
    !! a pole node, where dtheta/dx carries the pole convention.
    !!
    !! @param[in] x_scan      The scan's nodes, descending
    !! @param[in] theta_scan  theta on the scan (`theta_on_scan_s`), ascending
    pure subroutine radius_at_theta_s(params, res, x_scan, theta_scan, theta, r, dr_dtheta)
        real(kind = rk), intent(in) :: params(:)
        type(resolve_t), intent(in) :: res
        real(kind = rk), intent(in) :: x_scan(:), theta_scan(:), theta
        real(kind = rk), intent(out) :: r, dr_dtheta
        integer(kind = ik), parameter :: MAX_ITER = 60_ik
        real(kind = rk) :: x_lo, x_hi, f_lo, f_hi, x, x_new, f, th, dth, rho, z, drho_dx, dz_dx, dr_dx
        integer(kind = ik) :: lo, hi, mid, iter
        logical :: converged

        dr_dtheta = 0.0_rk
        if (theta <= 0.0_rk) then
            r = res%r_north
            return
        end if
        if (theta >= PI_C) then
            r = res%r_south
            return
        end if

        lo = 1_ik
        hi = size(theta_scan, kind = ik)
        do while (hi - lo > 1_ik)
            mid = (lo + hi) / 2_ik
            if (theta_scan(mid) <= theta) then
                lo = mid
            else
                hi = mid
            end if
        end do
        ! theta decreases in x: f = theta(x) - theta is >= 0 at x_lo and <= 0 at x_hi
        x_lo = x_scan(hi)
        f_lo = theta_scan(hi) - theta
        x_hi = x_scan(lo)
        f_hi = theta_scan(lo) - theta

        x = x_hi
        if (abs(f_hi) > 0.0_rk) x = x_lo + f_lo * (x_hi - x_lo) / (f_lo - f_hi)
        converged = .false.
        do iter = 1_ik, MAX_ITER
            call theta_point_s(params, res, x, th, dth, rho, z, drho_dx, dz_dx)
            f = th - theta
            if (.not. (abs(f) > 0.0_rk)) then
                converged = .true.
                exit
            end if
            if (f > 0.0_rk) then
                x_lo = x
            else
                x_hi = x
            end if
            x_new = 0.5_rk * (x_lo + x_hi)
            if (abs(dth) > 0.0_rk) x_new = x - f / dth
            if (.not. (x_new > x_lo .and. x_new < x_hi)) x_new = 0.5_rk * (x_lo + x_hi)
            if (abs(x_new - x) <= X_SOLVE_TOL .or. x_hi - x_lo <= X_SOLVE_TOL) then
                x = x_new
                converged = .true.
                exit
            end if
            x = x_new
        end do
        if (.not. converged) x = 0.5_rk * (x_lo + x_hi)

        call theta_point_s(params, res, x, th, dth, rho, z, drho_dx, dz_dx)
        r = hypot(rho, z)
        dr_dx = (rho * drho_dx + z * dz_dx) / r
        if (abs(dth) > 0.0_rk) dr_dtheta = dr_dx / dth
    end subroutine radius_at_theta_s

    !> `radius_at_theta_s` over a set of thetas.
    pure subroutine radius_at_thetas_s(params, res, x_scan, theta_scan, thetas, radii, dr_dtheta)
        real(kind = rk), intent(in) :: params(:)
        type(resolve_t), intent(in) :: res
        real(kind = rk), intent(in) :: x_scan(:), theta_scan(:), thetas(:)
        real(kind = rk), intent(out) :: radii(:), dr_dtheta(:)
        integer(kind = ik) :: k
        do k = 1_ik, size(thetas, kind = ik)
            call radius_at_theta_s(params, res, x_scan, theta_scan, thetas(k), radii(k), dr_dtheta(k))
        end do
    end subroutine radius_at_thetas_s

    !> The final-frame curve at one x with theta and dtheta/dx.
    pure subroutine theta_point_s(params, res, x, theta, dtheta_dx, rho, z, drho_dx, dz_dx)
        real(kind = rk), intent(in) :: params(:)
        type(resolve_t), intent(in) :: res
        real(kind = rk), intent(in) :: x
        real(kind = rk), intent(out) :: theta, dtheta_dx, rho, z, drho_dx, dz_dx
        real(kind = rk) :: r, dr, r2
        call series_at_s(params, res%n_active, x, r, dr)
        call lemniscate_point_s(x, r, dr, res%eps, rho, z, drho_dx, dz_dx)
        rho = rho / res%c
        z = z / res%c - res%z_shift
        drho_dx = drho_dx / res%c
        dz_dx = dz_dx / res%c
        theta = atan2(rho, z)
        r2 = rho**2 + z**2
        dtheta_dx = 0.0_rk
        if (r2 > 0.0_rk) dtheta_dx = (z * drho_dx - rho * dz_dx) / r2
    end subroutine theta_point_s

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
