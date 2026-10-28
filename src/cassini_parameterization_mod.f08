!> Public surface of the Cassini-oval nuclear shape library.
!!
!! The library's front door: `cache_t`, the usage checks, the gates and
!! every cached and one-shot routine. It does no floating-point arithmetic on
!! an output value; the code here checks, calls kernels, applies gates and
!! copies results. Every kernel lives in `cassini_parameterization_workers_mod`
!! (trailing zeros are trimmed there, by the resolve, after these checks).
!!
!! ## The shape
!!
!! Pashkevich's Cassini ovals (Nucl. Phys. A 169 (1971) 275): in lemniscate
!! coordinates (R, x) the surface is R(x) = 1 + sum_n alpha_n P_n(x), and
!! `params = [alpha, alpha1, ..., alphaN]` with eps from alpha by eq. 10.
!! [alpha] alone is the bare oval, [0] the sphere. Everything is in R0 = 1
!! units with the centre of mass at the origin; the library has no notion of
!! A, Z or N.
!!
!! ## Two tiers
!!
!! 1. **One-shot** -- `compute_*_standalone_s`: builds a local cache
!!    (max_params = size(params); the radius forms on the caller's thetas,
!!    the others on a one-node set; the non-scanning forms on the floor scan
!!    grid, which they never consult), calls the cached routine and lets the
!!    cache go out of scope.
!! 2. **Read-only cache** -- `cache_*_s` on a `cache_t` built once by
!!    `cache_init_s`: the Legendre tables (orders 0 .. max_params - 1) on the
!!    quadrature, check and scan grids, the quadrature weights, and the
!!    primary theta set. Nothing derived from a shape parameter is stored.
!!
!! One-shot and cached calls return the same status and bitwise-identical
!! outputs: there is one pipeline, and the one-shot forms are wrappers.
!!
!! ## Thread and lifetime rules (normative)
!!
!! - A cache is immutable after `cache_init_s`; every compute takes it
!!   `intent(in)` and is `pure`. Any number of threads may compute on one
!!   cache concurrently.
!! - `cache_init_s` and `cache_free_s` on a cache must not race with any
!!   other call on it.
!! - No hidden state: the library holds no mutable module-level state.
!! - A `cache_t` has allocatable components only: assignment is a deep copy
!!   and scope exit releases it. `cache_free_s` is available, not mandatory.
!!
!! ## Checks, in order (normative)
!!
!! The first failing check sets the status. Usage codes depend on sizes,
!! handles and grids only, never on parameter values.
!!
!!   cache_init_s   max_params > 21 -> 1; < 1 -> 5; n_x_grid below
!!                  CASSINI_N_X_GRID_FLOOR, an empty theta set, a theta
!!                  outside [0, pi] or a failed allocation -> 3
!!   cached call    uninitialized cache -> 2; size(params) outside
!!                  1 .. max_params -> 4; output extent -> 3; at-thetas only:
!!                  theta set -> 3
!!   one-shot call  empty params -> 4; size(params) > 21 -> 1; output extent,
!!                  scan grid and theta set -> 3
!!
!! Then the value codes: 100 (R(x) <= 0 on the check grid), 101 (a pole off
!! the axis), then the gates each output applies -- 103 (separated), then
!! 102 (pinched, or d_min below STAR_CONVEXITY_MARGIN).
!!
!!   radius_grid, radius_and_derivative, at-thetas, shape   100 101 103 102
!!   rho_profile                                            100 101 103
!!   rho_profile_unchecked, neck, geometry, quadrupole      100 101
!!
!! On any nonzero status every output argument is zero-filled over its
!! entire extent, `found` is .false., and a geometry record is zero-filled
!! with every flag .false.
module cassini_parameterization_mod

    use precision_utilities_mod, only: ik, rk
    use mathematical_and_physical_constants_mod, only: PI_C
    use shape_core_mod, only: SHAPE_VALID, SHAPE_ERROR_TOO_MANY_PARAMS, &
            SHAPE_ERROR_CACHE_NOT_INITIALIZED, SHAPE_ERROR_INVALID_GRID, &
            SHAPE_ERROR_WRONG_PARAM_COUNT, SHAPE_ERROR_INVALID_INIT, &
            SHAPE_MAX_PARAMS
    use cassini_parameterization_workers_mod, only: CASSINI_MAX_PARAMS, CASSINI_MAX_ORDER, &
            CASSINI_N_GL_HALF, CASSINI_N_R_CHECK, CASSINI_N_X_GRID_FLOOR, &
            CASSINI_N_POINTS_FLOOR, STAR_CONVEXITY_MARGIN, NECK_TIE_TOL, X_SOLVE_TOL, &
            STATUS_MESSAGE_LEN, CASSINI_ERROR_R_NONPOSITIVE, CASSINI_ERROR_TOROIDAL, &
            CASSINI_ERROR_NOT_STAR_CONVEX, CASSINI_ERROR_SEPARATED, tables_t, resolve_t, &
            build_tables_s, quadrature_nodes_s, cosine_nodes_s, resolve_s, quadrupole_f, &
            curve_on_tables_s, star_margin_s, neck_s, profile_s, theta_on_scan_s, &
            radius_at_thetas_s, n_active_f, epsilon_f

    implicit none

    private

    !---------------------------------------------------------------------------
    ! Read-only cache
    !---------------------------------------------------------------------------
    public :: cache_init_s
    public :: cache_free_s
    public :: cache_max_params_f
    public :: cache_n_x_grid_f
    public :: cache_n_thetas_f
    public :: cache_is_initialized_f

    !---------------------------------------------------------------------------
    ! Cached tier
    !---------------------------------------------------------------------------
    public :: cache_radius_grid_s
    public :: cache_radius_and_derivative_s
    public :: cache_radius_and_derivative_at_thetas_s
    public :: cache_shape_s
    public :: cache_rho_profile_s
    public :: cache_rho_profile_unchecked_s
    public :: cache_neck_s
    public :: cache_geometry_s
    public :: cache_quadrupole_s

    !---------------------------------------------------------------------------
    ! One-shot tier
    !---------------------------------------------------------------------------
    public :: compute_radius_grid_standalone_s
    public :: compute_radius_and_derivative_standalone_s
    public :: compute_shape_standalone_s
    public :: compute_rho_profile_standalone_s
    public :: compute_rho_profile_unchecked_standalone_s
    public :: compute_neck_standalone_s
    public :: compute_geometry_standalone_s
    public :: compute_quadrupole_standalone_s
    public :: compute_epsilon_f

    public :: status_message

    !---------------------------------------------------------------------------
    ! Constants, owned by the workers module and re-exported so the public
    ! surface is self-contained
    !---------------------------------------------------------------------------
    public :: CASSINI_MAX_PARAMS
    public :: CASSINI_MAX_ORDER
    public :: CASSINI_N_GL_HALF
    public :: CASSINI_N_R_CHECK
    public :: CASSINI_N_X_GRID_FLOOR
    public :: CASSINI_N_POINTS_FLOOR
    public :: STAR_CONVEXITY_MARGIN
    public :: NECK_TIE_TOL
    public :: X_SOLVE_TOL
    public :: STATUS_MESSAGE_LEN
    public :: SHAPE_MAX_PARAMS

    !> Largest params length any call accepts: min(SHAPE_MAX_PARAMS, CASSINI_MAX_PARAMS).
    integer(kind = ik), parameter, public :: CASSINI_PARAM_LIMIT = &
            min(SHAPE_MAX_PARAMS, CASSINI_MAX_PARAMS)

    !---------------------------------------------------------------------------
    ! Status codes
    !---------------------------------------------------------------------------
    !! Two disjoint ranges share one integer space. 0-99 are the shared
    !! shape_core codes, re-exported so a caller never has to use
    !! shape_core_mod directly; 100+ are cassini-specific. CASSINI_VALID is an
    !! alias of SHAPE_VALID, so success compares equal across libraries.
    public :: SHAPE_VALID
    public :: SHAPE_ERROR_TOO_MANY_PARAMS
    public :: SHAPE_ERROR_CACHE_NOT_INITIALIZED
    public :: SHAPE_ERROR_INVALID_GRID
    public :: SHAPE_ERROR_WRONG_PARAM_COUNT
    public :: SHAPE_ERROR_INVALID_INIT

    public :: CASSINI_ERROR_R_NONPOSITIVE
    public :: CASSINI_ERROR_TOROIDAL
    public :: CASSINI_ERROR_NOT_STAR_CONVEX
    public :: CASSINI_ERROR_SEPARATED

    integer(kind = ik), parameter, public :: CASSINI_VALID = SHAPE_VALID

    !> The one-node theta set of the theta-less one-shot forms.
    real(kind = rk), parameter :: ONE_THETA(1) = [0.0_rk]

    !> Read-only cache: Legendre tables (orders 0 .. max_params - 1) on the
    !! quadrature, check and scan grids, the quadrature weights and the
    !! primary theta set.
    !!
    !! Opaque: the components are private and may grow. Immutable after
    !! `cache_init_s`, shared freely across threads.
    type, public :: cache_t
        private
        logical :: initialized = .false.
        integer(kind = ik) :: max_params = 0_ik
        integer(kind = ik) :: n_x_grid = 0_ik
        real(kind = rk), allocatable :: thetas(:)
        real(kind = rk), allocatable :: quad_w(:)
        type(tables_t) :: quad_tab
        type(tables_t) :: check_tab
        type(tables_t) :: scan_tab
    end type cache_t

    !> The resolved shape as one record, in R0 = 1 units, centre of mass at
    !! the origin: eps (= s), the volume scale c, z_shift (the raw centre of
    !! mass over c), the pole distances, the gap (separated shapes only), the
    !! star margin d_min (0 for separated shapes) and the three flags.
    !! Succeeds for every shape that resolves, pinched and separated included.
    type, public :: cassini_geometry_t
        real(kind = rk) :: eps = 0.0_rk
        real(kind = rk) :: c = 0.0_rk
        real(kind = rk) :: z_shift = 0.0_rk
        real(kind = rk) :: r_north = 0.0_rk
        real(kind = rk) :: r_south = 0.0_rk
        real(kind = rk) :: gap = 0.0_rk
        real(kind = rk) :: d_min = 0.0_rk
        logical :: is_separated = .false.
        logical :: is_pinched = .false.
        logical :: is_folded = .false.
    end type cassini_geometry_t

contains

    !===========================================================================
    ! STATUS
    !===========================================================================

    !> Fixed diagnostic string for a status code, shared or cassini-specific.
    !!
    !! @param[in] code  Any integer; unrecognised values get a fallback message
    !! @return          Blank-padded message, never empty
    pure function status_message(code) result(msg)
        integer(kind = ik), intent(in) :: code
        character(len = STATUS_MESSAGE_LEN) :: msg
        select case (code)
        case (SHAPE_VALID);                        msg = 'valid'
        case (SHAPE_ERROR_TOO_MANY_PARAMS);        msg = 'too many parameters'
        case (SHAPE_ERROR_CACHE_NOT_INITIALIZED);  msg = 'cache not initialized'
        case (SHAPE_ERROR_INVALID_GRID);           msg = 'invalid grid: thetas, scan or profile size, or buffer extent'
        case (SHAPE_ERROR_WRONG_PARAM_COUNT);      msg = 'params length is 0 or exceeds max_params'
        case (SHAPE_ERROR_INVALID_INIT);           msg = 'invalid init arguments'
        case (CASSINI_ERROR_R_NONPOSITIVE);        msg = 'R(x) <= 0 on the check grid'
        case (CASSINI_ERROR_TOROIDAL);             msg = 'toroidal: a pole has left the symmetry axis'
        case (CASSINI_ERROR_NOT_STAR_CONVEX);      msg = 'shape not star-convex from the centre of mass'
        case (CASSINI_ERROR_SEPARATED);            msg = 'shape separated into two fragments'
        case default;                              msg = 'unknown status code'
        end select
    end function status_message

    !===========================================================================
    ! CACHE LIFECYCLE
    !===========================================================================

    !> Builds a read-only cache.
    !!
    !! Calling this on a live cache releases the old contents through
    !! `intent(out)`.
    !!
    !! @param[out] cache       Ready on success; uninitialized otherwise
    !! @param[in]  max_params  1 .. 21; sets the Legendre table order max_params - 1
    !! @param[in]  n_x_grid    Scan nodes, >= CASSINI_N_X_GRID_FLOOR
    !! @param[in]  thetas      Primary polar angles in [0, pi], at least one, any order
    !! @param[out] status      0, 1 (> 21), 5 (< 1) or 3 (grid, thetas, allocation)
    pure subroutine cache_init_s(cache, max_params, n_x_grid, thetas, status)

        type(cache_t), intent(out) :: cache
        integer(kind = ik), intent(in) :: max_params, n_x_grid
        real(kind = rk), intent(in) :: thetas(:)
        integer(kind = ik), intent(out) :: status

        real(kind = rk), allocatable :: quad_x(:), check_x(:), scan_x(:)
        integer :: alloc_stat

        if (max_params > CASSINI_PARAM_LIMIT) then
            status = SHAPE_ERROR_TOO_MANY_PARAMS
            return
        end if
        if (max_params < 1_ik) then
            status = SHAPE_ERROR_INVALID_INIT
            return
        end if
        if (n_x_grid < CASSINI_N_X_GRID_FLOOR) then
            status = SHAPE_ERROR_INVALID_GRID
            return
        end if
        call check_thetas_s(thetas, status)
        if (status /= SHAPE_VALID) return

        allocate(cache%thetas(size(thetas)), cache%quad_w(2_ik * CASSINI_N_GL_HALF), &
                quad_x(2_ik * CASSINI_N_GL_HALF), check_x(CASSINI_N_R_CHECK), scan_x(n_x_grid), &
                stat = alloc_stat)
        if (alloc_stat /= 0) then
            call cache_free_s(cache)
            status = SHAPE_ERROR_INVALID_GRID
            return
        end if

        cache%thetas = thetas
        call quadrature_nodes_s(quad_x, cache%quad_w)
        call cosine_nodes_s(CASSINI_N_R_CHECK, check_x)
        call cosine_nodes_s(n_x_grid, scan_x)
        call build_tables_s(quad_x, max_params - 1_ik, cache%quad_tab)
        call build_tables_s(check_x, max_params - 1_ik, cache%check_tab)
        call build_tables_s(scan_x, max_params - 1_ik, cache%scan_tab)
        if (cache%quad_tab%order < 0_ik .or. cache%check_tab%order < 0_ik &
                .or. cache%scan_tab%order < 0_ik) then
            call cache_free_s(cache)
            status = SHAPE_ERROR_INVALID_GRID
            return
        end if

        cache%max_params = max_params
        cache%n_x_grid = n_x_grid
        cache%initialized = .true.
        status = SHAPE_VALID

    end subroutine cache_init_s

    !> Releases the cache. Infallible, and safe on an already-free instance.
    pure subroutine cache_free_s(cache)
        type(cache_t), intent(inout) :: cache
        if (allocated(cache%thetas)) deallocate(cache%thetas)
        if (allocated(cache%quad_w)) deallocate(cache%quad_w)
        cache%quad_tab = tables_t()
        cache%check_tab = tables_t()
        cache%scan_tab = tables_t()
        cache%max_params = 0_ik
        cache%n_x_grid = 0_ik
        cache%initialized = .false.
    end subroutine cache_free_s

    !> The cache's max_params; 0 when uninitialized.
    pure function cache_max_params_f(cache) result(max_params)
        type(cache_t), intent(in) :: cache
        integer(kind = ik) :: max_params
        max_params = cache%max_params
    end function cache_max_params_f

    !> The cache's scan resolution; 0 when uninitialized.
    pure function cache_n_x_grid_f(cache) result(n_x_grid)
        type(cache_t), intent(in) :: cache
        integer(kind = ik) :: n_x_grid
        n_x_grid = cache%n_x_grid
    end function cache_n_x_grid_f

    !> Number of primary thetas of the cache; 0 when uninitialized.
    pure function cache_n_thetas_f(cache) result(n_thetas)
        type(cache_t), intent(in) :: cache
        integer(kind = ik) :: n_thetas
        n_thetas = 0_ik
        if (cache%initialized) n_thetas = size(cache%thetas, kind = ik)
    end function cache_n_thetas_f

    !> .true. after a successful `cache_init_s`, .false. after `cache_free_s`.
    pure function cache_is_initialized_f(cache) result(is_initialized)
        type(cache_t), intent(in) :: cache
        logical :: is_initialized
        is_initialized = cache%initialized
    end function cache_is_initialized_f

    !===========================================================================
    ! CACHED TIER
    !===========================================================================

    !> R(theta) at the cache's primary thetas.
    !!
    !! @param[out] radii   cache_n_thetas_f(cache) long
    !! @param[out] status  0, 2, 4, 3, 100, 101, 103 or 102
    pure subroutine cache_radius_grid_s(cache, params, radii, status)
        type(cache_t), intent(in) :: cache
        real(kind = rk), intent(in) :: params(:)
        real(kind = rk), intent(out) :: radii(:)
        integer(kind = ik), intent(out) :: status
        real(kind = rk), allocatable :: dr_sink(:)
        integer :: alloc_stat
        radii = 0.0_rk
        call check_call_s(cache, params, status)
        if (status /= SHAPE_VALID) return
        if (size(radii) /= size(cache%thetas)) then
            status = SHAPE_ERROR_INVALID_GRID
            return
        end if
        allocate(dr_sink(size(radii)), stat = alloc_stat)
        if (alloc_stat /= 0) then
            status = SHAPE_ERROR_INVALID_GRID
            return
        end if
        call rtheta_s(cache, params, cache%thetas, radii, dr_sink, status)
    end subroutine cache_radius_grid_s

    !> R(theta) and dR/dtheta at the cache's primary thetas.
    !!
    !! @param[out] radii, dr_dtheta  cache_n_thetas_f(cache) long each
    !! @param[out] status            0, 2, 4, 3, 100, 101, 103 or 102
    pure subroutine cache_radius_and_derivative_s(cache, params, radii, dr_dtheta, status)
        type(cache_t), intent(in) :: cache
        real(kind = rk), intent(in) :: params(:)
        real(kind = rk), intent(out) :: radii(:), dr_dtheta(:)
        integer(kind = ik), intent(out) :: status
        radii = 0.0_rk
        dr_dtheta = 0.0_rk
        call check_call_s(cache, params, status)
        if (status /= SHAPE_VALID) return
        if (size(radii) /= size(cache%thetas) .or. size(dr_dtheta) /= size(cache%thetas)) then
            status = SHAPE_ERROR_INVALID_GRID
            return
        end if
        call rtheta_s(cache, params, cache%thetas, radii, dr_dtheta, status)
    end subroutine cache_radius_and_derivative_s

    !> R(theta) and dR/dtheta at caller thetas (the WMMM hot path).
    !!
    !! @param[in]  thetas            At least one, each in [0, pi]
    !! @param[out] radii, dr_dtheta  size(thetas) long each
    !! @param[out] status            0, 2, 4, 3, 100, 101, 103 or 102
    pure subroutine cache_radius_and_derivative_at_thetas_s(cache, params, thetas, radii, &
            dr_dtheta, status)
        type(cache_t), intent(in) :: cache
        real(kind = rk), intent(in) :: params(:), thetas(:)
        real(kind = rk), intent(out) :: radii(:), dr_dtheta(:)
        integer(kind = ik), intent(out) :: status
        radii = 0.0_rk
        dr_dtheta = 0.0_rk
        call check_call_s(cache, params, status)
        if (status /= SHAPE_VALID) return
        if (size(radii) /= size(thetas) .or. size(dr_dtheta) /= size(thetas)) then
            status = SHAPE_ERROR_INVALID_GRID
            return
        end if
        call check_thetas_s(thetas, status)
        if (status /= SHAPE_VALID) return
        call rtheta_s(cache, params, thetas, radii, dr_dtheta, status)
    end subroutine cache_radius_and_derivative_at_thetas_s

    !> z_shift and the pole distances, with the R(theta)-path verdict.
    !!
    !! @param[out] status  0, 2, 4, 3, 100, 101, 103 or 102; on 102 and 103
    !!                     `cache_geometry_s` still has the poles
    pure subroutine cache_shape_s(cache, params, z_shift, r_north, r_south, status)
        type(cache_t), intent(in) :: cache
        real(kind = rk), intent(in) :: params(:)
        real(kind = rk), intent(out) :: z_shift, r_north, r_south
        integer(kind = ik), intent(out) :: status
        type(resolve_t) :: res
        real(kind = rk) :: d_min
        logical :: is_folded
        z_shift = 0.0_rk
        r_north = 0.0_rk
        r_south = 0.0_rk
        call check_call_s(cache, params, status)
        if (status /= SHAPE_VALID) return
        call run_resolve_s(cache, params, res, status)
        if (status /= SHAPE_VALID) return
        call run_margin_s(cache, params, res, d_min, is_folded, status)
        if (status /= SHAPE_VALID) return
        call gate_rtheta_s(res, d_min, status)
        if (status /= SHAPE_VALID) return
        z_shift = res%z_shift
        r_north = res%r_north
        r_south = res%r_south
    end subroutine cache_shape_s

    !> The parametric profile (z, rho, drho/dz), checked: 103 when separated.
    !!
    !! @param[out] z, rho, drho_dz  Equal sizes >= CASSINI_N_POINTS_FLOOR
    !! @param[out] status           0, 2, 4, 3, 100, 101 or 103
    pure subroutine cache_rho_profile_s(cache, params, z, rho, drho_dz, status)
        type(cache_t), intent(in) :: cache
        real(kind = rk), intent(in) :: params(:)
        real(kind = rk), intent(out) :: z(:), rho(:), drho_dz(:)
        integer(kind = ik), intent(out) :: status
        call profile_out_s(cache, params, .true., z, rho, drho_dz, status)
    end subroutine cache_rho_profile_s

    !> The parametric profile, unchecked: draws separated shapes (two lobes).
    !!
    !! @param[out] status  0, 2, 4, 3, 100 or 101
    pure subroutine cache_rho_profile_unchecked_s(cache, params, z, rho, drho_dz, status)
        type(cache_t), intent(in) :: cache
        real(kind = rk), intent(in) :: params(:)
        real(kind = rk), intent(out) :: z(:), rho(:), drho_dz(:)
        integer(kind = ik), intent(out) :: status
        call profile_out_s(cache, params, .false., z, rho, drho_dz, status)
    end subroutine cache_rho_profile_unchecked_s

    !> The neck: found with z_neck and rho_neck, or not found with both 0.
    !! Pinched and separated shapes: found, rho_neck = 0, z_neck = -z_shift.
    !!
    !! @param[out] status  0, 2, 4, 3, 100 or 101
    pure subroutine cache_neck_s(cache, params, z_neck, rho_neck, found, status)
        type(cache_t), intent(in) :: cache
        real(kind = rk), intent(in) :: params(:)
        real(kind = rk), intent(out) :: z_neck, rho_neck
        logical, intent(out) :: found
        integer(kind = ik), intent(out) :: status
        type(resolve_t) :: res
        real(kind = rk), allocatable :: rho(:), z(:), drho_dx(:), dz_dx(:)
        integer :: alloc_stat
        z_neck = 0.0_rk
        rho_neck = 0.0_rk
        found = .false.
        call check_call_s(cache, params, status)
        if (status /= SHAPE_VALID) return
        call run_resolve_s(cache, params, res, status)
        if (status /= SHAPE_VALID) return
        allocate(rho(cache%n_x_grid), z(cache%n_x_grid), drho_dx(cache%n_x_grid), &
                dz_dx(cache%n_x_grid), stat = alloc_stat)
        if (alloc_stat /= 0) then
            status = SHAPE_ERROR_INVALID_GRID
            return
        end if
        call curve_on_tables_s(params, res, cache%scan_tab, rho, z, drho_dx, dz_dx)
        call neck_s(params, res, cache%scan_tab, rho, found, z_neck, rho_neck)
    end subroutine cache_neck_s

    !> The resolved geometry record.
    !!
    !! @param[out] status  0, 2, 4, 3, 100 or 101
    pure subroutine cache_geometry_s(cache, params, geometry, status)
        type(cache_t), intent(in) :: cache
        real(kind = rk), intent(in) :: params(:)
        type(cassini_geometry_t), intent(out) :: geometry
        integer(kind = ik), intent(out) :: status
        type(resolve_t) :: res
        real(kind = rk) :: d_min
        logical :: is_folded
        call check_call_s(cache, params, status)
        if (status /= SHAPE_VALID) return
        call run_resolve_s(cache, params, res, status)
        if (status /= SHAPE_VALID) return
        call run_margin_s(cache, params, res, d_min, is_folded, status)
        if (status /= SHAPE_VALID) return
        geometry%eps = res%eps
        geometry%c = res%c
        geometry%z_shift = res%z_shift
        geometry%r_north = res%r_north
        geometry%r_south = res%r_south
        geometry%gap = res%gap
        geometry%d_min = d_min
        geometry%is_separated = res%is_separated
        geometry%is_pinched = res%is_pinched
        geometry%is_folded = is_folded
    end subroutine cache_geometry_s

    !> The quadrupole moment q2 = <r^2 P2(cos theta)> in R0 = 1 units.
    !!
    !! @param[out] status  0, 2, 4, 100 or 101
    pure subroutine cache_quadrupole_s(cache, params, q2, status)
        type(cache_t), intent(in) :: cache
        real(kind = rk), intent(in) :: params(:)
        real(kind = rk), intent(out) :: q2
        integer(kind = ik), intent(out) :: status
        type(resolve_t) :: res
        q2 = 0.0_rk
        call check_call_s(cache, params, status)
        if (status /= SHAPE_VALID) return
        call run_resolve_s(cache, params, res, status)
        if (status /= SHAPE_VALID) return
        q2 = quadrupole_f(params, res, cache%quad_tab, cache%quad_w)
    end subroutine cache_quadrupole_s

    !===========================================================================
    ! ONE-SHOT TIER
    !===========================================================================

    !> One-shot `cache_radius_grid_s` at caller thetas.
    pure subroutine compute_radius_grid_standalone_s(params, n_x_grid, thetas, radii, status)
        real(kind = rk), intent(in) :: params(:), thetas(:)
        integer(kind = ik), intent(in) :: n_x_grid
        real(kind = rk), intent(out) :: radii(:)
        integer(kind = ik), intent(out) :: status
        type(cache_t) :: cache
        radii = 0.0_rk
        call check_one_shot_s(params, status)
        if (status /= SHAPE_VALID) return
        if (size(radii) /= size(thetas)) then
            status = SHAPE_ERROR_INVALID_GRID
            return
        end if
        call cache_init_s(cache, size(params, kind = ik), n_x_grid, thetas, status)
        if (status /= SHAPE_VALID) return
        call cache_radius_grid_s(cache, params, radii, status)
    end subroutine compute_radius_grid_standalone_s

    !> One-shot `cache_radius_and_derivative_s` at caller thetas.
    pure subroutine compute_radius_and_derivative_standalone_s(params, n_x_grid, thetas, radii, &
            dr_dtheta, status)
        real(kind = rk), intent(in) :: params(:), thetas(:)
        integer(kind = ik), intent(in) :: n_x_grid
        real(kind = rk), intent(out) :: radii(:), dr_dtheta(:)
        integer(kind = ik), intent(out) :: status
        type(cache_t) :: cache
        radii = 0.0_rk
        dr_dtheta = 0.0_rk
        call check_one_shot_s(params, status)
        if (status /= SHAPE_VALID) return
        if (size(radii) /= size(thetas) .or. size(dr_dtheta) /= size(thetas)) then
            status = SHAPE_ERROR_INVALID_GRID
            return
        end if
        call cache_init_s(cache, size(params, kind = ik), n_x_grid, thetas, status)
        if (status /= SHAPE_VALID) return
        call cache_radius_and_derivative_s(cache, params, radii, dr_dtheta, status)
    end subroutine compute_radius_and_derivative_standalone_s

    !> One-shot `cache_shape_s`.
    pure subroutine compute_shape_standalone_s(params, n_x_grid, z_shift, r_north, r_south, status)
        real(kind = rk), intent(in) :: params(:)
        integer(kind = ik), intent(in) :: n_x_grid
        real(kind = rk), intent(out) :: z_shift, r_north, r_south
        integer(kind = ik), intent(out) :: status
        type(cache_t) :: cache
        z_shift = 0.0_rk
        r_north = 0.0_rk
        r_south = 0.0_rk
        call check_one_shot_s(params, status)
        if (status /= SHAPE_VALID) return
        call cache_init_s(cache, size(params, kind = ik), n_x_grid, ONE_THETA, status)
        if (status /= SHAPE_VALID) return
        call cache_shape_s(cache, params, z_shift, r_north, r_south, status)
    end subroutine compute_shape_standalone_s

    !> One-shot `cache_rho_profile_s`; the profile size is the arrays' size.
    pure subroutine compute_rho_profile_standalone_s(params, z, rho, drho_dz, status)
        real(kind = rk), intent(in) :: params(:)
        real(kind = rk), intent(out) :: z(:), rho(:), drho_dz(:)
        integer(kind = ik), intent(out) :: status
        call one_shot_profile_s(params, .true., z, rho, drho_dz, status)
    end subroutine compute_rho_profile_standalone_s

    !> One-shot `cache_rho_profile_unchecked_s`.
    pure subroutine compute_rho_profile_unchecked_standalone_s(params, z, rho, drho_dz, status)
        real(kind = rk), intent(in) :: params(:)
        real(kind = rk), intent(out) :: z(:), rho(:), drho_dz(:)
        integer(kind = ik), intent(out) :: status
        call one_shot_profile_s(params, .false., z, rho, drho_dz, status)
    end subroutine compute_rho_profile_unchecked_standalone_s

    !> One-shot `cache_neck_s`.
    pure subroutine compute_neck_standalone_s(params, n_x_grid, z_neck, rho_neck, found, status)
        real(kind = rk), intent(in) :: params(:)
        integer(kind = ik), intent(in) :: n_x_grid
        real(kind = rk), intent(out) :: z_neck, rho_neck
        logical, intent(out) :: found
        integer(kind = ik), intent(out) :: status
        type(cache_t) :: cache
        z_neck = 0.0_rk
        rho_neck = 0.0_rk
        found = .false.
        call check_one_shot_s(params, status)
        if (status /= SHAPE_VALID) return
        call cache_init_s(cache, size(params, kind = ik), n_x_grid, ONE_THETA, status)
        if (status /= SHAPE_VALID) return
        call cache_neck_s(cache, params, z_neck, rho_neck, found, status)
    end subroutine compute_neck_standalone_s

    !> One-shot `cache_geometry_s`.
    pure subroutine compute_geometry_standalone_s(params, n_x_grid, geometry, status)
        real(kind = rk), intent(in) :: params(:)
        integer(kind = ik), intent(in) :: n_x_grid
        type(cassini_geometry_t), intent(out) :: geometry
        integer(kind = ik), intent(out) :: status
        type(cache_t) :: cache
        call check_one_shot_s(params, status)
        if (status /= SHAPE_VALID) return
        call cache_init_s(cache, size(params, kind = ik), n_x_grid, ONE_THETA, status)
        if (status /= SHAPE_VALID) return
        call cache_geometry_s(cache, params, geometry, status)
    end subroutine compute_geometry_standalone_s

    !> One-shot `cache_quadrupole_s`.
    pure subroutine compute_quadrupole_standalone_s(params, q2, status)
        real(kind = rk), intent(in) :: params(:)
        real(kind = rk), intent(out) :: q2
        integer(kind = ik), intent(out) :: status
        type(cache_t) :: cache
        q2 = 0.0_rk
        call check_one_shot_s(params, status)
        if (status /= SHAPE_VALID) return
        call cache_init_s(cache, size(params, kind = ik), CASSINI_N_X_GRID_FLOOR, ONE_THETA, status)
        if (status /= SHAPE_VALID) return
        call cache_quadrupole_s(cache, params, q2, status)
    end subroutine compute_quadrupole_standalone_s

    !> eps from params, infallible. Accepts 0 .. 21 parameters (the empty
    !! vector is the sphere, eps = 0); a longer vector is a caller error and
    !! only its first 21 entries are read.
    pure function compute_epsilon_f(params) result(eps)
        real(kind = rk), intent(in) :: params(:)
        real(kind = rk) :: eps
        eps = epsilon_f(params, min(n_active_f(params), CASSINI_MAX_PARAMS))
    end function compute_epsilon_f

    !===========================================================================
    ! INTERNALS
    !===========================================================================

    !> Cached-call usage checks: 2, then 4.
    pure subroutine check_call_s(cache, params, status)
        type(cache_t), intent(in) :: cache
        real(kind = rk), intent(in) :: params(:)
        integer(kind = ik), intent(out) :: status
        if (.not. cache%initialized) then
            status = SHAPE_ERROR_CACHE_NOT_INITIALIZED
            return
        end if
        if (size(params, kind = ik) < 1_ik .or. size(params, kind = ik) > cache%max_params) then
            status = SHAPE_ERROR_WRONG_PARAM_COUNT
            return
        end if
        status = SHAPE_VALID
    end subroutine check_call_s

    !> One-shot usage checks: an empty vector (4), then one longer than the
    !! library limit (1) -- never silently truncated.
    pure subroutine check_one_shot_s(params, status)
        real(kind = rk), intent(in) :: params(:)
        integer(kind = ik), intent(out) :: status
        if (size(params, kind = ik) < 1_ik) then
            status = SHAPE_ERROR_WRONG_PARAM_COUNT
            return
        end if
        if (size(params, kind = ik) > CASSINI_PARAM_LIMIT) then
            status = SHAPE_ERROR_TOO_MANY_PARAMS
            return
        end if
        status = SHAPE_VALID
    end subroutine check_one_shot_s

    !> A theta set: at least one node, every node in [0, pi].
    pure subroutine check_thetas_s(thetas, status)
        real(kind = rk), intent(in) :: thetas(:)
        integer(kind = ik), intent(out) :: status
        status = SHAPE_VALID
        if (size(thetas) < 1) then
            status = SHAPE_ERROR_INVALID_GRID
            return
        end if
        if (any(thetas < 0.0_rk) .or. any(thetas > PI_C)) status = SHAPE_ERROR_INVALID_GRID
    end subroutine check_thetas_s

    !> Profile arrays: equal sizes, at least CASSINI_N_POINTS_FLOOR.
    pure subroutine check_profile_s(z, rho, drho_dz, status)
        real(kind = rk), intent(in) :: z(:), rho(:), drho_dz(:)
        integer(kind = ik), intent(out) :: status
        status = SHAPE_VALID
        if (size(rho) /= size(z) .or. size(drho_dz) /= size(z) &
                .or. size(z, kind = ik) < CASSINI_N_POINTS_FLOOR) status = SHAPE_ERROR_INVALID_GRID
    end subroutine check_profile_s

    !> The resolve on the cache's tables: 0, 100 or 101.
    pure subroutine run_resolve_s(cache, params, res, status)
        type(cache_t), intent(in) :: cache
        real(kind = rk), intent(in) :: params(:)
        type(resolve_t), intent(out) :: res
        integer(kind = ik), intent(out) :: status
        call resolve_s(params, cache%quad_tab, cache%quad_w, cache%check_tab, res, status)
    end subroutine run_resolve_s

    !> The star margin and fold flag on the cache's scan; 3 if the scan
    !! buffers cannot be allocated.
    pure subroutine run_margin_s(cache, params, res, d_min, is_folded, status)
        type(cache_t), intent(in) :: cache
        real(kind = rk), intent(in) :: params(:)
        type(resolve_t), intent(in) :: res
        real(kind = rk), intent(out) :: d_min
        logical, intent(out) :: is_folded
        integer(kind = ik), intent(out) :: status
        real(kind = rk), allocatable :: rho(:), z(:), drho_dx(:), dz_dx(:)
        integer :: alloc_stat
        d_min = 0.0_rk
        is_folded = .false.
        allocate(rho(cache%n_x_grid), z(cache%n_x_grid), drho_dx(cache%n_x_grid), &
                dz_dx(cache%n_x_grid), stat = alloc_stat)
        if (alloc_stat /= 0) then
            status = SHAPE_ERROR_INVALID_GRID
            return
        end if
        call curve_on_tables_s(params, res, cache%scan_tab, rho, z, drho_dx, dz_dx)
        call star_margin_s(params, res, cache%scan_tab%x, rho, z, drho_dx, dz_dx, d_min, is_folded)
        status = SHAPE_VALID
    end subroutine run_margin_s

    !> The R(theta)-path gate: 103 when separated, then 102 when pinched or
    !! when d_min < STAR_CONVEXITY_MARGIN.
    pure subroutine gate_rtheta_s(res, d_min, status)
        type(resolve_t), intent(in) :: res
        real(kind = rk), intent(in) :: d_min
        integer(kind = ik), intent(out) :: status
        status = SHAPE_VALID
        if (res%is_separated) then
            status = CASSINI_ERROR_SEPARATED
        else if (res%is_pinched .or. d_min < STAR_CONVEXITY_MARGIN) then
            status = CASSINI_ERROR_NOT_STAR_CONVEX
        end if
    end subroutine gate_rtheta_s

    !> The R(theta) pipeline after the usage checks: resolve, scan, gate,
    !! then the root-find at every theta. The outputs are written only on
    !! success (the callers zero-filled them).
    pure subroutine rtheta_s(cache, params, thetas, radii, dr_dtheta, status)
        type(cache_t), intent(in) :: cache
        real(kind = rk), intent(in) :: params(:), thetas(:)
        real(kind = rk), intent(inout) :: radii(:), dr_dtheta(:)
        integer(kind = ik), intent(out) :: status
        type(resolve_t) :: res
        real(kind = rk), allocatable :: rho(:), z(:), drho_dx(:), dz_dx(:), theta_scan(:), dtheta_dx(:)
        real(kind = rk) :: d_min
        logical :: is_folded
        integer :: alloc_stat
        call run_resolve_s(cache, params, res, status)
        if (status /= SHAPE_VALID) return
        allocate(rho(cache%n_x_grid), z(cache%n_x_grid), drho_dx(cache%n_x_grid), &
                dz_dx(cache%n_x_grid), theta_scan(cache%n_x_grid), dtheta_dx(cache%n_x_grid), &
                stat = alloc_stat)
        if (alloc_stat /= 0) then
            status = SHAPE_ERROR_INVALID_GRID
            return
        end if
        call curve_on_tables_s(params, res, cache%scan_tab, rho, z, drho_dx, dz_dx)
        call star_margin_s(params, res, cache%scan_tab%x, rho, z, drho_dx, dz_dx, d_min, is_folded)
        call gate_rtheta_s(res, d_min, status)
        if (status /= SHAPE_VALID) return
        call theta_on_scan_s(rho, z, drho_dx, dz_dx, theta_scan, dtheta_dx)
        call radius_at_thetas_s(params, res, cache%scan_tab%x, theta_scan, thetas, radii, dr_dtheta)
    end subroutine rtheta_s

    !> Both profile forms after zero-filling: 2, 4, 3, 100, 101, then 103
    !! when checked and separated.
    pure subroutine profile_out_s(cache, params, checked, z, rho, drho_dz, status)
        type(cache_t), intent(in) :: cache
        real(kind = rk), intent(in) :: params(:)
        logical, intent(in) :: checked
        real(kind = rk), intent(out) :: z(:), rho(:), drho_dz(:)
        integer(kind = ik), intent(out) :: status
        type(resolve_t) :: res
        z = 0.0_rk
        rho = 0.0_rk
        drho_dz = 0.0_rk
        call check_call_s(cache, params, status)
        if (status /= SHAPE_VALID) return
        call check_profile_s(z, rho, drho_dz, status)
        if (status /= SHAPE_VALID) return
        call run_resolve_s(cache, params, res, status)
        if (status /= SHAPE_VALID) return
        if (checked .and. res%is_separated) then
            status = CASSINI_ERROR_SEPARATED
            return
        end if
        call profile_s(params, res, size(z, kind = ik), z, rho, drho_dz)
    end subroutine profile_out_s

    !> Both one-shot profile forms: a local cache on the floor scan grid and
    !! the one-node theta set (neither is consulted).
    pure subroutine one_shot_profile_s(params, checked, z, rho, drho_dz, status)
        real(kind = rk), intent(in) :: params(:)
        logical, intent(in) :: checked
        real(kind = rk), intent(out) :: z(:), rho(:), drho_dz(:)
        integer(kind = ik), intent(out) :: status
        type(cache_t) :: cache
        z = 0.0_rk
        rho = 0.0_rk
        drho_dz = 0.0_rk
        call check_one_shot_s(params, status)
        if (status /= SHAPE_VALID) return
        call check_profile_s(z, rho, drho_dz, status)
        if (status /= SHAPE_VALID) return
        call cache_init_s(cache, size(params, kind = ik), CASSINI_N_X_GRID_FLOOR, ONE_THETA, status)
        if (status /= SHAPE_VALID) return
        call profile_out_s(cache, params, checked, z, rho, drho_dz, status)
    end subroutine one_shot_profile_s

end module cassini_parameterization_mod
