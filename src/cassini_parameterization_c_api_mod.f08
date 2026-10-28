!> C-interop layer for the Cassini-oval shape library (`cassini_parameterization.h`).
!!
!! One opaque handle: `cassini_param_cache_create` returns the `c_loc` of a
!! heap-allocated `cache_t`. The cache is immutable after creation, so every
!! compute takes the handle `const` and any number of threads may compute on
!! one handle concurrently. Create and destroy must not race with any other
!! call on the same handle.
!!
!! `cassini_param_cache_create` returns a null handle on failure and reports
!! the rejecting code through its nullable `status`. The cached computes
!! return the status directly; the one-shot entries report through a
!! trailing nullable `int32_t* status` (an absent `optional` dummy when C
!! passes NULL). Diagnostics are the static strings behind
!! `cassini_param_status_message`.
!!
!! This layer marshals. The parameter-count, extent and grid checks live in
!! the Fortran tiers and their codes pass through untouched. The checks that
!! live HERE are the ones the Fortran layer cannot see: the NULL-handle guard,
!! the negative-size guard on the C integer arguments (a negative size in a
!! cached call clamps to an empty extent, which the tier rejects as 4 for
!! params and 3 for an output; in a one-shot call a negative n_params,
!! n_thetas or n_z is a bad argument, 5; a negative n_x_grid anywhere reaches
!! the tier, which rejects it as a grid below the floor, 3), and the
!! allocation guard.
!!
!! Arrays arrive assumed-size and are read and written only over their stated
!! extent, so a NULL pointer with a zero or negative size is never touched.
!! Marshalling buffers are heap allocatables with `stat =`, never automatic
!! arrays: a wrong but huge size argument then fails recoverably (3 in a
!! cached call, 5 in a one-shot call) instead of overflowing the stack before
!! the tier can judge it. A stated size must be the ACTUAL extent of the
!! caller's buffer: outputs are zero-filled over exactly the stated extent.
module cassini_parameterization_c_api_mod

    use, intrinsic :: iso_c_binding, only: c_ptr, c_loc, c_f_pointer, c_associated, &
            c_null_ptr, c_char, c_null_char
    use c_bindings_mod, only: ik_c, rk_c
    use precision_utilities_mod, only: ik, rk
    use cassini_parameterization_mod, only: cache_t, cassini_geometry_t, cache_init_s, &
            cache_free_s, cache_radius_grid_s, cache_radius_and_derivative_s, &
            cache_radius_and_derivative_at_thetas_s, cache_shape_s, cache_rho_profile_s, &
            cache_rho_profile_unchecked_s, cache_neck_s, cache_geometry_s, cache_quadrupole_s, &
            compute_radius_grid_standalone_s, compute_radius_and_derivative_standalone_s, &
            compute_shape_standalone_s, compute_rho_profile_standalone_s, &
            compute_rho_profile_unchecked_standalone_s, compute_neck_standalone_s, &
            compute_geometry_standalone_s, compute_quadrupole_standalone_s, compute_epsilon_f, &
            STATUS_MESSAGE_LEN, CASSINI_PARAM_LIMIT, SHAPE_VALID, SHAPE_ERROR_TOO_MANY_PARAMS, &
            SHAPE_ERROR_CACHE_NOT_INITIALIZED, SHAPE_ERROR_INVALID_GRID, &
            SHAPE_ERROR_WRONG_PARAM_COUNT, SHAPE_ERROR_INVALID_INIT, CASSINI_ERROR_R_NONPOSITIVE, &
            CASSINI_ERROR_TOROIDAL, CASSINI_ERROR_NOT_STAR_CONVEX, CASSINI_ERROR_SEPARATED

    implicit none

    private

    public :: cassini_param_status_message
    public :: cassini_param_epsilon
    public :: cassini_param_cache_create
    public :: cassini_param_cache_destroy
    public :: cassini_param_cache_radius_grid
    public :: cassini_param_cache_radius_and_derivative
    public :: cassini_param_cache_radius_and_derivative_at_thetas
    public :: cassini_param_cache_shape
    public :: cassini_param_cache_rho_profile
    public :: cassini_param_cache_rho_profile_unchecked
    public :: cassini_param_cache_neck
    public :: cassini_param_cache_geometry
    public :: cassini_param_cache_quadrupole
    public :: cassini_param_radius_grid
    public :: cassini_param_radius_and_derivative
    public :: cassini_param_shape
    public :: cassini_param_rho_profile
    public :: cassini_param_rho_profile_unchecked
    public :: cassini_param_neck
    public :: cassini_param_geometry
    public :: cassini_param_quadrupole

    !> The C `cassini_param_geometry_t`: the components of `cassini_geometry_t`
    !! in the same order, flags as 0/1 ints.
    type, bind(c), public :: cassini_param_geometry_t
        real(rk_c) :: eps
        real(rk_c) :: c
        real(rk_c) :: z_shift
        real(rk_c) :: r_north
        real(rk_c) :: r_south
        real(rk_c) :: gap
        real(rk_c) :: d_min
        integer(ik_c) :: is_separated
        integer(ik_c) :: is_pinched
        integer(ik_c) :: is_folded
    end type cassini_param_geometry_t

    !---------------------------------------------------------------------------
    ! Static status strings
    !---------------------------------------------------------------------------
    !> One null-terminated string per status code, laid out as fixed columns of
    !! a module array initialized at compile time and never written again, so
    !! the pointer `cassini_param_status_message` hands out is valid forever and
    !! every thread reads the same bytes. The texts are the `status_message`
    !! texts: keep the two in sync.
    integer(kind = ik), parameter :: MSG_LEN = STATUS_MESSAGE_LEN + 1_ik
    integer(kind = ik), parameter :: N_MSG = 11_ik

    !> Codes in column order; column N_MSG is the unknown-code fallback.
    integer(kind = ik), parameter :: MSG_CODE(N_MSG - 1_ik) = [ &
            SHAPE_VALID, SHAPE_ERROR_TOO_MANY_PARAMS, SHAPE_ERROR_CACHE_NOT_INITIALIZED, &
            SHAPE_ERROR_INVALID_GRID, SHAPE_ERROR_WRONG_PARAM_COUNT, SHAPE_ERROR_INVALID_INIT, &
            CASSINI_ERROR_R_NONPOSITIVE, CASSINI_ERROR_TOROIDAL, CASSINI_ERROR_NOT_STAR_CONVEX, &
            CASSINI_ERROR_SEPARATED]

    character(kind = c_char, len = MSG_LEN), parameter :: MSG_TEXT(N_MSG) = &
            [character(kind = c_char, len = MSG_LEN) :: &
                    'valid' // c_null_char, &
                    'too many parameters' // c_null_char, &
                    'cache not initialized' // c_null_char, &
                    'invalid grid: thetas, scan or profile size, or buffer extent' // c_null_char, &
                    'params length is 0 or exceeds max_params' // c_null_char, &
                    'invalid init arguments' // c_null_char, &
                    'R(x) <= 0 on the check grid' // c_null_char, &
                    'toroidal: a pole has left the symmetry axis' // c_null_char, &
                    'shape not star-convex from the centre of mass' // c_null_char, &
                    'shape separated into two fragments' // c_null_char, &
                    'unknown status code' // c_null_char]

    character(kind = c_char), save, target :: MSG_CHARS(MSG_LEN, N_MSG) = &
            reshape(transfer(MSG_TEXT, c_null_char, MSG_LEN * N_MSG), [MSG_LEN, N_MSG])

contains

    !===========================================================================
    ! INTERNAL HELPERS
    !===========================================================================

    !> Non-negative Fortran extent for a C size argument: a negative size
    !! becomes an empty extent the Fortran tier rejects.
    pure function extent_f(n) result(extent)
        integer(ik_c), intent(in) :: n
        integer(kind = ik) :: extent
        extent = max(int(n, ik), 0_ik)
    end function extent_f

    !> Copy of the first `extent_f(n)` C doubles; ok = .false. if the heap refused.
    subroutine copy_in_s(values, n, values_f, ok)
        real(rk_c), intent(in) :: values(*)
        integer(ik_c), intent(in) :: n
        real(kind = rk), allocatable, intent(out) :: values_f(:)
        logical, intent(out) :: ok
        integer :: alloc_stat
        integer(kind = ik) :: m
        m = extent_f(n)
        allocate(values_f(m), stat = alloc_stat)
        ok = alloc_stat == 0
        if (ok .and. m > 0_ik) values_f = real(values(1:m), rk)
    end subroutine copy_in_s

    !> One to three output work buffers of extent_f(n); ok = .false. if the
    !! heap refused any of them.
    subroutine work_buffers_s(n, ok, a, b, c)
        integer(ik_c), intent(in) :: n
        logical, intent(out) :: ok
        real(kind = rk), allocatable, intent(out) :: a(:)
        real(kind = rk), allocatable, intent(out), optional :: b(:), c(:)
        integer :: stat_a, stat_b, stat_c
        stat_b = 0
        stat_c = 0
        allocate(a(extent_f(n)), stat = stat_a)
        if (present(b)) allocate(b(extent_f(n)), stat = stat_b)
        if (present(c)) allocate(c(extent_f(n)), stat = stat_c)
        ok = stat_a == 0 .and. stat_b == 0 .and. stat_c == 0
    end subroutine work_buffers_s

    !> Field-by-field copy into the C record; flags as 0/1.
    pure subroutine geometry_to_c(g, gc)
        type(cassini_geometry_t), intent(in) :: g
        type(cassini_param_geometry_t), intent(out) :: gc
        gc%eps = real(g%eps, rk_c)
        gc%c = real(g%c, rk_c)
        gc%z_shift = real(g%z_shift, rk_c)
        gc%r_north = real(g%r_north, rk_c)
        gc%r_south = real(g%r_south, rk_c)
        gc%gap = real(g%gap, rk_c)
        gc%d_min = real(g%d_min, rk_c)
        gc%is_separated = merge(1_ik_c, 0_ik_c, g%is_separated)
        gc%is_pinched = merge(1_ik_c, 0_ik_c, g%is_pinched)
        gc%is_folded = merge(1_ik_c, 0_ik_c, g%is_folded)
    end subroutine geometry_to_c

    !> The all-zero C record of a failed call.
    pure subroutine zero_geometry_c(gc)
        type(cassini_param_geometry_t), intent(out) :: gc
        gc = cassini_param_geometry_t(0.0_rk_c, 0.0_rk_c, 0.0_rk_c, 0.0_rk_c, 0.0_rk_c, 0.0_rk_c, &
                0.0_rk_c, 0_ik_c, 0_ik_c, 0_ik_c)
    end subroutine zero_geometry_c

    !===========================================================================
    ! DIAGNOSTICS AND EPSILON
    !===========================================================================

    !> Pointer to the static, null-terminated description of `status`.
    !! Unknown codes get the fallback column. Never returns NULL.
    function cassini_param_status_message(status) result(msg) &
            bind(c, name = 'cassini_param_status_message')
        integer(ik_c), value, intent(in) :: status
        type(c_ptr) :: msg
        integer(kind = ik) :: i, idx
        idx = N_MSG
        do i = 1_ik, N_MSG - 1_ik
            if (status == int(MSG_CODE(i), ik_c)) then
                idx = i
                exit
            end if
        end do
        msg = c_loc(MSG_CHARS(1_ik, idx))
    end function cassini_param_status_message

    !> eps from params (0 .. 21 of them; the empty vector is the sphere). A
    !! negative n_params, more than 21, or a buffer the heap refuses returns 0.0.
    function cassini_param_epsilon(params, n_params) result(eps) &
            bind(c, name = 'cassini_param_epsilon')
        integer(ik_c), value, intent(in) :: n_params
        real(rk_c), intent(in) :: params(*)
        real(rk_c) :: eps
        real(kind = rk), allocatable :: params_f(:)
        logical :: ok
        eps = 0.0_rk_c
        if (n_params < 0_ik_c .or. int(n_params, ik) > CASSINI_PARAM_LIMIT) return
        call copy_in_s(params, n_params, params_f, ok)
        if (.not. ok) return
        eps = real(compute_epsilon_f(params_f), rk_c)
    end function cassini_param_epsilon

    !===========================================================================
    ! CACHE LIFECYCLE
    !===========================================================================

    !> Create a read-only cache. Null on failure; `status` (nullable) receives
    !! the rejecting code, or 0. The tier's order: max_params > 21 -> 1, < 1 ->
    !! 5, then n_x_grid below the floor -> 3, then the theta set -> 3 (a
    !! negative n_thetas is an empty set). A theta buffer the heap cannot
    !! satisfy is judged as an empty set; failure to allocate the handle object
    !! itself reports 5.
    function cassini_param_cache_create(max_params, n_x_grid, thetas, n_thetas, status) &
            result(handle) bind(c, name = 'cassini_param_cache_create')
        integer(ik_c), value, intent(in) :: max_params, n_x_grid, n_thetas
        real(rk_c), intent(in) :: thetas(*)
        integer(ik_c), intent(out), optional :: status
        type(c_ptr) :: handle

        type(cache_t), pointer :: p
        type(cache_t) :: probe
        integer(kind = ik) :: st
        integer :: alloc_stat
        real(kind = rk), allocatable :: thetas_f(:)
        real(kind = rk) :: no_thetas(0)
        logical :: ok

        handle = c_null_ptr

        call copy_in_s(thetas, n_thetas, thetas_f, ok)
        if (.not. ok) then
            call cache_init_s(probe, int(max_params, ik), int(n_x_grid, ik), no_thetas, st)
            if (present(status)) status = int(st, ik_c)
            return
        end if

        allocate(p, stat = alloc_stat)
        if (alloc_stat /= 0) then
            if (present(status)) status = int(SHAPE_ERROR_INVALID_INIT, ik_c)
            return
        end if

        call cache_init_s(p, int(max_params, ik), int(n_x_grid, ik), thetas_f, st)
        if (st == SHAPE_VALID) then
            handle = c_loc(p)
        else
            deallocate(p)
        end if
        if (present(status)) status = int(st, ik_c)
    end function cassini_param_cache_create

    !> Release a cache. NULL-safe.
    subroutine cassini_param_cache_destroy(handle) bind(c, name = 'cassini_param_cache_destroy')
        type(c_ptr), value, intent(in) :: handle
        type(cache_t), pointer :: p
        if (.not. c_associated(handle)) return
        call c_f_pointer(handle, p)
        call cache_free_s(p)
        deallocate(p)
    end subroutine cassini_param_cache_destroy

    !===========================================================================
    ! CACHED COMPUTES
    !===========================================================================
    ! Each entry zero-fills its outputs over the stated extent, guards the
    ! handle (2), copies in (a heap refusal is 3: the stated size is what is
    ! implicated) and hands the tier the copies.

    function cassini_param_cache_radius_grid(cache, params, n_params, radii, n_radii) &
            result(status) bind(c, name = 'cassini_param_cache_radius_grid')
        type(c_ptr), value, intent(in) :: cache
        integer(ik_c), value, intent(in) :: n_params, n_radii
        real(rk_c), intent(in) :: params(*)
        real(rk_c), intent(out) :: radii(*)
        integer(ik_c) :: status

        type(cache_t), pointer :: p
        integer(kind = ik) :: st, n
        real(kind = rk), allocatable :: params_f(:), radii_f(:)
        logical :: ok_p, ok_w

        n = extent_f(n_radii)
        radii(1:n) = 0.0_rk_c
        if (.not. c_associated(cache)) then
            status = int(SHAPE_ERROR_CACHE_NOT_INITIALIZED, ik_c)
            return
        end if
        call c_f_pointer(cache, p)
        call copy_in_s(params, n_params, params_f, ok_p)
        call work_buffers_s(n_radii, ok_w, radii_f)
        if (.not. (ok_p .and. ok_w)) then
            status = int(SHAPE_ERROR_INVALID_GRID, ik_c)
            return
        end if
        call cache_radius_grid_s(p, params_f, radii_f, st)
        radii(1:n) = real(radii_f, rk_c)
        status = int(st, ik_c)
    end function cassini_param_cache_radius_grid

    function cassini_param_cache_radius_and_derivative(cache, params, n_params, radii, dr_dtheta, &
            n_radii) result(status) bind(c, name = 'cassini_param_cache_radius_and_derivative')
        type(c_ptr), value, intent(in) :: cache
        integer(ik_c), value, intent(in) :: n_params, n_radii
        real(rk_c), intent(in) :: params(*)
        real(rk_c), intent(out) :: radii(*), dr_dtheta(*)
        integer(ik_c) :: status

        type(cache_t), pointer :: p
        integer(kind = ik) :: st, n
        real(kind = rk), allocatable :: params_f(:), radii_f(:), dr_f(:)
        logical :: ok_p, ok_w

        n = extent_f(n_radii)
        radii(1:n) = 0.0_rk_c
        dr_dtheta(1:n) = 0.0_rk_c
        if (.not. c_associated(cache)) then
            status = int(SHAPE_ERROR_CACHE_NOT_INITIALIZED, ik_c)
            return
        end if
        call c_f_pointer(cache, p)
        call copy_in_s(params, n_params, params_f, ok_p)
        call work_buffers_s(n_radii, ok_w, radii_f, dr_f)
        if (.not. (ok_p .and. ok_w)) then
            status = int(SHAPE_ERROR_INVALID_GRID, ik_c)
            return
        end if
        call cache_radius_and_derivative_s(p, params_f, radii_f, dr_f, st)
        radii(1:n) = real(radii_f, rk_c)
        dr_dtheta(1:n) = real(dr_f, rk_c)
        status = int(st, ik_c)
    end function cassini_param_cache_radius_and_derivative

    function cassini_param_cache_radius_and_derivative_at_thetas(cache, params, n_params, thetas, &
            n_thetas, radii, dr_dtheta) result(status) &
            bind(c, name = 'cassini_param_cache_radius_and_derivative_at_thetas')
        type(c_ptr), value, intent(in) :: cache
        integer(ik_c), value, intent(in) :: n_params, n_thetas
        real(rk_c), intent(in) :: params(*), thetas(*)
        real(rk_c), intent(out) :: radii(*), dr_dtheta(*)
        integer(ik_c) :: status

        type(cache_t), pointer :: p
        integer(kind = ik) :: st, n
        real(kind = rk), allocatable :: params_f(:), thetas_f(:), radii_f(:), dr_f(:)
        logical :: ok_p, ok_t, ok_w

        n = extent_f(n_thetas)
        radii(1:n) = 0.0_rk_c
        dr_dtheta(1:n) = 0.0_rk_c
        if (.not. c_associated(cache)) then
            status = int(SHAPE_ERROR_CACHE_NOT_INITIALIZED, ik_c)
            return
        end if
        call c_f_pointer(cache, p)
        call copy_in_s(params, n_params, params_f, ok_p)
        call copy_in_s(thetas, n_thetas, thetas_f, ok_t)
        call work_buffers_s(n_thetas, ok_w, radii_f, dr_f)
        if (.not. (ok_p .and. ok_t .and. ok_w)) then
            status = int(SHAPE_ERROR_INVALID_GRID, ik_c)
            return
        end if
        call cache_radius_and_derivative_at_thetas_s(p, params_f, thetas_f, radii_f, dr_f, st)
        radii(1:n) = real(radii_f, rk_c)
        dr_dtheta(1:n) = real(dr_f, rk_c)
        status = int(st, ik_c)
    end function cassini_param_cache_radius_and_derivative_at_thetas

    function cassini_param_cache_shape(cache, params, n_params, z_shift, r_north, r_south) &
            result(status) bind(c, name = 'cassini_param_cache_shape')
        type(c_ptr), value, intent(in) :: cache
        integer(ik_c), value, intent(in) :: n_params
        real(rk_c), intent(in) :: params(*)
        real(rk_c), intent(out) :: z_shift, r_north, r_south
        integer(ik_c) :: status

        type(cache_t), pointer :: p
        integer(kind = ik) :: st
        real(kind = rk), allocatable :: params_f(:)
        real(kind = rk) :: zs, rn, rs
        logical :: ok

        z_shift = 0.0_rk_c
        r_north = 0.0_rk_c
        r_south = 0.0_rk_c
        if (.not. c_associated(cache)) then
            status = int(SHAPE_ERROR_CACHE_NOT_INITIALIZED, ik_c)
            return
        end if
        call c_f_pointer(cache, p)
        call copy_in_s(params, n_params, params_f, ok)
        if (.not. ok) then
            status = int(SHAPE_ERROR_INVALID_GRID, ik_c)
            return
        end if
        call cache_shape_s(p, params_f, zs, rn, rs, st)
        z_shift = real(zs, rk_c)
        r_north = real(rn, rk_c)
        r_south = real(rs, rk_c)
        status = int(st, ik_c)
    end function cassini_param_cache_shape

    function cassini_param_cache_rho_profile(cache, params, n_params, n_z, z, rho, drho_dz) &
            result(status) bind(c, name = 'cassini_param_cache_rho_profile')
        type(c_ptr), value, intent(in) :: cache
        integer(ik_c), value, intent(in) :: n_params, n_z
        real(rk_c), intent(in) :: params(*)
        real(rk_c), intent(out) :: z(*), rho(*), drho_dz(*)
        integer(ik_c) :: status
        status = cached_profile_f(cache, params, n_params, n_z, .true., z, rho, drho_dz)
    end function cassini_param_cache_rho_profile

    function cassini_param_cache_rho_profile_unchecked(cache, params, n_params, n_z, z, rho, &
            drho_dz) result(status) bind(c, name = 'cassini_param_cache_rho_profile_unchecked')
        type(c_ptr), value, intent(in) :: cache
        integer(ik_c), value, intent(in) :: n_params, n_z
        real(rk_c), intent(in) :: params(*)
        real(rk_c), intent(out) :: z(*), rho(*), drho_dz(*)
        integer(ik_c) :: status
        status = cached_profile_f(cache, params, n_params, n_z, .false., z, rho, drho_dz)
    end function cassini_param_cache_rho_profile_unchecked

    !> Both cached profile forms.
    function cached_profile_f(cache, params, n_params, n_z, checked, z, rho, drho_dz) result(status)
        type(c_ptr), intent(in) :: cache
        integer(ik_c), intent(in) :: n_params, n_z
        real(rk_c), intent(in) :: params(*)
        logical, intent(in) :: checked
        real(rk_c), intent(out) :: z(*), rho(*), drho_dz(*)
        integer(ik_c) :: status

        type(cache_t), pointer :: p
        integer(kind = ik) :: st, n
        real(kind = rk), allocatable :: params_f(:), z_f(:), rho_f(:), drho_f(:)
        logical :: ok_p, ok_w

        n = extent_f(n_z)
        z(1:n) = 0.0_rk_c
        rho(1:n) = 0.0_rk_c
        drho_dz(1:n) = 0.0_rk_c
        if (.not. c_associated(cache)) then
            status = int(SHAPE_ERROR_CACHE_NOT_INITIALIZED, ik_c)
            return
        end if
        call c_f_pointer(cache, p)
        call copy_in_s(params, n_params, params_f, ok_p)
        call work_buffers_s(n_z, ok_w, z_f, rho_f, drho_f)
        if (.not. (ok_p .and. ok_w)) then
            status = int(SHAPE_ERROR_INVALID_GRID, ik_c)
            return
        end if
        if (checked) then
            call cache_rho_profile_s(p, params_f, z_f, rho_f, drho_f, st)
        else
            call cache_rho_profile_unchecked_s(p, params_f, z_f, rho_f, drho_f, st)
        end if
        z(1:n) = real(z_f, rk_c)
        rho(1:n) = real(rho_f, rk_c)
        drho_dz(1:n) = real(drho_f, rk_c)
        status = int(st, ik_c)
    end function cached_profile_f

    function cassini_param_cache_neck(cache, params, n_params, z_neck, rho_neck, found) &
            result(status) bind(c, name = 'cassini_param_cache_neck')
        type(c_ptr), value, intent(in) :: cache
        integer(ik_c), value, intent(in) :: n_params
        real(rk_c), intent(in) :: params(*)
        real(rk_c), intent(out) :: z_neck, rho_neck
        integer(ik_c), intent(out) :: found
        integer(ik_c) :: status

        type(cache_t), pointer :: p
        integer(kind = ik) :: st
        real(kind = rk), allocatable :: params_f(:)
        real(kind = rk) :: zn, rn
        logical :: ok, fnd

        z_neck = 0.0_rk_c
        rho_neck = 0.0_rk_c
        found = 0_ik_c
        if (.not. c_associated(cache)) then
            status = int(SHAPE_ERROR_CACHE_NOT_INITIALIZED, ik_c)
            return
        end if
        call c_f_pointer(cache, p)
        call copy_in_s(params, n_params, params_f, ok)
        if (.not. ok) then
            status = int(SHAPE_ERROR_INVALID_GRID, ik_c)
            return
        end if
        call cache_neck_s(p, params_f, zn, rn, fnd, st)
        z_neck = real(zn, rk_c)
        rho_neck = real(rn, rk_c)
        found = merge(1_ik_c, 0_ik_c, fnd)
        status = int(st, ik_c)
    end function cassini_param_cache_neck

    function cassini_param_cache_geometry(cache, params, n_params, geometry) result(status) &
            bind(c, name = 'cassini_param_cache_geometry')
        type(c_ptr), value, intent(in) :: cache
        integer(ik_c), value, intent(in) :: n_params
        real(rk_c), intent(in) :: params(*)
        type(cassini_param_geometry_t), intent(out) :: geometry
        integer(ik_c) :: status

        type(cache_t), pointer :: p
        type(cassini_geometry_t) :: g
        integer(kind = ik) :: st
        real(kind = rk), allocatable :: params_f(:)
        logical :: ok

        call zero_geometry_c(geometry)
        if (.not. c_associated(cache)) then
            status = int(SHAPE_ERROR_CACHE_NOT_INITIALIZED, ik_c)
            return
        end if
        call c_f_pointer(cache, p)
        call copy_in_s(params, n_params, params_f, ok)
        if (.not. ok) then
            status = int(SHAPE_ERROR_INVALID_GRID, ik_c)
            return
        end if
        call cache_geometry_s(p, params_f, g, st)
        call geometry_to_c(g, geometry)
        status = int(st, ik_c)
    end function cassini_param_cache_geometry

    function cassini_param_cache_quadrupole(cache, params, n_params, q2) result(status) &
            bind(c, name = 'cassini_param_cache_quadrupole')
        type(c_ptr), value, intent(in) :: cache
        integer(ik_c), value, intent(in) :: n_params
        real(rk_c), intent(in) :: params(*)
        real(rk_c), intent(out) :: q2
        integer(ik_c) :: status

        type(cache_t), pointer :: p
        integer(kind = ik) :: st
        real(kind = rk), allocatable :: params_f(:)
        real(kind = rk) :: q
        logical :: ok

        q2 = 0.0_rk_c
        if (.not. c_associated(cache)) then
            status = int(SHAPE_ERROR_CACHE_NOT_INITIALIZED, ik_c)
            return
        end if
        call c_f_pointer(cache, p)
        call copy_in_s(params, n_params, params_f, ok)
        if (.not. ok) then
            status = int(SHAPE_ERROR_INVALID_GRID, ik_c)
            return
        end if
        call cache_quadrupole_s(p, params_f, q, st)
        q2 = real(q, rk_c)
        status = int(st, ik_c)
    end function cassini_param_cache_quadrupole

    !===========================================================================
    ! ONE-SHOT COMPUTES
    !===========================================================================
    ! Each entry zero-fills its outputs over the stated extent, rejects a
    ! negative size argument (n_params, n_thetas, n_z) as 5, copies in (a heap
    ! refusal is 5) and hands the tier the copies; n_x_grid passes through.

    subroutine cassini_param_radius_grid(params, n_params, n_x_grid, thetas, n_thetas, radii, &
            status) bind(c, name = 'cassini_param_radius_grid')
        integer(ik_c), value, intent(in) :: n_params, n_x_grid, n_thetas
        real(rk_c), intent(in) :: params(*), thetas(*)
        real(rk_c), intent(out) :: radii(*)
        integer(ik_c), intent(out), optional :: status

        integer(kind = ik) :: st, n
        real(kind = rk), allocatable :: params_f(:), thetas_f(:), radii_f(:)
        logical :: ok_p, ok_t, ok_w

        n = extent_f(n_thetas)
        radii(1:n) = 0.0_rk_c
        if (n_params < 0_ik_c .or. n_thetas < 0_ik_c) then
            if (present(status)) status = int(SHAPE_ERROR_INVALID_INIT, ik_c)
            return
        end if
        call copy_in_s(params, n_params, params_f, ok_p)
        call copy_in_s(thetas, n_thetas, thetas_f, ok_t)
        call work_buffers_s(n_thetas, ok_w, radii_f)
        if (.not. (ok_p .and. ok_t .and. ok_w)) then
            if (present(status)) status = int(SHAPE_ERROR_INVALID_INIT, ik_c)
            return
        end if
        call compute_radius_grid_standalone_s(params_f, int(n_x_grid, ik), thetas_f, radii_f, st)
        radii(1:n) = real(radii_f, rk_c)
        if (present(status)) status = int(st, ik_c)
    end subroutine cassini_param_radius_grid

    subroutine cassini_param_radius_and_derivative(params, n_params, n_x_grid, thetas, n_thetas, &
            radii, dr_dtheta, status) bind(c, name = 'cassini_param_radius_and_derivative')
        integer(ik_c), value, intent(in) :: n_params, n_x_grid, n_thetas
        real(rk_c), intent(in) :: params(*), thetas(*)
        real(rk_c), intent(out) :: radii(*), dr_dtheta(*)
        integer(ik_c), intent(out), optional :: status

        integer(kind = ik) :: st, n
        real(kind = rk), allocatable :: params_f(:), thetas_f(:), radii_f(:), dr_f(:)
        logical :: ok_p, ok_t, ok_w

        n = extent_f(n_thetas)
        radii(1:n) = 0.0_rk_c
        dr_dtheta(1:n) = 0.0_rk_c
        if (n_params < 0_ik_c .or. n_thetas < 0_ik_c) then
            if (present(status)) status = int(SHAPE_ERROR_INVALID_INIT, ik_c)
            return
        end if
        call copy_in_s(params, n_params, params_f, ok_p)
        call copy_in_s(thetas, n_thetas, thetas_f, ok_t)
        call work_buffers_s(n_thetas, ok_w, radii_f, dr_f)
        if (.not. (ok_p .and. ok_t .and. ok_w)) then
            if (present(status)) status = int(SHAPE_ERROR_INVALID_INIT, ik_c)
            return
        end if
        call compute_radius_and_derivative_standalone_s(params_f, int(n_x_grid, ik), thetas_f, &
                radii_f, dr_f, st)
        radii(1:n) = real(radii_f, rk_c)
        dr_dtheta(1:n) = real(dr_f, rk_c)
        if (present(status)) status = int(st, ik_c)
    end subroutine cassini_param_radius_and_derivative

    subroutine cassini_param_shape(params, n_params, n_x_grid, z_shift, r_north, r_south, status) &
            bind(c, name = 'cassini_param_shape')
        integer(ik_c), value, intent(in) :: n_params, n_x_grid
        real(rk_c), intent(in) :: params(*)
        real(rk_c), intent(out) :: z_shift, r_north, r_south
        integer(ik_c), intent(out), optional :: status

        integer(kind = ik) :: st
        real(kind = rk), allocatable :: params_f(:)
        real(kind = rk) :: zs, rn, rs
        logical :: ok

        z_shift = 0.0_rk_c
        r_north = 0.0_rk_c
        r_south = 0.0_rk_c
        if (n_params < 0_ik_c) then
            if (present(status)) status = int(SHAPE_ERROR_INVALID_INIT, ik_c)
            return
        end if
        call copy_in_s(params, n_params, params_f, ok)
        if (.not. ok) then
            if (present(status)) status = int(SHAPE_ERROR_INVALID_INIT, ik_c)
            return
        end if
        call compute_shape_standalone_s(params_f, int(n_x_grid, ik), zs, rn, rs, st)
        z_shift = real(zs, rk_c)
        r_north = real(rn, rk_c)
        r_south = real(rs, rk_c)
        if (present(status)) status = int(st, ik_c)
    end subroutine cassini_param_shape

    subroutine cassini_param_rho_profile(params, n_params, n_z, z, rho, drho_dz, status) &
            bind(c, name = 'cassini_param_rho_profile')
        integer(ik_c), value, intent(in) :: n_params, n_z
        real(rk_c), intent(in) :: params(*)
        real(rk_c), intent(out) :: z(*), rho(*), drho_dz(*)
        integer(ik_c), intent(out), optional :: status
        integer(ik_c) :: st
        call one_shot_profile_s(params, n_params, n_z, .true., z, rho, drho_dz, st)
        if (present(status)) status = st
    end subroutine cassini_param_rho_profile

    subroutine cassini_param_rho_profile_unchecked(params, n_params, n_z, z, rho, drho_dz, status) &
            bind(c, name = 'cassini_param_rho_profile_unchecked')
        integer(ik_c), value, intent(in) :: n_params, n_z
        real(rk_c), intent(in) :: params(*)
        real(rk_c), intent(out) :: z(*), rho(*), drho_dz(*)
        integer(ik_c), intent(out), optional :: status
        integer(ik_c) :: st
        call one_shot_profile_s(params, n_params, n_z, .false., z, rho, drho_dz, st)
        if (present(status)) status = st
    end subroutine cassini_param_rho_profile_unchecked

    !> Both one-shot profile forms.
    subroutine one_shot_profile_s(params, n_params, n_z, checked, z, rho, drho_dz, status)
        integer(ik_c), intent(in) :: n_params, n_z
        real(rk_c), intent(in) :: params(*)
        logical, intent(in) :: checked
        real(rk_c), intent(out) :: z(*), rho(*), drho_dz(*)
        integer(ik_c), intent(out) :: status

        integer(kind = ik) :: st, n
        real(kind = rk), allocatable :: params_f(:), z_f(:), rho_f(:), drho_f(:)
        logical :: ok_p, ok_w

        n = extent_f(n_z)
        z(1:n) = 0.0_rk_c
        rho(1:n) = 0.0_rk_c
        drho_dz(1:n) = 0.0_rk_c
        if (n_params < 0_ik_c .or. n_z < 0_ik_c) then
            status = int(SHAPE_ERROR_INVALID_INIT, ik_c)
            return
        end if
        call copy_in_s(params, n_params, params_f, ok_p)
        call work_buffers_s(n_z, ok_w, z_f, rho_f, drho_f)
        if (.not. (ok_p .and. ok_w)) then
            status = int(SHAPE_ERROR_INVALID_INIT, ik_c)
            return
        end if
        if (checked) then
            call compute_rho_profile_standalone_s(params_f, z_f, rho_f, drho_f, st)
        else
            call compute_rho_profile_unchecked_standalone_s(params_f, z_f, rho_f, drho_f, st)
        end if
        z(1:n) = real(z_f, rk_c)
        rho(1:n) = real(rho_f, rk_c)
        drho_dz(1:n) = real(drho_f, rk_c)
        status = int(st, ik_c)
    end subroutine one_shot_profile_s

    subroutine cassini_param_neck(params, n_params, n_x_grid, z_neck, rho_neck, found, status) &
            bind(c, name = 'cassini_param_neck')
        integer(ik_c), value, intent(in) :: n_params, n_x_grid
        real(rk_c), intent(in) :: params(*)
        real(rk_c), intent(out) :: z_neck, rho_neck
        integer(ik_c), intent(out) :: found
        integer(ik_c), intent(out), optional :: status

        integer(kind = ik) :: st
        real(kind = rk), allocatable :: params_f(:)
        real(kind = rk) :: zn, rn
        logical :: ok, fnd

        z_neck = 0.0_rk_c
        rho_neck = 0.0_rk_c
        found = 0_ik_c
        if (n_params < 0_ik_c) then
            if (present(status)) status = int(SHAPE_ERROR_INVALID_INIT, ik_c)
            return
        end if
        call copy_in_s(params, n_params, params_f, ok)
        if (.not. ok) then
            if (present(status)) status = int(SHAPE_ERROR_INVALID_INIT, ik_c)
            return
        end if
        call compute_neck_standalone_s(params_f, int(n_x_grid, ik), zn, rn, fnd, st)
        z_neck = real(zn, rk_c)
        rho_neck = real(rn, rk_c)
        found = merge(1_ik_c, 0_ik_c, fnd)
        if (present(status)) status = int(st, ik_c)
    end subroutine cassini_param_neck

    subroutine cassini_param_geometry(params, n_params, n_x_grid, geometry, status) &
            bind(c, name = 'cassini_param_geometry')
        integer(ik_c), value, intent(in) :: n_params, n_x_grid
        real(rk_c), intent(in) :: params(*)
        type(cassini_param_geometry_t), intent(out) :: geometry
        integer(ik_c), intent(out), optional :: status

        type(cassini_geometry_t) :: g
        integer(kind = ik) :: st
        real(kind = rk), allocatable :: params_f(:)
        logical :: ok

        call zero_geometry_c(geometry)
        if (n_params < 0_ik_c) then
            if (present(status)) status = int(SHAPE_ERROR_INVALID_INIT, ik_c)
            return
        end if
        call copy_in_s(params, n_params, params_f, ok)
        if (.not. ok) then
            if (present(status)) status = int(SHAPE_ERROR_INVALID_INIT, ik_c)
            return
        end if
        call compute_geometry_standalone_s(params_f, int(n_x_grid, ik), g, st)
        call geometry_to_c(g, geometry)
        if (present(status)) status = int(st, ik_c)
    end subroutine cassini_param_geometry

    subroutine cassini_param_quadrupole(params, n_params, q2, status) &
            bind(c, name = 'cassini_param_quadrupole')
        integer(ik_c), value, intent(in) :: n_params
        real(rk_c), intent(in) :: params(*)
        real(rk_c), intent(out) :: q2
        integer(ik_c), intent(out), optional :: status

        integer(kind = ik) :: st
        real(kind = rk), allocatable :: params_f(:)
        real(kind = rk) :: q
        logical :: ok

        q2 = 0.0_rk_c
        if (n_params < 0_ik_c) then
            if (present(status)) status = int(SHAPE_ERROR_INVALID_INIT, ik_c)
            return
        end if
        call copy_in_s(params, n_params, params_f, ok)
        if (.not. ok) then
            if (present(status)) status = int(SHAPE_ERROR_INVALID_INIT, ik_c)
            return
        end if
        call compute_quadrupole_standalone_s(params_f, q, st)
        q2 = real(q, rk_c)
        if (present(status)) status = int(st, ik_c)
    end subroutine cassini_param_quadrupole

end module cassini_parameterization_c_api_mod
