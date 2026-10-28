!> Public surface of the Cassini-oval nuclear shape library.
!!
!! The library's front door: status codes and messages, the read-only cache,
!! the usage checks, trailing-zero trimming, the gates and both tiers. It does
!! no floating-point arithmetic on an output value; every kernel lives in
!! `cassini_parameterization_workers_mod`.
module cassini_parameterization_mod

    use precision_utilities_mod, only: ik
    use shape_core_mod, only: SHAPE_VALID, SHAPE_ERROR_TOO_MANY_PARAMS, &
            SHAPE_ERROR_CACHE_NOT_INITIALIZED, SHAPE_ERROR_INVALID_GRID, &
            SHAPE_ERROR_WRONG_PARAM_COUNT, SHAPE_ERROR_INVALID_INIT, &
            SHAPE_MAX_PARAMS
    use cassini_parameterization_workers_mod, only: CASSINI_MAX_PARAMS, CASSINI_MAX_ORDER, &
            CASSINI_N_GL_HALF, CASSINI_N_R_CHECK, CASSINI_N_X_GRID_FLOOR, &
            CASSINI_N_POINTS_FLOOR, STAR_CONVEXITY_MARGIN, NECK_TIE_TOL, X_SOLVE_TOL, &
            STATUS_MESSAGE_LEN

    implicit none

    private

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

    integer(kind = ik), parameter, public :: CASSINI_VALID = SHAPE_VALID
    !> R(x) <= 0 at a node of the fixed check grid.
    integer(kind = ik), parameter, public :: CASSINI_ERROR_R_NONPOSITIVE = 100_ik
    !> R(1)^2 + s <= 0 or R(-1)^2 + s <= 0: a pole has left the axis.
    integer(kind = ik), parameter, public :: CASSINI_ERROR_TOROIDAL = 101_ik
    !> R(theta) outputs and shape only: pinched, or d_min below the margin.
    integer(kind = ik), parameter, public :: CASSINI_ERROR_NOT_STAR_CONVEX = 102_ik
    !> R(theta) outputs, shape and the checked profile: two fragments.
    integer(kind = ik), parameter, public :: CASSINI_ERROR_SEPARATED = 103_ik

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

end module cassini_parameterization_mod
