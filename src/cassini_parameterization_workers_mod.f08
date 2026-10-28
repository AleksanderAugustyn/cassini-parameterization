!> Arithmetic of the Cassini-oval nuclear shape library.
!!
!! Owns the constants of the library and every piece of arithmetic: the
!! Legendre tables and the series R(x), eps from alpha, the lemniscate map
!! with its conditioning rules and derivatives, the node sets, the
!! Gauss-Legendre resolve in t (volume, centre of mass, quadrupole), the scan
!! with the star margin, the neck, the profile and the R(theta) root-find.
module cassini_parameterization_workers_mod

    use precision_utilities_mod, only: ik, rk

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

end module cassini_parameterization_workers_mod
