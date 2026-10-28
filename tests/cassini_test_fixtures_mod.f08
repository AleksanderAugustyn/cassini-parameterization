!> Named fixtures of the cassini test suites. Values are [alpha, alpha1, ...,
!! alphaN] in their natural lengths (the library pads); the classification
!! and cross-check numbers live in the plan's fixture table.
module cassini_test_fixtures_mod

    use precision_utilities_mod, only: rk

    implicit none

    private

    public :: SPHERE, BARE_05, BARE_095, OBLATE, GS, SADDLE1, SD_MIN, SADDLE2, ASYM, &
            DIMPLED, FOLD_BOX, TWO_NECKS, STAR_FAIL, PINCHED, PINCHED_ASYM, SEPARATED, &
            NEAR_PINCH, R_NEG, TOROIDAL, TOROIDAL_ONE, THETA_ANCHORS

    real(kind = rk), parameter :: SPHERE(1) = [0.0_rk]
    real(kind = rk), parameter :: BARE_05(1) = [0.5_rk]
    real(kind = rk), parameter :: BARE_095(1) = [0.95_rk]
    real(kind = rk), parameter :: OBLATE(1) = [-0.25_rk]
    !> 236U ground state, first saddle, secondary minimum, second saddle (J fig. 6).
    real(kind = rk), parameter :: GS(5) = [0.25_rk, 0.0_rk, 0.0_rk, 0.0_rk, 0.075_rk]
    real(kind = rk), parameter :: SADDLE1(5) = [0.35_rk, 0.0_rk, 0.0_rk, 0.0_rk, -0.075_rk]
    real(kind = rk), parameter :: SD_MIN(5) = [0.525_rk, 0.0_rk, 0.0_rk, 0.0_rk, 0.025_rk]
    real(kind = rk), parameter :: SADDLE2(5) = [0.65_rk, 0.2_rk, 0.0_rk, 0.025_rk, 0.05_rk]
    !> Necked, asymmetric, folded; R'(0) = alpha1 - 1.5 alpha3 = 0 puts the neck at x = 0.
    real(kind = rk), parameter :: ASYM(5) = [0.8_rk, 0.15_rk, -0.05_rk, 0.1_rk, -0.05_rk]
    real(kind = rk), parameter :: DIMPLED(5) = [0.5_rk, 0.0_rk, 0.0_rk, 0.0_rk, -0.2_rk]
    real(kind = rk), parameter :: FOLD_BOX(5) = [0.125_rk, 0.25_rk, 0.0_rk, 0.2_rk, -0.175_rk]
    !> Symmetric with a central maximum of rho and two necks near x = +-0.846.
    real(kind = rk), parameter :: TWO_NECKS(5) = [0.25_rk, 0.0_rk, -0.4_rk, 0.0_rk, 0.5_rk]
    !> Connected, not folded, d_min < 0.
    real(kind = rk), parameter :: STAR_FAIL(5) = [0.95_rk, -0.2_rk, 0.0_rk, -0.05_rk, -0.05_rk]
    real(kind = rk), parameter :: PINCHED(1) = [1.0_rk]
    real(kind = rk), parameter :: PINCHED_ASYM(5) = [1.0_rk, 0.2_rk, 0.0_rk, 0.025_rk, 0.05_rk]
    real(kind = rk), parameter :: SEPARATED(1) = [1.05_rk]
    real(kind = rk), parameter :: NEAR_PINCH(1) = [1.0_rk - 1.0e-12_rk]
    !> R(1) = -0.1: code 100.
    real(kind = rk), parameter :: R_NEG(2) = [0.0_rk, -1.1_rk]
    !> R(+-1)^2 + s = -0.1: code 101.
    real(kind = rk), parameter :: TOROIDAL(1) = [-1.1_rk]
    !> R(-1)^2 + s < 0 < R(1)^2 + s: code 101, one-sided.
    real(kind = rk), parameter :: TOROIDAL_ONE(2) = [-0.9_rk, 0.3_rk]

    !> pi/8, pi/2, 7 pi/8, pre-rounded: the golden and cross-check thetas.
    real(kind = rk), parameter :: THETA_ANCHORS(3) = [0.39269908169872414_rk, &
            1.5707963267948966_rk, 2.748893571891069_rk]

end module cassini_test_fixtures_mod
