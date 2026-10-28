!> Prints golden values for the regression shapes as paste-ready Python.
!!
!! ## Recapture procedure
!!
!! The goldens in `python/tests/test_api.py` are comparisons against this
!! library at GOLDEN_ATOL = 1e-13. They must be regenerated whenever the
!! arithmetic behind the C API changes deliberately -- a new quadrature, a new
!! derivative form, a new root-find -- and NOT otherwise: a golden that drifts
!! without a deliberate recapture is a regression, not a refresh.
!!
!!   1. Build the capture program (Debug):
!!        `cmake --build build/debug --target cassini_param_golden_capture`
!!   2. Run it and keep the output:
!!        `./build/debug/cassini_param_golden_capture`
!!      Every line but the closing `Tests:` summary is valid Python. A
!!      non-zero exit or a FAIL line means a golden shape stopped being
!!      accepted -- investigate, do not paste.
!!   3. Replace the generated block of `python/tests/test_api.py` (between the
!!      `--- golden block` markers) with the printed lines.
!!   4. Run the Python suite against THIS tree, not an installed package:
!!        `CASSINI_PARAM_LIB=$PWD/build/release/libcassini_parameterization.so \
!!           PYTHONPATH=$PWD/python python -m pytest python/tests -q -p no:cacheprovider`
!!      Both variables are required: a bare `pytest` imports whatever
!!      `cassini_parameterization` is installed and loads whatever library
!!      that resolves. The suite must pass with GOLDEN_ATOL UNCHANGED;
!!      loosening it to make a recapture fit is never the fix.
!!
!! ## The tolerance budget
!!
!! Capture is Debug, but the goldens are compared against local Release
!! (-O3 -ffast-math) and the manylinux wheel (GCC 10.2) as well. The resolve
!! is a fixed-order quadrature and the root-finds stop at X_SOLVE_TOL = 1e-14
!! in x, so every value inherits the FMA/contraction differences between
!! those builds: a few ulp of values of order unity. GOLDEN_ATOL = 1e-13
!! leaves two orders of magnitude for that gap; an algorithm change that eats
!! it must be measured against the shipped configuration, never answered by a
!! looser tolerance.
!!
!! Every value is captured through the SAME C entries the Python bindings
!! call, with n_x_grid = 129 (the Python tests build Cache(5, 129, thetas)),
!! so the goldens include the C-level marshalling.
program golden_capture

    use precision_utilities_mod, only: ik, rk
    use c_bindings_mod, only: ik_c, rk_c
    use cassini_parameterization_c_api_mod, only: cassini_param_geometry_t, cassini_param_geometry, &
            cassini_param_neck, cassini_param_quadrupole, cassini_param_radius_and_derivative, &
            cassini_param_rho_profile
    use cassini_parameterization_mod, only: CASSINI_VALID
    use test_utils_mod, only: assert_int_eq, test_summary

    implicit none

    integer(kind = ik), parameter :: N_Z = 181_ik
    integer(kind = ik), parameter :: IDX(5) = [1_ik, 46_ik, 91_ik, 136_ik, 181_ik]
    integer(kind = ik_c), parameter :: N_X_GRID = 129_ik_c

    call capture('GS', [0.25_rk, 0.0_rk, 0.0_rk, 0.0_rk, 0.075_rk])
    call capture('SADDLE1', [0.35_rk, 0.0_rk, 0.0_rk, 0.0_rk, -0.075_rk])
    call capture('SD_MIN', [0.525_rk, 0.0_rk, 0.0_rk, 0.0_rk, 0.025_rk])
    call capture('SADDLE2', [0.65_rk, 0.2_rk, 0.0_rk, 0.025_rk, 0.05_rk])
    call capture('ASYM', [0.8_rk, 0.15_rk, -0.05_rk, 0.1_rk, -0.05_rk])
    call capture('DIMPLED', [0.5_rk, 0.0_rk, 0.0_rk, 0.0_rk, -0.2_rk])
    call test_summary()

contains

    subroutine print_list(name, suffix, values)
        character(len = *), intent(in) :: name, suffix
        real(kind = rk_c), intent(in) :: values(:)
        integer(kind = ik) :: i
        write(*, '(A,A,A)') name, suffix, ' = ['
        do i = 1_ik, size(values, kind = ik)
            write(*, '(A,ES24.16E3,A)') '    ', values(i), ','
        end do
        write(*, '(A)') ']'
    end subroutine print_list

    subroutine capture(name, params)
        character(len = *), intent(in) :: name
        real(kind = rk), intent(in) :: params(5)

        real(kind = rk_c) :: p(5), geo(7), nk(2), q2(1), z(N_Z), rho(N_Z), drho(N_Z)
        real(kind = rk_c) :: eval_thetas(3), eval_r(3), eval_dr(3), picked(size(IDX))
        type(cassini_param_geometry_t) :: g
        integer(kind = ik_c) :: code, found

        p = real(params, rk_c)

        call cassini_param_geometry(p, 5_ik_c, N_X_GRID, g, code)
        call assert_int_eq(int(code, ik), CASSINI_VALID, name // ': geometry valid')
        if (int(code, ik) /= CASSINI_VALID) return
        geo(1) = g%eps
        geo(2) = g%c
        geo(3) = g%z_shift
        geo(4) = g%r_north
        geo(5) = g%r_south
        geo(6) = g%gap
        geo(7) = g%d_min
        call print_list(name, '_GEOMETRY', geo)
        write(*, '(8A)') name, '_FLAGS = (', trim(merge('True ', 'False', g%is_separated == 1_ik_c)), ', ', &
                trim(merge('True ', 'False', g%is_pinched == 1_ik_c)), ', ', &
                trim(merge('True ', 'False', g%is_folded == 1_ik_c)), ')'

        call cassini_param_neck(p, 5_ik_c, N_X_GRID, nk(1), nk(2), found, code)
        call assert_int_eq(int(code, ik), CASSINI_VALID, name // ': neck valid')
        write(*, '(A,A,A)') name, '_NECK_FOUND = ', trim(merge('True ', 'False', found == 1_ik_c))
        call print_list(name, '_NECK', nk)

        call cassini_param_quadrupole(p, 5_ik_c, q2(1), code)
        call assert_int_eq(int(code, ik), CASSINI_VALID, name // ': quadrupole valid')
        call print_list(name, '_Q2', q2)

        ! Thetas as pre-rounded literals (pi/8, pi/2, 7pi/8), NOT arithmetic:
        ! Release builds may reassociate 7*PI_C/8 during constant folding;
        ! literals parse identically here and in Python.
        eval_thetas = [0.39269908169872414_rk_c, 1.5707963267948966_rk_c, &
                2.748893571891069_rk_c]
        call cassini_param_radius_and_derivative(p, 5_ik_c, N_X_GRID, eval_thetas, 3_ik_c, eval_r, &
                eval_dr, code)
        call assert_int_eq(int(code, ik), CASSINI_VALID, name // ': radius_and_derivative valid')
        call print_list(name, '_DERIV_R', eval_r)
        call print_list(name, '_DERIV_DR', eval_dr)

        call cassini_param_rho_profile(p, 5_ik_c, int(N_Z, ik_c), z, rho, drho, code)
        call assert_int_eq(int(code, ik), CASSINI_VALID, name // ': rho_profile valid')
        picked = z(IDX)
        call print_list(name, '_Z', picked)
        picked = rho(IDX)
        call print_list(name, '_RHO', picked)
    end subroutine capture

end program golden_capture
